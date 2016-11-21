##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: CancerIT <cgpit@sanger.ac.uk>
#
# This file is part of AscatNGS.
#
# AscatNGS is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

package Sanger::CGP::Ascat::Implement;

use strict;

use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);
use Cwd qw(abs_path getcwd);
use File::Basename;
use File::Spec;
use File::Which qw(which);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile);
use File::Copy qw(copy move);
use Capture::Tiny qw(capture);
use FindBin qw($Bin);
use List::Util qw(first);

use File::ShareDir qw(module_dir);

use Sanger::CGP::Ascat;

use PCAP::Threaded;
use PCAP::Bam;

const my $COUNT_READS => q{%s view -c %s %s};

const my $FAILED_SAMPLE_STATISTICS => qq{## WARNING ASCAT failed to %s ##\nNormalContamination 0.3\nPloidy ?\nrho 0\npsi 0\ngoodnessOfFit 0\n};

const my $ALLELE_COUNT_GENDER => ' -b %s -l %s -r %s -g ';

const my $ALLELE_COUNT_PARA => ' -b %s -o %s -l %s -c %s -r %s ';

const my $GREP_ALLELE_COUNTS => q{grep -v '^#' %s >> %s};

const my @ASCAT_RESULT_FILES => qw(
                                    %s.ASPCF.png
                                    %s.germline.png
                                    %s.tumour.png
                                    %s.copynumber.caveman.csv
                                    %s.copynumber.txt
                                    %s.samplestatistics.txt
                                  );
const my @ASCAT_OPTIONAL_PNG => qw(
                                    %s.ASCATprofile.png
                                    %s.rawprofile.png
                                    %s.sunrise.png
                                  );

const my $GENDER_MIN => 5;

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;
  return $self;
}

sub allele_count {
 my ($index_in, $options) = @_;

  my $tmp = $options->{'tmp'};

	# first handle the easy bit, skip if limit not set
	return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});

  my @seqs = snpLociChrs($options);
  my $chrCount = scalar @seqs;
  my @indicies = limited_indices($options, $index_in, $options->{'lociChrsBySample'});

  my $tmp_cover = File::Spec->catdir($tmp, 'allele_count');
  make_path($tmp_cover) unless(-e $tmp_cover);

  for my $index(@indicies) {
    my $samp_type = 'tumour';
    my $seq_idx = $index-1;
    if($index > $chrCount) {
      $samp_type = 'normal';
      $seq_idx = ($index-$chrCount)-1;
    }
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

    my $ac_out = File::Spec->catdir($tmp, 'allele_count');
    make_path($ac_out) unless(-e $ac_out);
    my $loci_files = File::Spec->catdir($tmp, 'allele_count', 'loci_files');
    make_path($loci_files) unless(-e $loci_files);

    my $chr = $seqs[$seq_idx-1];
    my $sname = sanitised_sample_from_bam($options->{$samp_type});
    my $alleleCountOut = File::Spec->catfile($ac_out,sprintf '%s.%s.allct', $sname, $chr);

    # first we need a loci file, generate from the gc file:
    my $loci_file = "$loci_files/$sname.$chr";
    my $gc_to_loci = qq{cut -f 2,3 $options->{snp_gc} | grep -vP '^Chr\\tPosition' | grep -P '^$chr\t' > $loci_file};

    my $allc_exe = _which('alleleCounter');
    my $allc_lib = dirname($allc_exe);

    my $command = $allc_exe;
    $command .= sprintf $ALLELE_COUNT_PARA, $options->{$samp_type}, $alleleCountOut, $loci_file, $chr, $options->{'reference'};
    $command .= '-m '.$options->{'minbasequal'} if exists $options->{'minbasequal'};

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), [$gc_to_loci, $command], $index);

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
}

sub ascat {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  $tmp = abs_path($tmp);
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tum_name = sanitised_sample_from_bam($options->{'tumour'});
  my $norm_name = sanitised_sample_from_bam($options->{'normal'});

  my $ascat_out = File::Spec->catdir($tmp, 'ascat');
  make_path($ascat_out) unless(-e $ascat_out);

  my $rdata = File::Spec->catfile($ascat_out,$tum_name.'.Rdata');

  my $tumcountfile = $tum_name . '.count';
  my $normcountfile = $norm_name . '.count';

  my $tumcount = File::Spec->catfile($ascat_out,$tumcountfile);
  my $normcount = File::Spec->catfile($ascat_out,$normcountfile);

  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_counts_mt', 0)) {
    merge_counts($options, $tmp, $tum_name, $tumcount);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_counts_mt', 0);
  }

  unless(PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_counts_wt', 0)) {
    merge_counts($options, $tmp, $norm_name, $normcount);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_counts_wt', 0);
  }

  my $clean_snp_gc = File::Spec->catfile($tmp,'SnpGcCorrections.tsv');
  # need to strip leading chr so ascat can work
  clean_snp_gc($options->{'snp_gc'}, $clean_snp_gc);

  my $snp_pos = File::Spec->catfile($tmp,'SnpPositions.tsv');
  my $gc_to_snppos = qq{cut -f 1-3 $clean_snp_gc > $snp_pos};

  my @chr_set = snpLociChrs($options);
  my $core_chrs = @chr_set;
  $core_chrs++ unless(first {$_ eq 'Y'  || $_ eq 'chrY'} @chr_set); # incase Y is not in the SNP set

  my $command = "cd $ascat_out; "._which('Rscript');

  my $mod_path = dirname(abs_path($0)).'/../share';
  $mod_path = module_dir('Sanger::CGP::Ascat::Implement') unless(-e File::Spec->catdir($mod_path, 'ascat'));

  my $ascat_path = File::Spec->catdir($mod_path, 'ascat');
  my $ascat_exe = File::Spec->catfile($ascat_path,'runASCAT.R');

  $command .= " $ascat_exe";
  $command .= " $ascat_path";
  $command .= ' '.$snp_pos;
  $command .= ' '.$clean_snp_gc;
  $command .= ' '.$options->{'tumour_name'};
  $command .= ' '.$tumcountfile;
  $command .= ' '.$options->{'normal_name'};
  $command .= ' '.$normcountfile;
  $command .= ' '.$options->{'gender'};
  $command .= ' '.$core_chrs;
  $command .= ' '.$rdata;

  if(defined($options->{'ploidy'}) && defined($options->{'purity'})){
    $command .= ' '.$options->{'purity'};
    $command .= ' '.$options->{'ploidy'};
  }

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), [$gc_to_snppos, $command], 0);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

}

sub finalise {
  my $options = shift;

  my $outdir = $options->{'outdir'};
  my $tmp = $options->{'tmp'};

  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $tum_name = sanitised_sample_from_bam($options->{'tumour'});
  my $ascat_out = File::Spec->catdir($tmp, 'ascat');
  my $cave_cn;
  my $force_complete = 0;
  my @commands;
  foreach my $f(@ASCAT_RESULT_FILES, @ASCAT_OPTIONAL_PNG){
    my $file = sprintf($f, $tum_name);
    my $from = File::Spec->catfile($ascat_out,$file);
    if(!-e $from) {
      if(exists $options->{'force'} && defined $options->{'force'}) {
        $force_complete = 1;
      }
      else {
        die "Expected ASCAT output file missing: $from\n";
      }
    }
    else {
      my $to = File::Spec->catfile($options->{'outdir'},$file);
      copy $from,$to;
      $cave_cn = $to if($to =~ m/copynumber.caveman.csv$/);
    }
  }

  my $cn_txt = File::Spec->catfile($options->{'outdir'}, sprintf('%s.copynumber.txt', $tum_name));
  my $bw_stub = File::Spec->catfile($options->{'outdir'}, sprintf('%s.copynumber', $tum_name));
  my $sunrise = File::Spec->catfile($options->{'outdir'}, sprintf('%s.sunrise.png', $tum_name));

  if($force_complete == 1) {

    my $fail_reason = 'generate a solution';
    $fail_reason = 'determine any purity/ploidy profile' unless (-e $sunrise);

    my $fake_file = sprintf '%s/%s.copynumber.caveman.csv', $options->{'outdir'}, $tum_name;
    my $fake_csv = "$^X ";
    $fake_csv .= _which('ascatFailedCnCsv.pl');
    $fake_csv .= sprintf ' -r %s -o %s', $options->{'reference'}, $fake_file;
    push @commands, $fake_csv;
    my $samp_stat_file = sprintf '%s/%s.samplestatistics.txt', $options->{'outdir'}, $tum_name;
    open my $STAT, '>', $samp_stat_file;
    printf $STAT $FAILED_SAMPLE_STATISTICS, $fail_reason or die "Failed to write line to $samp_stat_file";
    print $STAT "GenderChr $options->{genderChr}\n";
    print $STAT "GenderChrFound $options->{genderIsMale}\n";
    close $STAT;
    $cave_cn = $fake_file;

    my $share_path = dirname(abs_path($0)).'/../share';
    $share_path = module_dir('Sanger::CGP::Ascat::Implement') unless(-e File::Spec->catdir($share_path, 'images'));
    my $img_path = File::Spec->catdir($share_path, 'images');

    for my $f(@ASCAT_OPTIONAL_PNG) {
      my $file = sprintf($f, $tum_name);
      my $to = File::Spec->catfile($options->{'outdir'},$file);
      next if(-e $to);
      if($f eq '%s.sunrise.png') {
        copy(File::Spec->catfile($img_path, 'NoSunrise.png'), $to);
      }
      else {
        copy(File::Spec->catfile($img_path, 'NoSolution.png'), $to);
      }
    }

    unless(-e $cn_txt) {
      open my $CN, '>', $cn_txt;
      print $CN join "\t", (q{}, 'Chromosome','Position','Log R', 'segmented LogR', 'BAF', 'segmented BAF', 'Copy number', 'Minor allele', 'Raw copy number');
      print "\n";
      close $CN;
    }

  }
  else {
    my $samp_stat_file = sprintf '%s/%s.samplestatistics.txt', $options->{'outdir'}, $tum_name;
    open my $STAT, '>>', $samp_stat_file;
    print $STAT "GenderChr $options->{genderChr}\n";
    print $STAT "GenderChrFound $options->{genderIsMale}\n";
    close $STAT;
  }

  my $new_vcf = $cave_cn;
  $new_vcf =~ s/\.csv$/\.vcf/;
  my $command = "$^X ";
  $command .= _which('ascatCnToVCF.pl');
  $command .= " -o $new_vcf";
  $command .= " -r $options->{reference}";
  $command .= " -i $cave_cn";
  $command .= " -sbm $options->{tumour}";
  $command .= " -sbw $options->{normal}";
  $command .= " -ra $options->{assembly}" if(defined $options->{'assembly'});
  $command .= " -rs $options->{species}" if(defined $options->{'species'});
  $command .= " -msq $options->{protocol} -wsq $options->{protocol}" if(defined $options->{'protocol'});
  $command .= " -msp $options->{platform} -wsp $options->{platform}" if(defined $options->{'platform'});

  my $vcf_gz = $new_vcf.'.gz';
  my $sort_gz = sprintf q{(grep '^#' %s ; grep -v '^#' %s | sort -k 1,1 -k 2,2n) | %s -c > %s}, $new_vcf, $new_vcf, _which('bgzip'), $vcf_gz;

  my $tabix = _which('tabix');
  $tabix .= sprintf ' -p vcf %s', $vcf_gz;

  push @commands, $command, $sort_gz, $tabix;

  my $cn_txt_gz = qq{gzip -c $cn_txt > $cn_txt.gz};
  push @commands, $cn_txt_gz;

  unless($options->{'nb'}) {
    my $aliases = sex_chr_mapping($options->{'snp_gc'});

    my $cn_to_bw = "$^X ";
    $cn_to_bw .= _which('ascatToBigWig.pl');
    $cn_to_bw .= " -f $options->{reference}.fai";
    $cn_to_bw .= " -i $cn_txt.gz";
    $cn_to_bw .= " -o $bw_stub";
    $cn_to_bw .= " -a $aliases" if(defined $aliases && length $aliases > 0);

    push @commands, $cn_to_bw;
  }

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);

  unlink $new_vcf, $cn_txt;

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub clean_snp_gc {
  my ($snp_gc, $cleaned) = @_;
  open my $GC_IN, '<', $snp_gc or die "ERROR: Failed to open $snp_gc\n\t$!\n";
  open my $GC_OUT,'>', $cleaned or die "ERROR: Failed to create $cleaned\n\t$!\n";
  while (my $l = <$GC_IN>) {
    my @F = split /\t/, $l;
    $F[1] =~ s/^chr//i unless($. == 1);
    print $GC_OUT join qq{\t}, @F;
  }
  close $GC_IN;
  close $GC_OUT;
}

sub merge_counts {
  my ($options, $tmp, $sample, $outfile) = @_;
  my @chrs = snpLociChrs($options);

  open my $OFH, '>', $outfile or die "ERROR: Failed to create $outfile\n\t$!\n";

  my $first = 1;
  for my $chr(@chrs) {
    my $split = File::Spec->catfile($tmp, 'allele_count' ,sprintf '%s.%s.allct', $sample, $chr);
    open my $SF, '<', $split or die "ERROR: Failed to open $split\n\t$!\n";
    while(my $l = <$SF>) {
      if($l =~ m/^#/) {
        if($first == 1) {
          print $OFH $l;
          $first = 0;
        }
        next;
      }
      $l =~ s/^chr//i;
      print $OFH $l;
    }
    close $SF;
  }

  close $OFH;
  return 1;
}

sub get_allele_count_file_path {
  my ($tmp,$sample_name) = @_;
  return File::Spec->catfile(File::Spec->catdir($tmp, $sample_name),'sample.allele_count');
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_\-.]/_/ig; # sanitise sample name
  return $sample;
}

sub prepare {
  my $options = shift;
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumour'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normal'}))[0];
  return 1;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in PATH or local bin folder" unless(defined $path);
  return $path;
}

sub determine_gender {
  my $options = shift;
  my $gender_loci;
  if(defined $options->{'locus'}) {
    $gender_loci = $options->{'locus'};
  }
  else {
    my $mod_path = dirname(abs_path($0)).'/../share';
    $mod_path = module_dir('Sanger::CGP::Ascat::Implement') unless(-e File::Spec->catdir($mod_path, 'gender'));

    my $gender_path = File::Spec->catdir($mod_path, 'gender');
    $gender_loci = File::Spec->catfile($gender_path,'GRCh37d5_Y.loci');
  }

  my $idx = 0;
  $idx = $options->{'index'} if(exists $options->{'index'} && defined $options->{'index'});
  my $outfile = File::Spec->catfile($options->{'tmp'}, $idx.'.normal_gender.tsv');

  my $command = _which('alleleCounter.pl');
  $command .= sprintf $ALLELE_COUNT_GENDER, $options->{'normal'}, $gender_loci, $options->{'reference'};
  $command .= '-m '.$options->{'minbasequal'} if exists $options->{'minbasequal'};
  my ($stdout, $stderr, $exit) = capture { system($command); };
  die "Command failed: $command\nERROR: $stderr\n" if($exit);
  chomp $stdout;
  my ($gender_chr, $norm_gender) = split /\t/, $stdout;
  return ($norm_gender, $gender_chr);
}

=head

Updated to run from SnpGcPositions as this file is the only one technically needed

=cut

sub snpLociChrs {
  my $options = shift;
  my %chr_set;
  my @chrs;
  open my $IN, '<', $options->{'snp_gc'};
  while(my $line = <$IN>) {
    next if($line =~ m/Chr\tPosition/);
    my (undef, $chr, undef) = split /\t/, $line;
    unless(exists $chr_set{$chr}) {
      $chr_set{$chr} = 1;
      push @chrs, $chr;
    }
  }
  close $IN;
  return @chrs;
}

sub limited_indices {
	my ($options, $index_in, $count) = @_;
  my @indicies;
  if(exists $options->{'limit'}) {
    # main script checks index is not greater than limit or < 1
	  my $base = $index_in;
	  while($base <= $count) {
	    push @indicies, $base;
	    $base += $options->{'limit'};
	  }
	}
	else {
	  push @indicies, $index_in;
	}
	return @indicies;
}

sub sex_chr_mapping {
  my $snp_gc = shift;
  open my $I, '<', $snp_gc;
  my $idx = 0;
  my %chr_map;
  my @aliases;
  while (<$I>) {
    next if($. == 1);
    my $chr=(split /\t/)[1];
    unless(exists $chr_map{$chr}) {
      $chr_map{$chr} = 1;
      $idx++;
      next if($chr =~ m/^\d+$/);
      push @aliases, "$idx:$chr";
    }
  }
  return join ',', @aliases;
}


1;
