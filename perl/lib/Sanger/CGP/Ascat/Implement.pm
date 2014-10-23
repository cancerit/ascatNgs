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
use FindBin qw($Bin);

use File::ShareDir qw(module_dir);

use Sanger::CGP::Ascat;

use PCAP::Threaded;
use PCAP::Bam;

const my $COUNT_READS => q{%s view -c %s %s};

const my $FAILED_SAMPLE_STATISTICS => qq{## WARNING ASCAT failed to generate a solution ##\nNormalContamination 0.3\nPloidy ?\nrho 0\npsi 0goodnessOfFit 0\n};

const my $ALLELE_COUNT_PARA => ' -b %s -o %s -l %s ';

const my @ASCAT_RESULT_FILES => qw( %s.aberrationreliability.png
                                    %s.ASCATprofile.png
                                    %s.ASPCF.png
                                    %s.germline.png
                                    %s.rawprofile.png
                                    %s.sunrise.png
                                    %s.tumour.png
                                    %s.copynumber.caveman.csv
                                    %s.copynumber.txt
                                    %s.samplestatistics.csv
                                  );

const my $GENDER_MIN => 5;

sub allele_count {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my $ac_out = File::Spec->catdir($tmp, 'allele_count');
  make_path($ac_out) unless(-e $ac_out);

  my @inputs = ($options->{'tumour'}, $options->{'normal'});
  my $iter = 1;
  for my $input(@inputs) {
    next if($iter++ != $index); # skip to the relevant input in the list

    my $sname = sanitised_sample_from_bam($input);
    my $alleleCountOut = File::Spec->catfile($ac_out,$sname .'.allct');

    my $allc_exe = _which('alleleCounter');
    my $allc_lib = dirname($allc_exe);

    my $command = $allc_exe;
    $command .= sprintf $ALLELE_COUNT_PARA, $input, $alleleCountOut, $options->{'snp_loci'};
    $command .= '-m '.$options->{'minbasequal'} if exists $options->{'minbasequal'};

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
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

  unlink $tumcount if -l $tumcount;
  unlink $normcount if -l $normcount;

  my $lntc = 'ln -s '.File::Spec->catfile(File::Spec->catdir($tmp, 'allele_count'),$tum_name.'.allct')." $tumcount";
  my $lnnc =  'ln -s '.File::Spec->catfile(File::Spec->catdir($tmp, 'allele_count'),$norm_name.'.allct')." $normcount";

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $lntc, 0);
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $lnnc, 0);

  my $command = "cd $ascat_out; "._which('Rscript');

  my $mod_path = dirname(abs_path($0)).'/../share';
  $mod_path = module_dir('Sanger::CGP::Ascat::Implement') unless(-e File::Spec->catdir($mod_path, 'ascat'));

  my $ascat_path = File::Spec->catdir($mod_path, 'ascat');
  my $ascat_exe = File::Spec->catfile($ascat_path,'runASCAT.R');

  $command .= " $ascat_exe";
  $command .= " $ascat_path";
  $command .= ' '.$options->{'snp_pos'};
  $command .= ' '.$options->{'snp_gc'};
  $command .= ' '.$options->{'tumour_name'};
  $command .= ' '.$tumcountfile;
  $command .= ' '.$options->{'normal_name'};
  $command .= ' '.$normcountfile;
  $command .= ' '.$options->{'gender'};
  $command .= ' '.$rdata;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

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
  foreach my $f(@ASCAT_RESULT_FILES){
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
  if($force_complete == 1) {
    my $fake_file = sprintf '%s/%s.copynumber.caveman.csv', $options->{'outdir'}, $tum_name;
    my $fake_csv = "$^X ";
    $fake_csv .= _which('failed_cn_csv.pl');
    $fake_csv .= sprintf ' -r %s -o %s', $options->{'reference'}, $fake_file;
    push @commands, $fake_csv;
    my $samp_stat_file = sprintf '%s/%s.samplestatistics.csv', $options->{'outdir'}, $tum_name;
    open my $STAT, '>', $samp_stat_file;
    print $STAT $FAILED_SAMPLE_STATISTICS or die "Failed to write line to $samp_stat_file";
    close $STAT;
    $cave_cn = $fake_file;
  }
  my $new_vcf = $cave_cn;
  $new_vcf =~ s/\.csv$/\.vcf/;
  my $command = "$^X ";
  $command .= _which('CN_to_VCF.pl');
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
  my $bgzip = _which('bgzip');
  $bgzip .= sprintf ' -c %s > %s', $new_vcf, $vcf_gz;

  my $tabix = _which('tabix');
  $tabix .= sprintf ' -p vcf %s', $vcf_gz;

  push @commands, $command, $bgzip, $tabix;

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);

  unlink $new_vcf;

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub get_allele_count_file_path {
  my ($tmp,$sample_name) = @_;
  return File::Spec->catfile(File::Spec->catdir($tmp, $sample_name),'sample.allele_count');
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_-]/_/ig; # sanitise sample name
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
  my $gender = 'XX';
  my $samtools = _which('samtools');
  my $command = sprintf $COUNT_READS, $samtools, $options->{'tumour'}, $options->{'locus'};
  my $tumour_count = `$command`;
  chomp $tumour_count;
  $command = sprintf $COUNT_READS, $samtools, $options->{'normal'}, $options->{'locus'};
  my $normal_count = `$command`;
  chomp $normal_count;
  if($normal_count > $GENDER_MIN && $tumour_count > $GENDER_MIN) {
    $gender = 'XY'; # male
  }
  elsif($normal_count <= $GENDER_MIN && $tumour_count <= $GENDER_MIN) {
    $gender = 'XX'; # female
  }
  else {
    die "Gender guess inconclusive, tum depth: $tumour_count - norm depth: $normal_count\n";
  }
  return $gender;
}

1;
