#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2016 Genome Research Ltd.
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

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);
use File::Copy;

use Sanger::CGP::Ascat::Implement;

use PCAP::Cli;

const my @VALID_PROCESS => qw(allele_count ascat finalise);
my %index_max = ( 'allele_count'   => -1,
                  'ascat' => 1,
                  'finalise' => 1);
const my @VALID_GENDERS => qw(XX XY L);

{
  my $options = setup();
  Sanger::CGP::Ascat::Implement::prepare($options);
  my $threads = PCAP::Threaded->new($options->{'threads'});
  &PCAP::Threaded::disable_out_err if(exists $options->{'index'});

  # register any process that can run in parallel here
  $threads->add_function('allele_count', \&Sanger::CGP::Ascat::Implement::allele_count);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  if(!exists $options->{'process'} || $options->{'process'} eq 'allele_count') {
    my $jobs = $options->{'lociChrsBySample'};
    $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
    $threads->run($jobs, 'allele_count', $options);
  }

  Sanger::CGP::Ascat::Implement::ascat($options) if(!exists $options->{'process'} || $options->{'process'} eq 'ascat');
  if(!exists $options->{'process'} || $options->{'process'} eq 'finalise') {
    Sanger::CGP::Ascat::Implement::finalise($options);
    cleanup($options) unless($options->{'noclean'} == 1);
  }
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  remove_tree $tmpdir if(-e $tmpdir);
  return 0;
}

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'o|outdir=s' => \$opts{'outdir'},
              't|tumour=s' => \$opts{'tumour'},
              'n|normal=s' => \$opts{'normal'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
              'c|cpus=i' => \$opts{'threads'},
              'x|limit=i' => \$opts{'limit'},
              'q|minbasequal=i' => \$opts{'minbasequal'},
              's|snp_loci=s' => \$opts{'snp_loci'}, # still here so it doesn't break if defined
              'sp|snp_pos=s' => \$opts{'snp_pos'}, # still here so it doesn't break if defined
              'sg|snp_gc=s' => \$opts{'snp_gc'},
              'g|gender=s' => \$opts{'gender'},
              'gc|genderChr=s' => \$opts{'genderChr'},
              'l|locus=s' => \$opts{'locus'},
              'r|reference=s' => \$opts{'reference'},
              'rs|species=s' => \$opts{'species'},
              'ra|assembly=s' => \$opts{'assembly'},
              'pr|protocol=s' => \$opts{'protocol'},
              'pl|platform=s' => \$opts{'platform'},
              'pu|purity=s' => \$opts{'purity'},
              'pi|ploidy=s' => \$opts{'ploidy'},
              'f|force' => \$opts{'force'},
              'nc|noclean' => \$opts{'noclean'},
              'nb|nobigwig' => \$opts{'nobigwig'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1, -exitval => 0) if(defined $opts{'h'});
  pod2usage(-verbose => 2, -exitval => 0) if(defined $opts{'m'});

  if($opts{'v'}){
    print Sanger::CGP::Ascat->VERSION."\n";
    exit;
  }

  warn "Executing: $opts{cmd}\n";

  if(!defined($opts{'noclean'})){
    $opts{'noclean'} = 0;
  }

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($defined);

  for my $item (qw(tumour normal snp_gc reference outdir protocol gender)) {
    pod2usage(-msg  => "\nERROR: Option '-$item' must be defined.\n", -verbose => 1,  -output => \*STDERR) unless(defined $opts{$item});
  }

  if((defined($opts{'purity'}) && !defined($opts{'ploidy'})) || (!defined($opts{'purity'}) && defined($opts{'ploidy'}))){
    pod2usage(-msg  => "\nERROR: If one of purity/ploidy are defined, both should be defined.\n", -verbose => 1,  -output => \*STDERR);
  }

  for my $item(qw(tumour normal snp_gc reference outdir locus)) {
    $opts{$item} = File::Spec->rel2abs( $opts{$item} ) if(defined $opts{$item});
  }

  PCAP::Cli::file_for_reading('tumour', $opts{'tumour'});
  PCAP::Cli::file_for_reading('normal', $opts{'normal'});
  PCAP::Cli::file_for_reading('snp_gc', $opts{'snp_gc'});
  PCAP::Cli::file_for_reading('reference', $opts{'reference'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  if(defined $opts{'snp_loci'}) {
    warn qq{NOTE: '-snp_loci' file is no longer required\n};
    delete $opts{'snp_loci'};
  }
  if(defined $opts{'snp_pos'}) {
    warn qq{NOTE: '-snp_pos' file is no longer required\n};
    delete $opts{'snp_pos'};
  }

  my $final_logs = File::Spec->catdir($opts{'outdir'}, 'logs');
  if(-e $final_logs) {
    warn "NOTE: Presence of '$final_logs' directory suggests successful complete analysis, please delete to rerun\n";
    exit 0;
  }

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
  delete $opts{'limit'} unless(defined $opts{'limit'});
  $opts{'minbasequal'} = 20 unless(defined $opts{'minbasequal'});

  $opts{'lociChrsBySample'} = scalar Sanger::CGP::Ascat::Implement::snpLociChrs(\%opts) * 2;

  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    if(exists $opts{'index'}) {
      my $max = $index_max{$opts{'process'}};

      if($max == -1) {
        # can only be alleleCount
        if(exists $opts{'limit'}) {
          $max = $opts{'limit'};
        }
        else {
          $max = $opts{'lociChrsBySample'}
        }
      }

      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
      die "No max has been defined for this process type\n" if($max == 0);
      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max, 1);
    }
  }
  elsif(exists $opts{'index'}) {
    pod2usage(-msg  => "\nERROR: -index cannot be defined without -process.\n", -verbose => 1,  -output => \*STDERR);
  }

  $opts{'threads'} = 1 unless(defined $opts{'threads'});

  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpAscat');
  make_path($tmpdir) unless(-d $tmpdir);
  my $progress = File::Spec->catdir($tmpdir, 'progress');
  make_path($progress) unless(-d $progress);
  my $logs = File::Spec->catdir($tmpdir, 'logs');
  make_path($logs) unless(-d $logs);

  $opts{'tmp'} = $tmpdir;

  if(defined $opts{'gender'}){
    pod2usage(-message => "\nERROR: Unknown gender value: $opts{gender}\n", -verbose => 1) unless(first {$_ eq $opts{'gender'}} @VALID_GENDERS);
    if($opts{'gender'} eq 'L') {
      my ($is_male, $gender_chr) = Sanger::CGP::Ascat::Implement::determine_gender(\%opts);
      $opts{'genderChr'} = $gender_chr;
      $opts{'genderIsMale'} = $gender_chr;
      $opts{'gender'} = $is_male eq 'N' ? 'XX' : 'XY';
    }
    else {
      pod2usage(-message => "\nERROR: genderChr must be set when gender is XX/XY\n", -verbose => 1) if(!defined $opts{'genderChr'});
      pod2usage(-message => "\nERROR: gender must be XX, XY or L\n", -verbose => 1)if($opts{'gender'} !~ m/^X[XY]$/);
      $opts{'genderIsMale'} = $opts{'gender'} eq 'XX' ? 'N' : 'Y';
    }
  } else {
    pod2usage(-message => "\nERROR: gender not set\n", -verbose => 1);
  }

  return \%opts;
}



__END__

=head1 ascat.pl

Reference implementation of Cancer Genome Project Ascat
copy-number analysis pipeline.

=head1 SYNOPSIS

ascat.pl [options]

  Please define as many of the parameters as possible

  Required parameters

    -outdir       -o    Folder to output result to.
    -tumour       -t    Tumour BAM/CRAM file
    -normal       -n    Normal BAM/CRAM file
    -reference    -r    Reference fasta
    -snp_gc       -sg   Snp GC correction file
    -protocol     -pr   Sequencing protocol (e.g. WGS, WXS)
    -gender       -g    Sample gender (XX, XY, L)
                          For XX/XY see '-gc'
                          When 'L' see '-l'

  Targeted processing (further detail under OPTIONS):
    -process      -p    Only process this step then exit, optionally set -index
    -index        -i    Optionally restrict '-p' to single job
    -limit        -x    Specifying 2 will balance processing between '-i 1 & 2'
                        Must be paired with '-p allele_count'

  Optional parameters
    -genderChr    -gc   Specify the 'Male' sex chromosome: Y,chrY...
    -species      -rs   Reference species [BAM HEADER]
    -assembly     -ra   Reference assembly [BAM HEADER]
    -platform     -pl   Seqeuncing platform [BAM HEADER]
    -minbasequal  -q    Minimum base quality required before allele is used. [20]
    -cpus         -c    Number of cores to use. [1]
                        - recommend max 2 during 'input' process.
    -locus        -l    Using a list of loci, default when '-L' [share/gender/GRCh37d5_Y.loci]
                        - these are loci that will not present at all in a female sample
    -force        -f    Force completion - solution not possible
                        - adding this will result in successful completion of analysis even
                          when ASCAT can't generate a solution.  A default copynumber of 5/2
                          (tumour/normal) and contamination of 30% will be set along with a
                          comment in '*.samplestatistics.csv' to indicate this has occurred.
    -purity       -pu   Purity (rho) setting for manual setting of sunrise plot location
    -ploidy       -pi   Ploidy (psi) setting for manual setting of sunrise plot location
    -noclean      -nc   Finalise results but don't clean up the tmp directory.
                        - Useful when including a manual check and restarting ascat with new pu and pi params.
    -nobigwig     -nb   Don't generate BigWig files.

  Other
    -help         -h    Brief help message
    -man          -m    Full documentation.
    -version      -v    Ascat version number
