#!/usr/bin/perl

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

use Data::Dumper;

use Sanger::CGP::Ascat::Implement;

use PCAP::Cli;

const my @VALID_PROCESS => qw(allele_count ascat finalise);
my %index_max = ( 'allele_count'   => 2,
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
  $threads->run(2, 'allele_count', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'allele_count');

  Sanger::CGP::Ascat::Implement::ascat($options) if(!exists $options->{'process'} || $options->{'process'} eq 'ascat');
  if(!exists $options->{'process'} || $options->{'process'} eq 'finalise') {
    Sanger::CGP::Ascat::Implement::finalise($options);
    cleanup($options);
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
              'q|minbasequal=i' => \$opts{'minbasequal'},
              's|snp_loci=s' => \$opts{'snp_loci'},
              'sp|snp_pos=s' => \$opts{'snp_pos'},
              'sg|snp_gc=s' => \$opts{'snp_gc'},
              'g|gender=s' => \$opts{'gender'},
              'l|locus=s' => \$opts{'locus'},
              'r|reference=s' => \$opts{'reference'},
              'rs|species=s' => \$opts{'species'},
              'ra|assembly=s' => \$opts{'assembly'},
              'pr|protocol=s' => \$opts{'protocol'},
              'pl|platform=s' => \$opts{'platform'},
  ) or pod2usage(2);

  pod2usage(-message => Sanger::CGP::Ascat::license, -verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if($opts{'v'}){
    print Sanger::CGP::Ascat->VERSION."\n";
    exit;
  }

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($defined);

  PCAP::Cli::file_for_reading('tumour', $opts{'tumour'});
  PCAP::Cli::file_for_reading('normal', $opts{'normal'});
  PCAP::Cli::file_for_reading('snp_loci', $opts{'snp_loci'});
  PCAP::Cli::file_for_reading('snp_pos', $opts{'snp_pos'});
  PCAP::Cli::file_for_reading('snp_gc', $opts{'snp_gc'});
  PCAP::Cli::file_for_reading('reference', $opts{'reference'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
  $opts{'minbasequal'} = 20 unless(defined $opts{'minbasequal'});

  if(defined $opts{'gender'}){
    pod2usage(-message => 'unknown gender value: '.$opts{'gender'}, -verbose => 1) unless(first {$_ eq $opts{'gender'}} @VALID_GENDERS);
    if($opts{'gender'} eq 'L') {
      $opts{'gender'} = Sanger::CGP::Ascat::Implement::determine_gender(\%opts);
    }
  } else {
    pod2usage(-message => 'gender not set', -verbose => 1);
  }


  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    if(exists $opts{'index'}) {
      my $max = $index_max{$opts{'process'}};
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
      die "No max has been defined for this process type\n" if($max == 0);
      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max, 1);
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

  $opts{'threads'} = 1 unless(defined $opts{'threads'});

  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpAscat');
  make_path($tmpdir) unless(-d $tmpdir);
  my $progress = File::Spec->catdir($tmpdir, 'progress');
  make_path($progress) unless(-d $progress);
  my $logs = File::Spec->catdir($tmpdir, 'logs');
  make_path($logs) unless(-d $logs);

  $opts{'tmp'} = $tmpdir;

  return \%opts;
}



__END__

=head1 ascat.pl

Reference implementation of Cancer Genome Project Ascat
copy-number analysis pipeline.

=head1 SYNOPSIS

ascat.pl [options]

  Please defined as many of the parameters as possible

  Required parameters

    -outdir       -o    Folder to output result to.
    -tumour       -t    Tumour BAM file
    -normal       -n    Normal BAM file
    -reference    -r    Reference fasta
    -snp_loci     -s    Snp locus file
    -snp_pos      -sp   Snp position file
    -snp_gc       -sg   Snp GC correction file
    -gender       -g    Sample gender (XX, XY, L)
                          When 'L' define '-l'

  Targeted processing (further detail under OPTIONS):
    -process      -p    Only process this step then exit, optionally set -index
    -index        -i    Optionally restrict '-p' to single job

  Optional parameters
    -species      -rs   Reference species [BAM HEADER]
    -assembly     -ra   Reference assembly [BAM HEADER]
    -protocol     -pr   Sequencing protocol (e.g. WGS, WXS)
    -platform     -pl   Seqeuncing platform [BAM HEADER]
    -minbasequal  -q    Minimum base quality required before allele is used. [20]
    -cpus         -c    Number of cores to use. [1]
                        - recommend max 2 during 'input' process.
    -locus        -l    Attempt to determine gender using a male specific locus.
                          e.g. Y:2654896-2655740 (GRCh37)


  Other
    -help         -h    Brief help message
    -man          -m    Full documentation.
    -version      -v    Ascat version number
