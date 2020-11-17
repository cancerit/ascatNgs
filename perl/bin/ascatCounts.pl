#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2020 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
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
use File::Which qw(which);
use FindBin qw($Bin);

use Sanger::CGP::Ascat::Implement;

use PCAP::Cli;

{
  my $options = setup();
  $options->{'tumour'} = $options->{'bam'}; #map input file to tumour
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'bam'}))[0];
  my $threads = PCAP::Threaded->new($options->{'threads'});

  # register any process that can run in parallel here
  $threads->add_function('allele_count', \&Sanger::CGP::Ascat::Implement::allele_count);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  my $jobs = $options->{'lociChrsBySample'};
  $jobs = $options->{'limit'} if(exists $options->{'limit'} && defined $options->{'limit'});
  $threads->run($jobs, 'allele_count', $options);

  Sanger::CGP::Ascat::Implement::merge_counts_and_index($options);
  cleanup($options) unless($options->{'noclean'} == 1);
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  move(File::Spec->catdir($tmpdir, 'ascatCounts'), File::Spec->catdir($options->{'outdir'}, 'ascatCounts')) || die $!;
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  remove_tree $tmpdir if(-e $tmpdir);
  return 0;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in PATH or local bin folder" unless(defined $path);
  return $path;
}

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'o|outdir=s' => \$opts{'outdir'},
              'b|bam=s' => \$opts{'bam'},
              'c|cpus=i' => \$opts{'threads'},
              'q|minbasequal=i' => \$opts{'minbasequal'},
              'sg|snp_gc=s' => \$opts{'snp_gc'},
              'r|reference=s' => \$opts{'reference'},
              'nc|noclean' => \$opts{'noclean'},
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

  for my $item (qw(bam snp_gc reference outdir)) {
    pod2usage(-msg  => "\nERROR: Option '-$item' must be defined.\n", -verbose => 1,  -output => \*STDERR) unless(defined $opts{$item});
  }

  for my $item(qw(bam reference outdir)) {
    $opts{$item} = File::Spec->rel2abs( $opts{$item} ) if(defined $opts{$item});
  }

  PCAP::Cli::file_for_reading('bam', $opts{'bam'});
  PCAP::Cli::file_for_reading('snp_gc', $opts{'snp_gc'});
  PCAP::Cli::file_for_reading('reference', $opts{'reference'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  my $final_logs = File::Spec->catdir($opts{'outdir'}, 'logs');
  if(-e $final_logs) {
    warn "NOTE: Presence of '$final_logs' directory suggests successful complete analysis, please delete to rerun\n";
    exit 0;
  }

  $opts{'minbasequal'} = 20 unless(defined $opts{'minbasequal'});

  $opts{'lociChrsBySample'} = scalar Sanger::CGP::Ascat::Implement::snpLociChrs(\%opts);

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

ascatCounts.pl [options]

  Please define as many of the parameters as possible

  Required parameters

    -outdir       -o    Folder to output result to.
    -bam          -b    Input BAM/CRAM file
    -reference    -r    Reference fasta
    -snp_gc       -sg   Snp GC correction file

  Optional parameters
    -minbasequal  -q    Minimum base quality required before allele is used. [20]
    -cpus         -c    Number of cores to use. [1]
    -noclean      -nc   Finalise results but don't clean up the tmp directory.

  Other
    -help         -h    Brief help message
    -man          -m    Full documentation.
    -version      -v    Ascat version number
