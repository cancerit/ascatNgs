#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2018 Genome Research Ltd.
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

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use Const::Fast qw(const);

use PCAP::Cli;
use Sanger::CGP::Ascat;

# chr, 0-start, 1-stop, value (use string to leave precision as source)
const my $BG_FRMT => "%s\t%d\t%d\t%s\n";
const my $DFLT_ALIAS => '23:X,24:Y';

my $options = setup();
my $chr_map = chr_lookup($options->{'alias'});

my $z = new IO::Uncompress::Gunzip $options->{'input'}, MultiStream=>1 or die "gunzip failed: $GunzipError\n";

my $bg_fhs = get_output_files($options->{'outpre'});

my $chr_last = q{};
my ($start_segbaf, $end_segbaf, $val_segbaf) = (0,0,'NA');
my ($start_seglogr, $end_seglogr, $val_seglogr) = (0,0,'NA');
my ($start_total_cn, $end_total_cn, $val_total_cn) = (0,0,'NA');
my ($start_minor_cn, $end_minor_cn, $val_minor_cn) = (0,0,'NA');

while(my $line = <$z>) {
  next if($line =~ /^\t/);
  chomp $line;
  my (undef, $chr, $pos, $logR, $segLogR, $baf, $segBaf, $total_cn, $minor_cn, $raw_cn) = split /\t/, $line;
  $chr = $chr_map->{$chr} if(exists $chr_map->{$chr});

  if($chr ne $chr_last) {
    if($chr_last ne q{}) {
      printf {$bg_fhs->{'segLogR'}} $BG_FRMT, $chr_last, $start_seglogr-1, $end_seglogr, $val_seglogr if($val_seglogr ne 'NA');
      printf {$bg_fhs->{'segBaf'}} $BG_FRMT, $chr_last, $start_segbaf-1, $end_segbaf, $val_segbaf if($val_segbaf ne 'NA');
      printf {$bg_fhs->{'totalCn'}} $BG_FRMT, $chr_last, $start_total_cn-1, $end_total_cn, $val_total_cn if($val_total_cn ne 'NA');
      printf {$bg_fhs->{'minorCn'}} $BG_FRMT, $chr_last, $start_minor_cn-1, $end_minor_cn, $val_minor_cn if($val_minor_cn ne 'NA');
    }
    $chr_last = $chr;
    ($start_segbaf, $end_segbaf, $val_segbaf) = (0,0,'NA');
    ($start_seglogr, $end_seglogr, $val_seglogr) = (0,0,'NA');
    ($start_total_cn, $end_total_cn, $val_total_cn) = (0,0,'NA');
    ($start_minor_cn, $end_minor_cn, $val_minor_cn) = (0,0,'NA');
  }

  printf {$bg_fhs->{'logR'}} $BG_FRMT, $chr, $pos-1, $pos, $logR if($logR ne 'NA');
  printf {$bg_fhs->{'baf'}} $BG_FRMT, $chr, $pos-1, $pos, $baf if($baf ne 'NA');
  printf {$bg_fhs->{'rawCn'}} $BG_FRMT, $chr, $pos-1, $pos, $raw_cn if($raw_cn ne 'NA');

  if($segLogR ne $val_seglogr) {
    printf {$bg_fhs->{'segLogR'}} $BG_FRMT, $chr_last, $start_seglogr-1, $end_seglogr, $val_seglogr if($val_seglogr ne 'NA');
    ($start_seglogr, $end_seglogr, $val_seglogr) = ($pos, $pos, $segLogR);
  }
  else {
    $end_seglogr = $pos;
  }

  if($segBaf ne $val_segbaf) {
    printf {$bg_fhs->{'segBaf'}} $BG_FRMT, $chr_last, $start_segbaf-1, $end_segbaf, $val_segbaf if($val_segbaf ne 'NA');
    ($start_segbaf, $end_segbaf, $val_segbaf) = ($pos, $pos, $segBaf);
  }
  else {
    $end_segbaf = $pos;
  }

  if($total_cn ne $val_total_cn) {
    printf {$bg_fhs->{'totalCn'}} $BG_FRMT, $chr_last, $start_total_cn-1, $end_total_cn, $val_total_cn if($val_total_cn ne 'NA');
    ($start_total_cn, $end_total_cn, $val_total_cn) = ($pos, $pos, $total_cn);
  }
  else {
    $end_total_cn = $pos;
  }

  if($minor_cn ne $val_minor_cn) {
    printf {$bg_fhs->{'minorCn'}} $BG_FRMT, $chr_last, $start_minor_cn-1, $end_minor_cn, $val_minor_cn if($val_minor_cn ne 'NA');
    ($start_minor_cn, $end_minor_cn, $val_minor_cn) = ($pos, $pos, $minor_cn);
  }
  else {
    $end_minor_cn = $pos;
  }
}
close $z;

printf {$bg_fhs->{'segLogR'}} $BG_FRMT, $chr_last, $start_seglogr-1, $end_seglogr, $val_seglogr if($val_seglogr ne 'NA');
printf {$bg_fhs->{'segBaf'}} $BG_FRMT, $chr_last, $start_segbaf-1, $end_segbaf, $val_segbaf if($val_segbaf ne 'NA');
printf {$bg_fhs->{'totalCn'}} $BG_FRMT, $chr_last, $start_total_cn-1, $end_total_cn, $val_total_cn if($val_total_cn ne 'NA');
printf {$bg_fhs->{'minorCn'}} $BG_FRMT, $chr_last, $start_minor_cn-1, $end_minor_cn, $val_minor_cn if($val_minor_cn ne 'NA');

close_output_files($options->{'outpre'}, $bg_fhs);
convert_bg_to_bw($options->{'outpre'}, $options->{'fai'});

sub chr_lookup {
  my $alias_lst = shift;
  my %chr_map = split /[:,]/, $alias_lst;
  return \%chr_map;
}

sub convert_bg_to_bw {
  my ($outpre, $chromlst) = @_;
  for (qw(logR segLogR baf segBaf totalCn minorCn rawCn)) {
    my $bg_file = sprintf '%s.%s.bg', $outpre, $_;
    my $bw_file = sprintf '%s.%s.bw', $outpre, $_;
    my $cmd = sprintf 'bg2bw -i %s -o %s -c %s', $bg_file, $bw_file, $chromlst;
    warn "Executing: $cmd\n";
    system($cmd) && die "Execution failed: $!";
    unlink $bg_file;
  }
}

sub get_output_files {
  my $outpre = shift;
  my %fhs;
  for (qw(logR segLogR baf segBaf totalCn minorCn rawCn)) {
    my $bg_file = sprintf '%s.%s.bg', $outpre, $_;
    open my $fh, '>', $bg_file or die "Failed to create $bg_file";
    $fhs{$_} = $fh;
  }
  return \%fhs;
}

sub close_output_files {
  my ($outpre, $fhs) = @_;
  for (keys %{$fhs}) {
    my $bg_file = sprintf '%s.%s.bg', $outpre, $_;
    close $fhs->{$_} or die "Failed to close $bg_file";
  }
}

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  warn "Executing: $opts{cmd}\n";
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'i|input=s' => \$opts{'input'},
              'o|outpre=s' => \$opts{'outpre'},
              'f|fai=s' => \$opts{'fai'},
              'a|alias:s' => \$opts{'alias'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1, -exitval => 0) if(defined $opts{'h'});
  pod2usage(-verbose => 2, -exitval => 0) if(defined $opts{'m'});

  if($opts{'v'}){
    print Sanger::CGP::Ascat->VERSION."\n";
    exit;
  }

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Parameters must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($defined);


  PCAP::Cli::file_for_reading('input', $opts{'input'});
  PCAP::Cli::file_for_reading('fai', $opts{'fai'});
  pod2usage(-msg  => "\nERROR: '-o' must be defined.\n", -verbose => 1,  -output => \*STDERR) unless(defined $opts{'outpre'});

  $opts{'alias'} = $DFLT_ALIAS unless(defined $opts{'alias'});

  return \%opts;
}

__END__

=head1 ascatToBigWig.pl

Converts raw data used in log and BAF plots to BigWig data files
suitable for use with JBrowse/MuiltBigWig.

=head1 SYNOPSIS

ascatToBigWig.pl [options]

  Required parameters

    -input    -i  Input *.ascat_ngs.cn.tsv[.gz]
    -outpre   -o  File name to use in output, suffixes applied are:
                   - *.logR.bw
                   - *.segLogR.bw
                   - *.baf.bw
                   - *.segBaf.bw
    -fai      -f  Path to chrom list or fai file
                   - a tsv file with first 2 columns chr, chr-length
    -alias    -a  csv for numeric to chr mapping applied to chr column of input [23:X,24:Y]


  Other
    -help     -h  Brief help message
    -man      -m  Full documentation.
    -version  -v  Ascat version number
