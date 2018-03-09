#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2018 Genome Research Ltd.
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of ascatNgs.
#
# cgpVcf is free software: you can redistribute it and/or modify it under
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
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use Getopt::Long;
use Pod::Usage qw(pod2usage);

use Try::Tiny;
use PCAP::Cli;

{
  my $opts = setup();

  my ($output_loc,$OUT_FH);
  #Open out file (OR STDOUT depending on what's been given at command line)
  if($opts->{'o'}) {
    open($OUT_FH, '>', $opts->{'o'});
    $output_loc = $opts->{'o'};
  }
  else{
    $OUT_FH = \*STDOUT;
    $output_loc = 'STDOUT';
  }

  open my $IN_FH, '<', $opts->{'r'}.'.fai' or die "Failed to read from $opts->{r}.fai: $!";
  my ($start,$wt_cn_tot,$wt_cn_min,$mt_cn_tot,$mt_cn_min) = (1,2,1,5,2);
  while(my $line = <$IN_FH>){
    my $seg_no = $.;
    chomp $line;
    my ($chr, $end) = (split "\t", $line)[0,1];
    print $OUT_FH join(q{,}, $seg_no, $chr, $start, $end, $wt_cn_tot, $wt_cn_min, $mt_cn_tot, $mt_cn_min),"\n";
  }
  close $IN_FH;
  close($OUT_FH) if($opts->{'o'});

}

sub setup{
  my %opts;
  #Store the command used to run this script.
  $opts{'cmd'} = join " ", $0, @ARGV;
  my @random_args;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'v|version' => \$opts{'v'},
					'r|reference=s' => \$opts{'r'},
					'o|out-file=s' => \$opts{'o'},
          '<>' => sub{push(@random_args,shift(@_));}
  ) or pod2usage(2);

  if(defined $opts{'v'}){
    print "Version: ".Sanger::CGP::Ascat->VERSION."\n";
    exit;
  }

  PCAP::Cli::file_for_reading('r', $opts{'r'});

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  return \%opts;
}


__END__

=head1 NAME

ascatFailedCnCsv.pl - Generate CopyNumber csv for a failed ASCAT run suitable for caveman, 5/2 (tumour/normal).

=head1 SYNOPSIS

ascatFailedCnCsv.pl [options]

    Required parameters:
      -reference  -r  Path to reference fasta (with corresponding fai)

    Optional parameters:
      -out-file   -o    Output file [STDOUT].

    Other:
      -help       -h   Brief help message.
      -man        -m   Full documentation.
      -version    -v   Version information.


=head1 PARAMETERS

=over 8

=item B<-reference>

Reference file

=item B<-out-file>

File to write VCF format output to. If not provided STDOUT is used.

=item B<-version>

Print version information.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<ascatFailedCnCsv.pl> generates CopyNumber csv for a failed ASCAT run suitable for caveman, 5/2 (tumour/normal).

=cut
