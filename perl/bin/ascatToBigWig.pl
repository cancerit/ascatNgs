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

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use FindBin qw($Bin);
use lib "$Bin/../lib";

use Getopt::Long;
use Pod::Usage qw(pod2usage);
use PCAP::Cli;


use Sanger::CGP::Ascat;

my $options = setup();




sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  warn "Executing: $opts{cmd}\n";
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'v|version' => \$opts{'v'},
              'i|input' => \$opts{'input'},
              'o|outdir' => \$opts{'outdir'},
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
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($defined);


  PCAP::Cli::file_for_reading('input', $opts{'input'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  return \%opts;
}
