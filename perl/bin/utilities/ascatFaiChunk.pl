#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014-2019 Genome Research Ltd.
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
use List::Util qw(first);

if(@ARGV < 2) {
  die "\tUSAGE:\n\tperl fai_chunk.pl genome.fa \$MB_PER_SEC [SKIP_CHR...] > genome_chunks.bed";
}

my $MB = 1_000_000;

my $fai = shift @ARGV;
my $mb = shift @ARGV;
my @skip = @ARGV;

$fai = $fai.'.fai' unless($fai =~ m/[.]fai$/);

my $step = $mb * $MB;

open my $IN, '<', $fai or die $!;
while(my $line = <$IN>) {
  chomp $line;
  my ($chr, $len) = (split /\t/, $line)[0..1];
  next if(first { $chr eq $_ } @skip);
  my $curr = 1;
  while($curr < $len) {
    if($curr + $step >= $len) {
      printf "%s\t%d\t%d\n", $chr, $curr - 1, $len;
      last;
    }
    else {
      printf "%s\t%d\t%d\n", $chr, $curr - 1, $curr+$step;
      $curr = $curr+$step;
    }
  }
}
close $IN;
