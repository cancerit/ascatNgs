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

use strict;
use Bio::DB::HTS;
use Bio::DB::HTS::Faidx;

use Const::Fast qw(const);

const my $HEADER => "\tChr\tPosition\tProbe\t200bp\t400bp\t800bp\t1600bp\t3200bp\t6400bp\t12800bp\t25600bp\t51200bp\t102400bp\t204800bp\t1M\t2M\t5M\t10M\n";

const my @GC_BINS => (24,
                      200,
                      400,
                      800,
                      1_600,
                      3_200,
                      6_400,
                      12_800,
                      25_600,
                      51_200,
                      102_400,
                      204_800,
                      1_000_000,
                      2_000_000,
                      5_000_000,
                      1_000_0000);

die "USAGE: ascatSnpPanelGcCorrections.pl genome.fa SnpPositions.tsv" unless(@ARGV == 2);

my $ref = shift @ARGV; # reference.fa
my $snps = shift @ARGV; # SnpPositions

my $index = Bio::DB::HTS::Faidx->new($ref); # has horrid memory leak

print $HEADER;

open my $SNPIN, '<', $snps or die $!;
while(my $line = <$SNPIN>) {
  next if($line =~ /^\t/);
  chomp $line;
  my ($id, $chr, $pos) = split /\t/, $line;

  my $seq_len = $index->length($chr);

  my @gc_set = ($id, $chr, $pos);
  for my $bin(@GC_BINS) {
    my $pad = $bin/2;
    my $start = $pos - $pad;
    my $end = $pos + $pad;
    $start = 1 if($start < 1);
    $end = $seq_len if($end > $seq_len);
    my $range = sprintf '%s:%d-%d', $chr, $start, $end;

    my ($seq, $len) = $index->get_sequence($range);

    my $gc = $seq =~ tr/GCgc/GCgc/;
    my $gc_pct = sprintf '%.2f', ($gc/$len) * 100;
    push @gc_set, $gc_pct;
  }
  print join("\t", @gc_set),"\n";
}
close $SNPIN;

__END__
probe = +/-14bp
200bp
400bp
800bp
1600bp
3200bp
6400bp
12800bp
25600bp
51200bp
102400bp
204800bp
1M
2M
5M
10M
