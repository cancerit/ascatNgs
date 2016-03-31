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
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $usage = <<'END_USAGE';

USAGE: ascatSnpPanelFromVcfs.pl snps.vcf[.gz]

  This script was created to generate a shared SNP panel for species having
  multiple strains which have become homozygos through in-breeding.

  The resulting output is likely to only be useful in experiments where crossing
  of strains has been performed to add hetrozygous SNPs back into the population.

  The initial use has been against the Mouse Genome Project outputs found here:
    ftp://ftp-mouse.sanger.ac.uk/REL-*-SNPs_Indels/mgp.*.merged.snps_all.*.vcf.gz

  As not all VCF files are created equal you may need to make modifications for
  other sources.

END_USAGE

my $infile = shift @ARGV;

die $usage unless(defined $infile);


my $z = new IO::Uncompress::Gunzip $infile, AutoClose => 1, MultiStream => 1
        or die "gunzip failed: $GunzipError\n";

my $no_id = 0;
my $sample_names;
my $sample_count;
my %sample_dist;
print "\tChr\tPosition\n";
my $last_kb =  -1;
my $cluster = {};
while(my $line = <$z>) {
  next if($line =~ m/^##/);
  chomp $line;
  if($line =~ m/^#/) {
    my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @l_samples) = split /\t/, $line;
    $sample_names = \@l_samples;
    $sample_count = @l_samples -1;
    next;
  }
  my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @sample_gts) = split /\t/, $line;

  my $kb = int ($pos / 2000);
  if($last_kb != $kb) {
    my $result = evaluate_set($cluster, \%sample_dist);
    $cluster = {};
    print $result if($result);

    $last_kb = $kb;
  }

  next if($id eq q{.}); # may want these later, but for now generate without

  next if($filter ne 'PASS');
  next if(length $ref != 1 || length $alt != 1);

  my $smpls_w_snp = sample_count(\@sample_gts, $sample_names, $sample_count);
  my $samp_count = scalar @{$smpls_w_snp};
  next if($samp_count < 3);

  push @{$cluster->{$samp_count}}, [$chr, $pos, $id, $smpls_w_snp];
}
my $result = evaluate_set($cluster, \%sample_dist);
print $result if($result);

warn Dumper(\%sample_dist);

sub evaluate_set {
  my ($cluster, $sample_dist) = @_;
  my $retval;
  my @s_counts = keys %{$cluster};
  if(@s_counts) {
    my $max = (sort {$a<=>$b} @s_counts)[-1];
    my ($l_chr, $l_pos, $l_id, $smpls_w_snp) = @{$cluster->{$max}->[0]};
    for my $n(@{$smpls_w_snp}) {
      $sample_dist{$n}++;
    }
    $retval = join("\t", $l_id, $l_chr, $l_pos)."\n";
  }
  return $retval;
}

sub sample_count {
  my ($gts, $samples, $iter_max) = @_;
  my @names;
  for my $i(0..$iter_max) {
    # sample must start '1/1' and end ':1'
    if($gts->[$i] =~ m|^1/1:.*:1$|) {
      push @names, $samples->[$i];
    }
  }
  return \@names;
}
