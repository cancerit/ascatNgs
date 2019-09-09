#!/usr/bin/perl

########## LICENCE ##########
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
########## LICENCE ##########

use strict;
use Bio::DB::HTS::Faidx;

my $usage = <<'USAGE_END';

USAGE: perl ascatSnpPanelMerge.pl genome.fa ascatSnpPanelGeneratorOut.1 [..N]

  Takes the output files from ascatSnpPanelGenerator.pl and identifies common
  SNPs from the set of normal samples.

  A HET SNP is generated if it exists in >66% of samples.
  A HOM SNP is generated if it exists in >33% of samples.

  Locations with more than 2 alleles expressed across the panel are excluded.
  Locations within 500bp of another potential SNP are excluded.

USAGE_END

die $usage if(@ARGV < 2);

my $ref = shift @ARGV;

my %files_by_index;
my @indexes;
for my $snp_file(@ARGV) {
  my ($idx) = $snp_file =~ m/[.]([[:digit:]]+)$/;
  push @indexes, $idx unless(exists $files_by_index{$idx});
  push @{$files_by_index{$idx}}, $snp_file;
}
@indexes = sort {$a<=>$b} @indexes;
my $total_indicies = scalar @indexes;
my $total_samples = scalar @{$files_by_index{1}};

my $min_het_samples = int $total_samples * 0.66;
my $min_hom_samples = int $total_samples * 0.33;

my %loc;
for my $idx(@indexes) {
  warn "Processing index $idx/$total_indicies\n";
  my $total_samples = scalar @{$files_by_index{$idx}};
  for my $file(@{$files_by_index{$idx}}) {
    my ($sample) = $file =~ m/^(.+)[.][[:digit:]]+$/;
    open my $IN, '<', $file or die $!;
    while(my $line = <$IN>) {
      next if($. == 1);
      chomp $line;
      my ($id, $chr, $pos, $alleleA, $fracA, $alleleB, $fracB, $used_depth, $total_depth) = split /\t/, $line;
      my $abAllele;
      if($alleleB eq q{-}) {
        $abAllele = $alleleB;
      }
      else {
        $abAllele = join(q{}, sort($alleleA,$alleleB));
      }
      push @{$loc{$chr}{$pos}{$abAllele}}, $sample;
    }
    close $IN;
  }
}

my $index = Bio::DB::HTS::Faidx->new($ref);
my @seq_ids = $index->get_all_sequence_ids();
print "\tChr\tPosition\n";

my $het_id = 1;
my $hom_id = 1;
for my $seq_id(@seq_ids) {
  my @positions = sort {$a<=>$b} keys %{$loc{$seq_id}};
  my $max = (scalar @positions)-1;
  for my $idx (0..$max) {
    my $pos = $positions[$idx];
    my @abAlleles = keys %{$loc{$seq_id}{$pos}};
    if(scalar @abAlleles > 1) {
      if(@abAlleles == 2 && exists $loc{$seq_id}{$pos}{'-'}) {
        # could be hom/het
        for my $al(@abAlleles) {
          next if($al eq q{-});
          push @{$loc{$seq_id}{$pos}{$al}}, @{$loc{$seq_id}{$pos}{'-'}};
          delete $loc{$seq_id}{$pos}{'-'};
          @abAlleles = ($al);
        }
      }
      else {
        # this skips positions with multiple combinations of alleles
        next;
      }
    }
    if($idx > 1) {
      next if($pos - 500 < $positions[$idx-1]);
    }
    elsif($idx+1 < $max) {
      next if($pos + 500 > $positions[$idx+1]);
    }
    my $ab = $abAlleles[0];
    my @samples = @{$loc{$seq_id}{$pos}{$ab}};
    if($ab eq q{-}) {
      next unless(@samples > $min_hom_samples);
      printf "HOM_%s\t%s\t%d\n", $hom_id++, $seq_id, $pos;
    }
    else {
      next unless(@samples > $min_het_samples);
      printf "HET_%s\t%s\t%d\n", $het_id++, $seq_id, $pos;
    }
  }
}
