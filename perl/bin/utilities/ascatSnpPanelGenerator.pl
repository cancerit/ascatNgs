#!/usr/bin/perl

########## LICENCE ##########
# Copyright (c) 2014-16 Genome Research Ltd.
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
########## LICENCE ##########

use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Pod::Usage;
use Bio::DB::HTS;
use Bio::DB::HTS::Tabix;
use Const::Fast qw(const);
use File::Path qw(make_path);

const my @ALLELES => (qw(A C G T));

# horrid globals used in callback, initialised just before it is executed.
my ($_ref_seq, $_chr, $_start, $_end, $_region, $OFH);
my $_row_id = 1;


my $options = option_builder();

my $hts = Bio::DB::HTS->new(-bam  => $options->{'hf'},
                            -fasta => $options->{'r'});

my $regions = get_regions($options, $hts);
my $outfile = prep_outfile($options, $hts);
find_snps($options, $hts, $outfile, $regions);

sub find_snps {
  my ($options, $hts, $outfile, $regions) = @_;

  my $fai = $hts->fai;

  open $OFH, '>', $outfile or die "Can't create $outfile: $!";
  print $OFH "\tChr\tPosition\tA\tAfrac\tB\tBfrac\tUsedDepth\tTotalDepth\n";

  for my $r(@{$regions}) {
    ($_region, $_chr, $_start, $_end) = @{$r};
    $_ref_seq = $fai->fetch($_region); # single request for whole region
    warn "Pileup of region: $_region\n";
    $hts->fast_pileup($_region, \&snp_search_cb);
  }

  close $OFH or die "Can't close $outfile: $!";
}

sub snp_search_cb {
  my ($seqid,$pos,$pileup) = @_;
  return unless($pos >= $_start && $pos <= $_end);

  my $total_depth = scalar @{$pileup};
  return if($total_depth < 15 || $total_depth > 50);

  my %calls = ( 'A' => 0,
                'C' => 0,
                'G' => 0,
                'T' => 0);
  for my $p(@{$pileup}) {
    # order of exclusion tests has been optimised
    # to ensure largest filter happens at each step
    # current order gave ~15% speedup

    next if( $p->indel != 0   # skip indel positions
          || $p->is_del       # skip padded del
          || $p->is_refskip); # skip ref-skip

    my $a = $p->alignment;
    next if($a->qual < 35); # only want high qual mappings
    my $flag = $a->flag;
    next if( ($flag & 3) != 3         # must be paired, mapped properly
          || ($flag & 3852) == 3852); # either unmapped, non-primary, QC-fail, Dup or secondary

    my $qpos = $p->qpos;
    my $qbase  = substr($a->qseq,$qpos,1);
    next unless($qbase =~ m/^[ACTGacgt]/);
    next if($a->qscore->[$qpos] < 20);
    $calls{uc $qbase}++;
  }

  my $depth = $calls{A} + $calls{C} + $calls{G} + $calls{T};
  return if($depth/$total_depth < 0.8);

  my @freqs;
  for my $b(@ALLELES) {
    next if($calls{$b} == 0);
    push @freqs, [$b, $calls{$b} / $depth];
  }
  return if((scalar @freqs) == 4);

  @freqs = reverse sort { $a->[1] <=> $b->[1] } @freqs;
  if($freqs[0][1] >= 0.95) {
    return if(substr($_ref_seq, $pos - $_start, 1) eq $freqs[0][0]);

    push @freqs, [q{-},0];
  }
  else {
    return if($freqs[0][1] + $freqs[1][1] <= 0.9);
    return if($freqs[1][1] < 0.25);
  }
  printf $OFH "CN_%d\t%s\t%d\t%s\t%.2f\t%s\t%.2f\t%d\t%d\n", $_row_id++, $_chr, $pos, $freqs[0][0], $freqs[0][1], $freqs[1][0], $freqs[1][1], $depth, $total_depth;
}

sub get_regions {
  my ($options, $hts) = @_;
  my @targets = $hts->seq_ids;

  my @top_ranges;
  if(defined $options->{'c'}) {
    open my $BED, '<', $options->{'c'} or die $!;
    while(my $line = <$BED>) {
      my @bits = split /\t/, $line;
      $bits[1]++;
      if(defined $options->{'i'}) {
        if($. == $options->{'i'}) {
          push @top_ranges, \@bits;
          last;
        }
      }
      else {
        push @top_ranges, \@bits;
      }
    }
    close $BED;
  }
  else {
    for(@targets) {
      push @top_ranges, [$_, 1, $hts->length($_)];
    }
  }

  my @refined_regions;
  my $exclude;
  $exclude = Bio::DB::HTS::Tabix->new(filename => $options->{'e'}) if(defined $options->{'e'});
  for my $t_range(@top_ranges) {
    my ($c_chr, $c_start, $c_end) = @{$t_range};
    my $iter;
    $iter = $exclude->query(sprintf "%s:%d-%d", $c_chr, $c_start, $c_end) if(defined $exclude);
    while(my $record = $iter->next){
      my ($start0, $end1) = (split /\t/, $record)[1..2];
      push @refined_regions, [sprintf("%s:%d-%d", $c_chr, $c_start, $start0), $c_chr, $c_start, $start0] unless($start0 == 0);
      $c_start = $end1;
    }
    push @refined_regions, [sprintf("%s:%d-%d", $c_chr, $c_start, $c_end), $c_chr, $c_start, $c_end];
  }
  return \@refined_regions;
}

sub prep_outfile {
  my ($options, $hts) = @_;
  my @raw_samples = $hts->header->text =~ m/\@RG[^\n]*[\t]SM\:([^\t\n]+)/gxms;
  die "ERROR: Input file has no 'SM' field in '\@RG' header lines\n" if(scalar @raw_samples == 0);

  my @uniq_samples = keys { map {$_ => 1} @raw_samples };
  die "ERROR: Input file header indicates it contains multiple samples, '\@RG' header tag 'SM' values:\n\t".
        join(q{,},@uniq_samples)."\n" if(scalar @uniq_samples > 1);

  my $sample = $uniq_samples[0];
  $sample =~ s/[^[:alpha:][:digit:]_\-.]/_/g; # sanitise the name

  my $idx = 0;
  $idx = $options->{'i'} if(defined $options->{'i'});
  return sprintf '%s/%s.tsv.%d', $options->{'o'}, $sample, $idx;
}


sub option_builder {
	my ($factory) = @_;

	my %opts = ();

	GetOptions (
		'h|help' => \$opts{'h'},
		'hf|hts-file=s' => \$opts{'hf'},
		'e|exclude=s' => \$opts{'e'},
		'r|ref=s' => \$opts{'r'},
		'i|index=n' => \$opts{'i'},
		'c|chunks=s' => \$opts{'c'},
		'o|outdir=s' => \$opts{'o'},
	);

	pod2usage(0) if($opts{'h'});
	pod2usage('No BAM(s) supplied') unless($opts{'hf'});
	pod2usage('No reference supplied') unless($opts{'r'});
	pod2usage('No outdir supplied') unless($opts{'o'});

	make_path($opts{'o'}) unless(-d $opts{'o'});

	return \%opts;
}

__END__

=head1 NAME

ascatSnpPanelGeneration.pl - Identifies likely HET/HOM locations for an individual sample.

=head1 SYNOPSIS

Attempts to generate a panel of obvious SNPs for use in copynumber.

snpPanelGeneration.pl [-h]

  Required Options:

    -hts-file (-hf) Indexed BAM file
                     - These should be normals (or as close to normal as you can get)

    -ref      (-r)  Relevant reference file (with *.fai co-located)

    -outdir   (-o)  Output directory
                      Files will be named as per BAM header sample name e.g.
                        SAMPLE_hets.tsv.<index>

    -chunks   (-c)  Bed file of regions to search

    -index    (-i)  Line from chunks file to process
                      See ascatFaiChunk.pl

    -exclude  (-e)  Bed file of regions to exclude (e.g. simple repeats)

    -help     (-h)  Brief documentation


  Examples:
    One job to generate complete HET candidates for sample (not recommended):
      ascatSnpPanelGeneration.pl -ref genome.fa -hf sampleA.bam > sampleA-hets.tsv.0

    One job per line of chunk file (i.e. parallel):
      ascatSnpPanelGeneration.pl -ref genome.fa -hf sampleA.bam -c chunks.bed -i 1 > sampleA-hets.tsv.1

=cut
