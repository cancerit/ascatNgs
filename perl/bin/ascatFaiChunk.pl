#!/usr/bin/perl

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
