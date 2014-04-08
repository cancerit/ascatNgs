#!/usr/bin/env perl

# Author: kr2 (Original)

### SVN INFO ###
# $LastChangedBy: kr2 $
################

use strict;
use Data::Dumper;
use Carp;
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use Getopt::Long 'GetOptions';
use Pod::Usage;

{
  my $options = option_builder();
  run($options);
}

sub run {
  my ($options) = @_;
  open my $FH, '>', $options->{'o'} or croak 'Failed to create '.$options->{'o'};
  my $geno_ob = Sanger::CGP::Ascat::Genotype->new();
  if($options->{'l'}) {
    $geno_ob->get_full_loci_profile($options->{'b'}, $FH, $options->{'l'}, $options->{'m'});
  }
  else {
    $geno_ob->get_full_snp6_profile($options->{'b'}, $FH, $options->{'m'});
  }
  close $FH;
}

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
		'h'    => \$opts{'h'},
		'b=s' => \$opts{'b'},
		'o=s' => \$opts{'o'},
		'l=s' => \$opts{'l'},
		'm=n' => \$opts{'m'},
	);

	pod2usage(0) if($opts{'h'});
	pod2usage(1) if(!$opts{'o'} || !$opts{'b'});

	return \%opts;
}

package Sanger::CGP::Ascat::Genotype;

use strict;
use Carp;
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use Bio::DB::Sam;
use Bio::DB::Bam::AlignWrapper;

use constant MAX_PILEUP_DEPTH => 1_000_000;
use constant MIN_MAPQ => 35;
use constant TAB => "\t";
use constant NL => "\n";

my $g_pu_data; # required for pileup;
my $g_pb_qual;
my $g_sam;

sub new {
  my ($class, $opts) = @_;
  my $self = { };
  bless $self, $class;
  if(defined $opts) {
    $self->{'species'} = $opts->{'species'};
    $self->{'build'} = $opts->{'build'};
  }
  return $self;
}

=item get_full_snp6_profile
  Writes tab seperated allelic counts and depth to specified FH
  Uses all snps defined in file used by ngs_cn (format slightly different)
=cut
sub get_full_snp6_profile {
  my ($self, $bam_file, $fh, $min_qual) = @_;
  $g_pb_qual = $min_qual || MIN_MAPQ;
  my $sam = Bio::DB::Sam->new(-bam => $bam_file);
  $sam->max_pileup_cnt(MAX_PILEUP_DEPTH);
  $g_sam = $sam;
  my $snp6_file = $self->ngs_cn_snps({'species'=>'HUMAN','build'=>37});
  my ($region, $chr, $pos, $allA, $allB);
  print $fh "#CHR\tPOS\tCount_Allele_A\tCount_Allele_B\tGood_depth\n" or croak "Failed to write line: $OS_ERROR\n";
  open my $SNP6, '<', $snp6_file or croak "Unable to open $snp6_file for reading: $OS_ERROR\n";
  while(my $line = <$SNP6>) {
    chomp $line;
    ($chr, $pos, undef, undef, $allA, $allB) = split /\s/, $line;
    $g_pu_data = Sanger::CGP::Ascat::Genotype::PileupData->new($chr, $pos, $allA, $allB);
    $region = $chr.':'.$pos.'-'.$pos;
    $sam->fast_pileup($region, \&allele_counts);
    print $fh $g_pu_data->chr,TAB,$g_pu_data->pos,TAB,$g_pu_data->count_A,TAB,$g_pu_data->count_B,TAB,$g_pu_data->depth,NL or croak "Failed to write line: $OS_ERROR\n";
  }
  close $SNP6;
  return 1;
}

=item get_full_loci_profile
  Writes tab seperated allelic counts and depth to specified FH
  Uses all loci defined in specified file
=cut
sub get_full_loci_profile {
  my ($self, $bam_file, $fh, $loci_file, $min_qual) = @_;
  $g_pb_qual = $min_qual || MIN_MAPQ;
  my $sam = Bio::DB::Sam->new(-bam => $bam_file);
  $sam->max_pileup_cnt(MAX_PILEUP_DEPTH);
  $g_sam = $sam;
  my ($region, $chr, $pos, $allA, $allB);
  print $fh "#CHR\tPOS\tCount_A\tCount_C\tCount_G\tCount_T\tGood_depth\tInputFile\n" or croak "Failed to write line: $OS_ERROR\n";
  open my $LOCI, '<', $loci_file or croak 'Unable to open '.$loci_file.' for reading';
  while(my $line = <$LOCI>) {
    chomp $line;
    ($chr, $pos) = split /\s/, $line;
    $g_pu_data = Sanger::CGP::Ascat::Genotype::PileupData->new($chr, $pos);
    $region = $chr.':'.$pos.'-'.$pos;
    $sam->fast_pileup($region, \&allele_counts);
    print $fh $g_pu_data->chr or croak "Failed to write line: $OS_ERROR\n";
    print $fh TAB,$g_pu_data->pos or croak "Failed to write line: $OS_ERROR\n";
    print $fh TAB,$g_pu_data->residue_count('A') or croak "Failed to write line: $OS_ERROR\n";
    print $fh TAB,$g_pu_data->residue_count('C') or croak "Failed to write line: $OS_ERROR\n";
    print $fh TAB,$g_pu_data->residue_count('G') or croak "Failed to write line: $OS_ERROR\n";
    print $fh TAB,$g_pu_data->residue_count('T') or croak "Failed to write line: $OS_ERROR\n";
    print $fh TAB,$g_pu_data->depth or croak "Failed to write line: $OS_ERROR\n";
    print $fh TAB,$bam_file,NL or croak "Failed to write line: $OS_ERROR\n";
  }
  close $LOCI;
  return 1;
}

sub allele_counts {
  my ($seqid, $pos, $pu) = @_;
  return if($pos != $g_pu_data->pos);
  foreach my $p (@{$pu}) {
    next if($p->indel || $p->is_refskip);
    my $a = $p->alignment;
    if($g_pu_data->readgrp) {
      next if($a->get_tag_values('RG') != $g_pu_data->readgrp);
    }
    next if(!$a->proper_pair);
    next if($a->unmapped || $a->munmapped);
    next if($a->qual < $g_pb_qual);

    if($g_pb_qual) {
      my $fa = Bio::DB::Bam::AlignWrapper->new($a, $g_sam);
      my $qual = ($fa->qscore)[$p->qpos];
      next if($qual <= $g_pb_qual);
    }

    # get the base at this pos
    my $qbase  = substr($a->qseq, $p->qpos, 1);
    $g_pu_data->register_allele($qbase);
  }
  return 1;
}


1;

package Sanger::CGP::Ascat::Genotype::PileupData;

use strict;
use Carp;
use English qw( -no_match_vars );
use warnings FATAL => 'all';

use Const::Fast qw(const);

const my $MIN_DEPTH => 4;

const my $BIT_COVERED => 1;
const my $BIT_ALLELE_A => 2;
const my $BIT_ALLELE_B => 4;

sub new {
  my ($class, $chr, $pos, $allA, $allB, $rg) = @_;
  my $self =  { 'chr' => $chr,
                'pos' => $pos,
                'allele_A' => $allA,
                'allele_B' => $allB,
                'count_A' => 0,
                'count_B' => 0,
                'depth' => 0,
                'readgroup_id' => $rg,
                'A' => 0,
                'C' => 0,
                'G' => 0,
                'T' => 0,
              };
  bless $self, $class;
  return $self;
}

sub chr {
  return shift->{'chr'};
}

sub pos {
  return shift->{'pos'};
}

sub residue_count {
  my ($self, $residue) = @_;
  return $self->{uc $residue};
}

sub readgrp {
  return shift->{'readgroup_id'};
}

sub allele_A {
  my ($self, $allele) = @_;
  $self->{'allele_A'} = $allele if($allele);
  return $self->{'allele_A'};
}

sub allele_B {
  my ($self, $allele) = @_;
  $self->{'allele_A'} = $allele if($allele);
  return $self->{'allele_A'};
}

sub inc_A {
  shift->{'count_A'}++;
  return 1;
}

sub inc_B {
  shift->{'count_B'}++;
  return 1;
}

sub count_A {
  return shift->{'count_A'};
}

sub count_B {
  return shift->{'count_B'};
}

sub depth {
  return shift->{'depth'}
}

sub register_allele {
  my ($self, $allele) = @_;
  $allele = uc $allele;
  if($self->{'allele_A'} && $allele eq $self->{'allele_A'}) {
    $self->inc_A;
  }
  elsif($self->{'allele_B'} && $allele eq $self->{'allele_B'}) {
    $self->inc_B;
  }
  $self->{$allele}++;
  $self->{'depth'}++;
  return 1;
}

sub encoded_snp_status {
  my $self = shift;
  my $bit_val = 0;
  if($self->depth >= $MIN_DEPTH) {
    $bit_val += $BIT_COVERED;
    $bit_val += $BIT_ALLELE_A if($self->count_A > 0);
    $bit_val += $BIT_ALLELE_B if($self->count_B > 0);
  }
  return $bit_val;
}

sub readable_snp_status {
  my $self = shift;
  my $status = q{};
  if($self->depth >= $MIN_DEPTH) {
    $status .= 'cov';
    $status .= 'A' if($self->count_A > 0);
    $status .= 'B' if($self->count_B > 0);
  }
  else {
    $status = 'nc';
  }
  return $status;
}

1;


__END__

=head1 NAME

snp6alleleCounts.pl - Generate tab seperated file with allelic counts and depth for each snp6 position.

=head1 SYNOPSIS

snp6alleleCounts.pl

  Required:

    -b      BWA bam file (expects co-located index)
    -o      Output file
    -m      Minimum base quality to include (integer)
    -l      Alternate loci file (just needs chr pos)
             - output is different, counts for each residue

  Optional:
    -h        This message

=cut

