#!/usr/bin/perl

##########LICENCE##########
# Copyright (c) 2014 Genome Research Ltd.
#
# Author: David Jones <cgpit@sanger.ac.uk>
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

use Bio::DB::Sam;
use Try::Tiny;

use Sanger::CGP::Vcf;
use Sanger::CGP::Vcf::VCFCNConverter;
use Sanger::CGP::Vcf::BamUtil;

use Sanger::CGP::Vcf;
use Sanger::CGP::Vcf::Contig;
use Sanger::CGP::Vcf::Sample;
use Sanger::CGP::Vcf::VcfProcessLog;

{
  my $opts = setup();

  my $mt_sam = Bio::DB::Sam->new(-bam => $opts->{'sbm'}, -fasta => $opts->{'r'});
  my $wt_sam = Bio::DB::Sam->new(-bam => $opts->{'sbw'}, -fasta => $opts->{'r'});
  my $fai = Bio::DB::Sam::Fai->load($opts->{'r'});

  #parse samples and contigs from the bam files.
  my $contigs = Sanger::CGP::Vcf::BamUtil->parse_contigs($mt_sam->header->text.$wt_sam->header->text,$opts->{'rs'},$opts->{'ra'});
  my $mt_samples = Sanger::CGP::Vcf::BamUtil->parse_samples($mt_sam->header->text,$opts->{'mss'},$opts->{'msq'},$opts->{'msa'},$opts->{'msc'},$opts->{'msd'});
  my $wt_samples = Sanger::CGP::Vcf::BamUtil->parse_samples($wt_sam->header->text,$opts->{'wss'},$opts->{'wsq'},$opts->{'wsa'},$opts->{'wsc'},$opts->{'wsd'});

  die "No samples found in normal bam file." if(scalar values %$wt_samples == 0);
  die "Multiple samples found in normal bam file." if(scalar values %$wt_samples > 1);
  die "No samples found in mutant bam file." if(scalar values %$mt_samples == 0);
  die "Multiple samples found in mutant bam file." if(scalar values %$mt_samples > 1);

  #Setup the converter with the parsed contigs
  my $record_converter = new Sanger::CGP::Vcf::VCFCNConverter(
    -contigs => [values %$contigs]
  );

  my ($input_loc,$output_loc,$IN_FH,$OUT_FH);
  try{
    #Open in and out files (OR STDIN/STDOUT depending on what's been given at command line)
    if($opts->{'f'}){
      open($IN_FH,'<',$opts->{'f'}) or die("Error opening input file $opts->{f}: $!");
      $input_loc = $opts->{'f'};
    }else{
      $IN_FH = \*STDIN;
      $input_loc = 'STDIN';
    }

    if($opts->{'o'}){
      open($OUT_FH, '>', $opts->{'o'}) or die("Error opening output file $opts->{o}: $!");
      $output_loc = $opts->{'o'};
    }else{
      $OUT_FH = \*STDOUT;
      $output_loc = 'STDOUT';
    }

    #Generate header.
    my @process_logs = ();
    my $wt_sample = (values(%$wt_samples))[0];
    my $mt_sample = (values(%$mt_samples))[0];
    push @process_logs, new Sanger::CGP::Vcf::VcfProcessLog(
	        -input_vcf_source => basename($0),
          -input_vcf_ver => Sanger::CGP::Vcf->VERSION,
          -input_vcf_param => $opts,
	  );
    my $source = basename($0). '_v'. Sanger::CGP::Vcf->VERSION;
    print $OUT_FH $record_converter->generate_header($wt_sample,$mt_sample,
                                   \@process_logs, $opts->{'r'}, $source)
                                          or croak("Unable to write VCF header: $!");


    #Iterate through input and create a record for each.
    while(<$IN_FH>){
      my $line = $_;
      chomp($line);
      my ($seg_no,$chr,$start,$end,$wt_cn_tot,$wt_cn_min,$mt_cn_tot,$mt_cn_min) = split(',',$line);
      die("Unrecognised format in CN segments passed: '$line'.") unless(defined($seg_no)
                                                                      && defined($chr)
                                                                      && defined($start)
                                                                      && defined($end)
                                                                      && defined($wt_cn_tot)
                                                                      && defined($wt_cn_min)
                                                                      && defined($mt_cn_tot)
                                                                      && defined($mt_cn_min));

      my $start_allele = $fai->fetch("$chr:$start-$start");
      print $OUT_FH $record_converter->generate_record($chr,$start,$end,$start_allele,$wt_cn_tot,$wt_cn_min,$mt_cn_tot,$mt_cn_min);
    }

  }catch{
    die "Error! Reading from: |$input_loc| Writing to: |$output_loc|\n$_";
  }finally{
    close($IN_FH) or die("Error trying to close $opts->{f}: $!") if($opts->{'f'});
    close($OUT_FH) or die("Error trying to close $opts->{o}: $!") if($opts->{'o'});
  }

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
					'o|out-file=s' => \$opts{'o'},
          'f|in-file=s' => \$opts{'f'},
          'mss|sample-study-mut=s' => \$opts{'mss'},
          'msa|sample-accession-mut=s' => \$opts{'msa'},
          'msc|sample-accession-source-mut=s' => \$opts{'msc'},
          'msp|sample-platform-mut=s' => \$opts{'msp'},
          'msq|sample-sequencing-protocol-mut=s' =>\$opts{'msq'},
          'msd|sample-description-mut=s' => \$opts{'msd'},
          'wss|sample-study-norm=s' => \$opts{'wss'},
          'wsa|sample-accession-norm=s' => \$opts{'wsa'},
          'wsc|sample-accession-source-norm=s' => \$opts{'wsc'},
          'wsp|sample-platform-norm=s' => \$opts{'wsp'},
          'wsq|sample-sequencing-protocol-norm=s' =>\$opts{'wsq'},
          'wsd|sample-description-norm=s' => \$opts{'wsd'},
          'sbm|sample-bam-mut=s' => \$opts{'sbm'},
          'sbw|sample-bam-norm=s' => \$opts{'sbw'},
          'rs|reference-species=s' => \$opts{'rs'},
          'ra|reference-assembly=s' => \$opts{'ra'},
          'r|reference=s' => \$opts{'r'},
          '<>' => sub{push(@random_args,shift(@_));}
  ) or pod2usage(2);

  if(defined $opts{'v'}){
    print "Version: ".Sanger::CGP::Vcf->VERSION."\n";
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});
  if($opts{'f'}){
    if(! -e $opts{'f'}){
      pod2usage(-message  => "\nERROR: f|in-file must an existing file if defined, otherwise STDIN is used: '".$opts{'f'}."'.\n", -verbose => 1,  -output => \*STDERR);
    }
  }
  pod2usage(-message  => "\nERROR: mss|sample-study-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'mss'});
  pod2usage(-message  => "\nERROR: msa|sample-accession-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'msa'});
  pod2usage(-message  => "\nERROR: msc|sample-accession-source-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'msc'});
  pod2usage(-message  => "\nERROR: msp|sample-platform-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'msp'});
  pod2usage(-message  => "\nERROR: msq|sample-sequencing-protocol-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'msq'});
  pod2usage(-message  => "\nERROR: msd|sample-description-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'msd'});
  pod2usage(-message  => "\nERROR: wss|sample-study-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'wss'});
  pod2usage(-message  => "\nERROR: wsa|sample-accession-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'wsa'});
  pod2usage(-message  => "\nERROR: wsc|sample-accession-source-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'wsc'});
  pod2usage(-message  => "\nERROR: wsp|sample-platform-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'wsp'});
  pod2usage(-message  => "\nERROR: wsq|sample-sequencing-protocol-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'wsq'});
  pod2usage(-message  => "\nERROR: wsd|sample-description-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'wsd'});
  pod2usage(-message  => "\nERROR: sbm|sample-bam-mut must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'sbm'});
  if(! -e $opts{'sbm'}){
    pod2usage(-message  => "\nERROR: sbm|sample-bam-mut must be an existing file.\n", -verbose => 1,  -output => \*STDERR);
  }
  pod2usage(-message  => "\nERROR: sbw|sample-bam-norm must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'sbw'});
  if(! -e $opts{'sbw'}){
    pod2usage(-message  => "\nERROR: sbw|sample-bam-norm must be an existing file.\n", -verbose => 1,  -output => \*STDERR);
  }
  pod2usage(-message  => "\nERROR: r|reference must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'r'});
  pod2usage(-message  => "\nERROR: rs|reference-species must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'rs'});
  pod2usage(-message  => "\nERROR: ra|reference-assembly must be defined.\n", -verbose => 1,  -output => \*STDERR) unless($opts{'ra'});
  return \%opts;
}


__END__

=head1 NAME

CN_to_VCF.pl - Take CGP copy number algorithm segment output and convert to VCF format.

=head1 SYNOPSIS

CN_to_VCF.pl [options]

  	Required parameters:
     -out-file                        -o   Output file.
     -in-file                         -f   Input file.
     -sample-study-mut                -mss  Sample study.
     -sample-accession-mut            -msa  Sample Accession.
     -sample-accession-source-mut     -msc  Sample Accession Source.
     -sample-platform-mut             -msp  Sample Platform.
     -sample-sequencing-protocol-mut  -msq  Sample Sequencing Protocol.
     -sample-description-mut          -msd  Sample Description.
     -sample-study-norm               -wss  Sample study.
     -sample-accession-norm           -wsa  Sample Accession.
     -sample-accession-source-norm    -wsc  Sample Accession Source.
     -sample-platform-norm            -wsp  Sample Platform.
     -sample-sequencing-protocol-norm -wsq  Sample Sequencing Protocol.
     -sample-description-norm         -wsd  Sample Description.
     -sample-bam-mut                  -sbm  Mutant sample bam file.
     -sample-bam-norm                 -sbw  Normal sample bam file.
     -reference-species               -rs   Reference species.
     -reference-assembly              -ra   Reference assembly.
     -reference                       -r    Reference file

    Other:
     -help     -h   Brief help message.
     -man      -m   Full documentation.
     -version  -v   Version information.


=head1 OPTIONS

=over 8

=item B<-out-file>

File to write VCF format output to. If not provided STDOUT is used.

=item B<-in-file>

File to read segment information from. If not provided STDIN is used.

=item B<-sample-study-mut>

Study of mutant sample

=item B<-sample-accession-mut>

Accession of mutant sample

=item B<-sample-accession-source-mut>

Source of mutant sample accession

=item B<-sample-platform-mut>

Platform used to sequence mutant sample

=item B<-sample-sequencing-protocol-mut>

Protocol used in sequencing mutant sample

=item B<-sample-description-mut>

Description of mutant sample

=item B<-sample-study-norm>

Study of normal sample

=item B<-sample-accession-norm>

Accession of normal sample

=item B<-sample-accession-source-norm>

Source of normal sample accession

=item B<-sample-platform-norm>

Platform used to sequence normal sample

=item B<-sample-sequencing-protocol-norm>

Protocol used in sequencing normal sample

=item B<-sample-description-norm>

Description of normal sample

=item B<-sample-bam-mut>

Mutant sample bam file

=item B<-sample-bam-norm>

Normal sample bam file

=item B<-reference-species>

Reference assembly

=item B<-reference-assembly>

Reference species

=item B<-reference>

Reference file

=item B<-version>

Print version information.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<CN_to_VCF.pl> will attempt convert the input segmented file into VCF copy number format.

=cut

