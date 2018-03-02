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

use strict;
use Test::More;
use Test::Fatal;
use Data::Dumper;
use File::Spec;

use FindBin qw($Bin);
use Sanger::CGP::Vcf::VcfUtil;
use IPC::Open3;

my $lib_path = "$Bin/../lib";
my $bin_path = "$Bin/../bin";

my $SCRIPT = "perl $Bin/../bin/ascatCnToVCF.pl";
my $test_data_path = "$Bin/../testData/";
my $outfile = 'output.vcf';

my $test_mut_bam = File::Spec->catfile($test_data_path,'TESTMUT.bam');
my $test_mut_bai = File::Spec->catfile($test_data_path,'TESTMUT.bam.bai');
my $test_norm_bam = File::Spec->catfile($test_data_path,'TESTNORM.bam');
my $test_norm_bai = File::Spec->catfile($test_data_path,'TESTNORM.bam.bai');
my $test_segments_good = File::Spec->catfile($test_data_path,'test_segments.txt');
my $test_segments_bad = File::Spec->catfile($test_data_path,'test_segments_bad.txt');
my $ref = File::Spec->catfile($test_data_path,'genome_22.fa');

#Setup all params...
my $mut_study = 6;
my $norm_study = 9;
my $mut_accession = 12345;
my $norm_accession = 54321;
my $mut_source = 'TUM_SRC';
my $norm_source = 'NORM_SRC';
my $mut_protocol = 'TUM_PROT';
my $norm_protocol = 'NORM_PROT';
my $mut_platform = 'HiSeq';
my $norm_platform = 'HiSeq';
my $mut_desc = 'TUM_SAMP';
my $norm_desc = 'NORM_SAMP';
my $ref_spp = 'HUMAN';
my $ref_ass = '37';

my $command = $SCRIPT." ";
  $command .= "-sbm $test_mut_bam ";
  $command .= "-sbw $test_norm_bam ";
  $command .= "-mss $mut_study ";
  $command .= "-wss $norm_study ";
  $command .= "-msa $mut_accession ";
  $command .= "-wsa $norm_accession ";
  $command .= "-msc $mut_source ";
  $command .= "-wsc $norm_source ";
  $command .= "-msq $mut_protocol ";
  $command .= "-wsq $norm_protocol ";
  $command .= "-msp $mut_platform ";
  $command .= "-wsp $norm_platform ";
  $command .= "-msd $mut_desc ";
  $command .= "-wsd $norm_desc ";
  $command .= "-rs $ref_spp ";
  $command .= "-ra $ref_ass ";
  $command .= "-r $ref ";

my $EXP_OUT = [ "##fileformat=VCFv4.1",
                "##fileDate=".Sanger::CGP::Vcf::VcfUtil->get_date(),
                "##source_".Sanger::CGP::Vcf::VcfUtil->get_date().".1=ascatCnToVCF.pl_v".Sanger::CGP::Vcf->VERSION,
                "##reference=$ref",
                "##contig=<ID=1,assembly=37,length=249250621,species=HUMAN>",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of this structural variant\">",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
                "##ALT=<ID=CNV,Description=\"Copy number variable region\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=TCN,Number=1,Type=Integer,Description=\"Total copy number\">",
                "##FORMAT=<ID=MCN,Number=1,Type=Integer,Description=\"Minor allele copy number\">",
                "##vcfProcessLog_".Sanger::CGP::Vcf::VcfUtil->get_date().".1=<InputVCFSource=<ascatCnToVCF.pl>,InputVCFVer=<".Sanger::CGP::Vcf->VERSION.">>",
                "##SAMPLE=<ID=NORMAL,Description=\"$norm_desc\",Accession=$norm_accession,Platform=$norm_platform,Protocol=$norm_protocol,SampleName=TESTNORM,Source=$norm_source,Study=$norm_study>",
                "##SAMPLE=<ID=TUMOUR,Description=\"$mut_desc\",Accession=$mut_accession,Platform=$mut_platform,Protocol=$mut_protocol,SampleName=TESTMUT,Source=$mut_source,Study=$mut_study>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOUR",
                "22\t10000009\t.\tN\t<CNV>\t.\t.\tSVTYPE=CNV;END=62896419\tGT:TCN:MCN\t./.:2:1\t./.:3:1",
                "22\t10000024\t.\tN\t<CNV>\t.\t.\tSVTYPE=CNV;END=63111528\tGT:TCN:MCN\t./.:2:1\t./.:5:2"];


subtest 'Good segments file test' => sub{
  my $newcommand = $command."-f $test_segments_good ";
  $newcommand .= "-o $outfile ";

  my $res = `$newcommand`;

  ok(!$res,'Check script run');

  my $exists = 0;
  $exists = 1 if(-e $outfile);
  ok($exists,'Check output exists');

  my $FH;
  my $result_out= "";
  open($FH,'<',$outfile) or die("Error opening $outfile to get test results: $!");
    while(<$FH>){
      $result_out .= $_;
    }
  close($FH);
  is_deeply(result_to_arrayref($result_out), $EXP_OUT, 'Compare output results for good segment file');

  unlink($outfile);
};

subtest 'Good segments STDOUT test' => sub{
  my $newcommand = $command."-f $test_segments_good ";

  my $res = `$newcommand`;

  is_deeply(result_to_arrayref($res),$EXP_OUT,'Compare output STDOUT results for good segment file');
};

subtest 'Good segments STDIN and STDOUT test' => sub {
  my $newcommand = "cat $test_segments_good | ".$command;

  my $res = `$newcommand`;

  is_deeply(result_to_arrayref($res),$EXP_OUT,'Compare output STDOUT results for good segment file');
};

subtest 'Bad segments check for error' => sub {
  my $newcommand = "cat $test_segments_bad | ".$command;

  #Using IPC::Open3 to capture STDERR
  my ($stdin,$stdout,$stderr);
  use Symbol 'gensym'; $stderr = gensym;
  my $pid = open3($stdin,$stdout,$stderr,$newcommand);
  waitpid( $pid, 0 );
  my $return = $? >> 8;

  my $cnt = 0;
  while(<$stderr>){
    $cnt++;
    my $line = $_;
    if($cnt == 1){
      like($stderr, qr/^Error! Reading from: |STDIN| Writing to: |STDOUT|/,'Check for error message');
    }elsif($cnt == 3){
      like($stderr, qr/^Unrecognised format in CN segments passed/,'Check for error message segments');
    }
  }
  isnt($return, 1,'Check for error return value');
};

done_testing();

sub result_to_arrayref {
  my @t = split "\n", shift;
  return \@t;
}
