#!/bin/bash

########## LICENCE ##########
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
########## LICENCE ##########

# v2.4.3
ASCAT_SRC="https://raw.githubusercontent.com/Crick-CancerGenomics/ascat/6d40e69a2919ddfc1cda870310203c772bf846ce/ASCAT/R/ascat.R"

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/ICGC"
  exit 0
fi

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo; echo
) >>$INIT_DIR/setup.log 2>&1

set -eu

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
PERLROOT=$INST_PATH/lib/perl5

# allows user to knowingly specify other PERL5LIB areas.
if [ -z ${CGP_PERLLIBS+x} ]; then
  export PERL5LIB="$PERLROOT"
else
  export PERL5LIB="$PERLROOT:$CGP_PERLLIBS"
fi

#add bin path for install tests
export PATH=$INST_PATH/bin:$PATH

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

cd $SETUP_DIR

## grab cpanm and stick in workspace, then do a self upgrade into bin:
get_file $SETUP_DIR/cpanm https://cpanmin.us/
perl $SETUP_DIR/cpanm -l $INST_PATH App::cpanminus
CPANM=`which cpanm`
echo $CPANM

PCAP=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' PCAP`
if [[ "x$PCAP" == "x" ]] ; then
  echo "PREREQUISITE: Please install PCAP-core before proceeding:"
  echo "  https://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases"
  exit 1;
fi

AC=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::AlleleCount`
if [[ "x$AC" == "x" ]] ; then
  echo "PREREQUISITE: Please install alleleCount before proceeding:"
  echo "  https://github.com/cancerit/alleleCount/releases"
  exit 1;
fi

VCF=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Sanger::CGP::Vcf`
if [[ "x$VCF" == "x" ]] ; then
  echo "PREREQUISITE: Please install cgpVcf before proceeding:"
  echo "  https://github.com/cancerit/cgpVcf/releases"
  exit 1;
fi

HTS=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' Bio::DB::HTS`
if [[ "x$HTS" == "x" ]] ; then
  echo "PREREQUISITE: Please install Bio::DB::HTS before proceeding:"
  echo "  https://github.com/Ensembl/Bio-HTS/releases"
  exit 1;
fi

perlmods=( "File::ShareDir" "File::ShareDir::Install" )

for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  $CPANM -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
  echo
  done_message "" "Failed during installation of $i."
done

cd $INIT_DIR/perl

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi

$CPANM -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null
done_message "" "Failed during installation of core dependencies."

echo -n "Installing ascatNgs ..."
get_file share/ascat/ascat.R $ASCAT_SRC && \
perl Makefile.PL INSTALL_BASE=$INST_PATH && \
make && \
make test && \
make install && \
rm share/ascat/ascat.R
done_message "" "ascatNgs install failed."

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

echo "!!! If you have pre-v1.8 installations of ascatNgs you may want to remove the following:";
echo -e "\t$INST_PATH/bin/CN_to_VCF.pl"
echo -e "\t$INST_PATH/bin/failed_cn_csv.pl"
echo "  These have been renamed to the following:";
echo -e "\t$INST_PATH/bin/ascatCnToVCF.pl"
echo -e "\t$INST_PATH/bin/ascatFailedCnCsv.pl"
echo


exit 0
