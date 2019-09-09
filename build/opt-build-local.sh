#! /bin/bash

set -xe

if [[ -z "${TMPDIR}" ]]; then
  TMPDIR=/tmp
fi

set -u

if [ "$#" -lt "1" ] ; then
  echo "Please provide an installation path such as /opt/ICGC"
  exit 1
fi

# get path to this script
SCRIPT_PATH=`dirname $0`;
SCRIPT_PATH=`(cd $SCRIPT_PATH && pwd)`

# get the location to install to
INST_PATH=$1
mkdir -p $1
INST_PATH=`(cd $1 && pwd)`
echo $INST_PATH

# get current directory
INIT_DIR=`pwd`

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
echo "Max compilation CPUs set to $CPU"

SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR/distro # don't delete the actual distro directory until the very end
mkdir -p $INST_PATH/bin
cd $SETUP_DIR

# make sure tools installed can see the install loc of libraries
set +u
export LD_LIBRARY_PATH=`echo $INST_PATH/lib:$LD_LIBRARY_PATH | perl -pe 's/:\$//;'`
export PATH=`echo $INST_PATH/bin:$PATH | perl -pe 's/:\$//;'`
export MANPATH=`echo $INST_PATH/man:$INST_PATH/share/man:$MANPATH | perl -pe 's/:\$//;'`
export PERL5LIB=`echo $INST_PATH/lib/perl5:$PERL5LIB | perl -pe 's/:\$//;'`
set -u

# Rlibs
if [ ! -e $SETUP_DIR/Rlib.success ]; then
  cd $INIT_DIR/Rsupport
  ./setupR.sh $OPT
  touch $SETUP_DIR/Rlib.success
fi

#add bin path for install tests
export PATH="$INST_PATH/bin:$PATH"

cd $INIT_DIR/perl

echo "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
cpanm --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null

echo "Installing ascatngs (perl)..."
cd $INIT_DIR/perl
perl perl/Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install

# cleanup all junk
rm -rf $SETUP_DIR
