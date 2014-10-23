LICENCE
=======
Copyright (c) 2014 Genome Research Ltd.

Author: CancerIT <cgpit@sanger.ac.uk>

This file is part of AscatNGS.

AscatNGS is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."

AscatNGS
========
AscatNGS contains the Cancer Genome Projects workflow implementation of the ASCAT copy number
algorithm for paired end sequencing.

For details of the underlying algorithm please see the [ASCAT](http://heim.ifi.uio.no/bioinf/Projects/ASCAT/) site.

---

###Dependencies/Install
Please install the following first:

* [PCAP-core](http://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases)
* [alleleCount](https://github.com/cancerit/alleleCount/releases)

Please see these for any child dependencies.

Once complete please run:

./setup.sh /some/install/location

Please be aware that this software requires the Rscript executable to be pre-installed.

---

##Creating a release
####Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
1. Update `perl/lib/Sanger/CGP/Ascat.pm` to the correct version (adding rc/beta to end if applicable).
2. Run `./prerelease.sh`
3. Check all tests and coverage reports are acceptable.
4. Commit the updated docs tree and updated module/version.
5. Push commits.
6. Use the GitHub tools to draft a release.
