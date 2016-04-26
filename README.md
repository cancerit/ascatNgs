AscatNGS
========
AscatNGS contains the Cancer Genome Projects workflow implementation of the ASCAT copy number
algorithm for paired end sequencing.

For details of the underlying algorithm please see the [ASCAT](http://heim.ifi.uio.no/bioinf/Projects/ASCAT/) site.

---

###Dependencies/Install
Please install the following first:

* [PCAP-core v2.0+](http://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases)
  * `ascatNgs` has other dependancies fulfilled by `PCAP-core`.
* [cgpVcf v2.0+](https://github.com/cancerit/cgpVcf/releases)
* [alleleCount v3.0+](https://github.com/cancerit/alleleCount/releases)

Please see these for any child dependencies.

Once complete please run:

./setup.sh /some/install/location

If you are installing each tool (PCAP-code, cgpVcf, alleleCount) to an independent area you should set the environment variable `CGP_PERLLIB` to include the relevant perl libraries from those.

Please be aware that this software requires the Rscript executable to be pre-installed.

---

##Creating a release
####Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
1. Update `perl/lib/Sanger/CGP/Ascat.pm` to the correct version.
2. Update `CHANGES.md` with key updates.
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.

---

LICENCE
========
Copyright (c) 2014-2016 Genome Research Ltd.

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