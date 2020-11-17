# AscatNGS

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Master Badge][travis-master]][travis-base] | [![Develop Badge][travis-develop]][travis-base] |

AscatNGS contains the Cancer Genome Projects workflow implementation of the ASCAT copy number
algorithm for paired end sequencing.

**_We do not offer support for WXS analysis, please contact the authors of ASCAT [here][ascat-web]_**.

For details of the underlying algorithm please see the [ASCAT][ascat-web] site.  The GitHub
repository can be found [here][ascat-gh].

Contents:

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [ASCAT core algorithm](#ascat-core-algorithm)
- [Docker, Singularity and Dockstore](#docker-singularity-and-dockstore)
- [Dependencies/Install](#dependenciesinstall)
	- [Reference files](#reference-files)
- [Creating a release](#creating-a-release)
	- [Preparation](#preparation)
	- [Cutting the release](#cutting-the-release)
- [LICENCE](#licence)

<!-- /TOC -->

## ASCAT core algorithm

The core `ascat.R` script is pulled from the primary ASCAT repository during installation.  The
linked version is currently [`v2.5.1`][ascat-release].

### Minimum sequence depth

The minimum depth that you are likely to get reliable results from is 15x genome coverage.  If you
have very good quality data with good insert size (little or no overlap of reads) you may be successful
with 12x.

## Docker, Singularity and Dockstore

There is a pre-built image containing this codebase on quay.io.

* [dockstore-cgpwgs][ds-cgpwgs-git]: Contains additional tools for WGS analysis.

This was primarily designed for use with dockstore.org but can be used as normal containers.

The docker images are know to work correctly after import into a singularity image.

## Dependencies/Install

Please install the following first:

* [PCAP-core v2.0+][pcap-core-rel]
  * `ascatNgs` has other dependancies fulfilled by `PCAP-core`.
* [cgpVcf v2.0+][cgpvcf-rel]
* alleleCount
  * [alleleCount v4.0.0+][allelecount-rel] - for fragment based counting (recommended)
  * [alleleCount v3.0.0 - v3.3.1][allelecount-v3] - for read based counting
* [RColorBrewer][rcolorbrewer]

Please see these for any child dependencies.

Once complete please run:

```
./setup.sh /some/install/location
```

If you are installing each tool (`PCAP-code`, `cgpVcf`, `alleleCount`) to an independent area you should
set the environment variable `CGP_PERLLIBS` to include the relevant perl libraries from those, e.g.

```
export CGP_PERLLIBS=/pcapInst/libs/perl5:/cgpvcfInst/libs/perl5:...
```

Please be aware that this software requires the Rscript executable to be pre-installed.

### Reference files

Please see the [wiki][ascatngs-wiki] for how to obtain/generate the
`SnpGcCorrections.tsv` file.

## Creating a release

### Preparation

* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

### Cutting the release

1. Update `perl/lib/Sanger/CGP/Ascat.pm` to the correct version.
1. Update `CHANGES.md` with key updates.
1. Run `./prerelease.sh`
1. Check all tests and coverage reports are acceptable.
1. Commit the updated docs tree and updated module/version.
1. Push commits.
1. Use the GitHub tools to draft a release.

## LICENCE

```
Copyright (c) 2014-2020 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

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
```

<!-- Travis -->
[travis-base]: https://travis-ci.org/cancerit/ascatNgs
[travis-master]: https://travis-ci.org/cancerit/ascatNgs.svg?branch=master
[travis-develop]: https://travis-ci.org/cancerit/ascatNgs.svg?branch=dev

<!-- refs -->
[allelecount-v3]: https://github.com/cancerit/alleleCount/releases/tag/v3.3.1
[allelecount-rel]: https://github.com/cancerit/alleleCount/releases
[cgpvcf-rel]: https://github.com/cancerit/cgpVcf/releases
[pcap-core-rel]: https://github.com/cancerit/PCAP-core/releases
[ascat-web]: https://www.crick.ac.uk/peter-van-loo/software/ASCAT
[ascat-gh]: https://github.com/Crick-CancerGenomics/ascat
[ascat-release]: https://github.com/Crick-CancerGenomics/ascat/releases/tag/v2.5
[ascatngs-wiki]: https://github.com/cancerit/ascatNgs/wiki
[ds-cgpwgs-git]: https://github.com/cancerit/dockstore-cgpwgs
[rcolorbrewer]: https://cran.r-project.org/web/packages/RColorBrewer/index.html
