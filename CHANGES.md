# Changes

## 4.4.1

* Update base image for performance from htslib 1.11.

## 4.4.0

* Added ascatCounts to produce counts files
* Modified ascat wrapper to handle count files as input

## 4.3.4

* Eliminated redundant logic from setup script
* Updated Dockerfile to use PCAP-core 5.2.1 as base (Ubuntu 20.04)

## 4.3.3

* Add Rprofile to fix issues with scientific notation

## 4.3.2

* Use lookup tables to translate contig names in case we've stripped out 'chr' in GRCh38
* Fix error in `*.samplestatistics.txt` - GenderChrFound was not being set to Y/N when found but rather the male contig name.

## 4.3.1

* Correct version in dockerfile
* Remove development libraries from final build
  * Drops image size from 1.33GB to ~1GB

## 4.3.0

* Added a Docker container

## 4.2.1

* Update tabix call to use query_full

## 4.2.0

* Upgrade to core R [ASCAT v2.5.1](https://github.com/Crick-CancerGenomics/ascat/releases/tag/v2.5)
  * applies patch to handle sparse data on a contig (usually Human Y)
* Recommend upgrading to [alleleCount v4.0.0+](https://github.com/cancerit/alleleCount/releases)
  * Switches to fragment based counts instead of read based to prevent double counting.

## 4.1.2

* Updated tabix->query to use tabix->query_full

## 4.1.1

* Fix error in setup.sh testing of alleleCounter version

## 4.1.0

* Modified to use new >= 3.3.0 allelecounter code with dense SNP functionality

## 2.1.2

* Minor fix to usage typo

## 2.0.0

* Now species agnostic provided you can generate a panel of SNPs, approx 1 HET per 2kb is required.
  * Tools are included in the package for guided and unguided SNP generation.
  * More details can be found in the [wiki](https://github.com/cancerit/ascatNgs/wiki).
* Supports CRAM
* Allele count phase is now threaded.

NOTE: Individual scripts called by `ascat.pl` have been renamed to prevent potential namespace clashes.
