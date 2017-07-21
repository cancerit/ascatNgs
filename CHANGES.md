### 4.1.0
* Modified to use new >= 3.3.0 allelecounter code with dense SNP functionality

### 2.1.2
* Minor fix to usage typo

### 2.0.0
* Now species agnostic provided you can generate a panel of SNPs, approx 1 HET per 2kb is required.  Tools are included in the package for guided and unguided SNP generation.  More details can be found in the [wiki](https://github.com/cancerit/ascatNgs/wiki).
* Supports CRAM
* Allele count phase is now threaded.

NOTE: Individual scripts called by `ascat.pl` have been renamed to prevent potential namespace clashes.
