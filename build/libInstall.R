instLib = commandArgs(T)[1]

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

ipak_bioc <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

if( (version$major == 3 && version$minor >=5) || version$major > 3) {
  # biocmanager versions of R
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(ask=FALSE, lib=instLib, lib.loc=instLib)
  ipak_bioc(c("data.table"))
  ipak_bioc(c("gam"))
  ipak_bioc(c("VGAM"))
  ipak_bioc(c("stringr"))
  ipak_bioc(c("mgcv"))
  ipak_bioc(c("poweRlaw"))
  ipak_bioc(c("zlibbioc"))
  ipak_bioc(c("RColorBrewer"))
} else {
  # OLD versions of R
  source("http://bioconductor.org/biocLite.R")
  ipak(c("data.table"))
  ipak(c("gam"))
  if ( version$major == 3 && version$minor < 2 ) {
    install.packages("VGAM_1.0-3.tar.gz", type="source", lib=instLib, lib.loc=instLib)
  } else {
    ipak(c("VGAM"))
  }
  ipak(c("stringr"))
  ipak(c("mgcv"))
  ipak(c("poweRlaw"))
  ipak(c("zlibbioc"))
  ipak(c("RColorBrewer"))
}

# works on old and new
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", lib=instLib)
library(devtools)
options(download.file.method = "auto")
install_github("Irrationone/copynumber", ref="87d2663fe6b11c03cf6006b4ee9ed70450eacb5a")
