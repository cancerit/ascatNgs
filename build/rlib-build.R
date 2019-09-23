instLib = commandArgs(T)[1]

r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)
source("http://bioconductor.org/biocLite.R")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    biocLite(new.pkg, ask=FALSE, lib=instLib, lib.loc=instLib)
  sapply(pkg, library, character.only = TRUE)
}

tmp <- c("devtools")
ipak(tmp)
library(devtools)
options(download.file.method = "auto")

# ASCAT and BRASS
ipak(c("data.table"))
ipak(c("mgcv"))
ipak(c("gam"))
ipak(c("VGAM"))
ipak(c("stringr"))
ipak(c("poweRlaw"))
ipak(c("zlibbioc"))
ipak(c("RColorBrewer"))

install_github("Irrationone/copynumber", ref="87d2663fe6b11c03cf6006b4ee9ed70450eacb5a")
