### master.R ---
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  5 2026 (10:09)
## Version:
## Last-Updated: Sep  4 2026 (13:58)
##           By: Helene
##     Update #: 24
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

run.master <- function(install.simevent = FALSE, file.path = "") {
  library(data.table)
  library(survival)
  library(ggplot2)
  library(xtable)
  library(gridExtra)
  library(zoo)
  library(nleqslv)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(devtools)
  if (install.simevent) {
    remotes::install_github("miclukacova/simevent")
  }
  library(simevent)

  devtools::load_all(paste0("./", file.path))
}

######################################################################
### master.R ends here
