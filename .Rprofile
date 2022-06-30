#Sys.setenv(TMPDIR="/home/ubuntu/mnt/tmp")
Sys.setenv(MAKEFLAGS = "-j4")
source("renv/activate.R")
renv::settings$use.cache(TRUE)
options("readr.num_threads" = 1)

# options(renv.config.pak.enabled=TRUE)

## For Linux and Windows users, we'll use RStudio Package Manager (RSPM).
#if (Sys.info()[['sysname']] %in% c('Linux', 'Windows')) {
#  options(repos = c(RSPM = "https://packagemanager.rstudio.com/all/latest"))
#} else {
  ## For Mac users, we'll default to installing from CRAN/MRAN instead, since
  ## RSPM does not yet support Mac binaries.
#  options(repos = c(CRAN = "https://cran.rstudio.com/"))
  # options(renv.config.mran.enabled = TRUE) ## TRUE by default
#}
#options(renv.config.repos.override = getOption("repos"))

