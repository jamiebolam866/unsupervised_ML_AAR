# Mission: unsupervised clustering analyses

# check and set working directory
getwd()

# call renv::restore() to install packages from lockfile
renv::restore()

# set up renv
renv::init()
# install packages
renv::install("tidyverse")
renv::install("here")
renv::install("usethis")
renv::install("pak")
renv::install("devtools")
pak::pkg_install("ddauber/r4np")

# movement & time-series data analysis
renv::install("dygraphs")
renv::install("xts")
renv::install("tidyr")
# machine learning
renv::install("umap")
pak::pak("YuHuiDeakin/rabc")
#plotting
renv::install("RColorBrewer")
renv::install("htmlwidgets")

# load packages
library(tidyverse)
library(here)
library(usethis)
library(pak)
library(devtools)
library(stringi)
library(r4np)
library(data.table)
# movement & time-series data analysis
library(dygraphs)
library(xts)
library(tidyr)
# machine learning
library(umap)
library(rabc)
# plotting
library(RColorBrewer)
library(htmlwidgets)

# save packages to lockfile
renv::snapshot()

# ensure variables don't save between sessions
usethis::use_blank_slate()

# set up folder system using r4np package - if not already
# ensure r4np is loaded before running the below code (hashed out)
# r4np::create_project_folder()

# set up git
# usethis::use_git()
# usethis::use_github() # use this code to create GitHub repo, but I already did this before manually
