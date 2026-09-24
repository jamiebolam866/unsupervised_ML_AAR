# Mission: unsupervised clustering analyses

# check and set working directory
getwd()

# call renv::restore() to install packages from lockfile
renv::restore()

# set up renv
renv::init()
# install core setup packages
renv::install("tidyverse")
renv::install("here")
renv::install("usethis")
renv::install("pak")
renv::install("devtools")
renv::install("ddauber/r4np")

# machine learning
renv::install("YuHuiDeakin/rabc")

# load packages
library(tidyverse)
library(here)
library(usethis)
library(pak)
library(devtools)
library(r4np)

# machine learning
library(rabc)


# save packages to lockfile
renv::snapshot()

# ensure variables don't save between sessions
usethis::use_blank_slate()

# set up folder system using r4np package - if not already
# ensure r4np is loaded before running the below code (hashed out)
# r4np::create_project_folder()

# set up git
usethis::use_git()
usethis::use_github() # use this code to create GitHub repo, but I already did this before manually
