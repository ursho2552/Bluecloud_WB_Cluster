# --- R Package
library(tidyverse)
library(RSQLite)
# library(RPostgreSQL) not working, need to check that !
library(abind) 
library(feather) 
library(reticulate) 
library(parallel)

library(raster)
library(virtualspecies) 
library(ncdf4)
library(vroom)

library(RColorBrewer) 
library(pastecs)


# --- Input / Output directories
project_wd <- "/net/meso/work/aschickele/Diversity"
setwd(project_wd)

# --- Custom functions
source("./function/sample_raster_NA.R")
source("./function/outlier_iqr_col.R")

# --- Other custom arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex

# --- Data specific parameters

# --- Model specific parameters