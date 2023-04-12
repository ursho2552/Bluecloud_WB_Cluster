# --- R Package
library(tidyverse)
library(tidymodels)
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

library(ecospat)

# --- Seed
set.seed(123)

# --- Input / Output directories
project_wd <- getwd()

# --- Necessary code steps
source(file = "./code/01a_list_bio.R")
source(file = "./code/01b_query_bio.R")
source(file = "./code/02_query_env.R")
source(file = "./code/03_query_check.R")
source(file = "./code/04_pseudo_abs.R")
source(file = "./code/05_folds.R")
source(file = "./code/06_hyperparameters.R")
source(file = "./code/07a_model_wrapper.R")
source(file = "./code/08a_eval_wrapper.R")
source(file = "./code/09a_proj_wrapper.R")
source(file = "./code/10a_standard_maps.R")

# --- Custom functions
source("./function/sample_raster_NA.R")
source("./function/outlier_iqr_col.R")
source("./function/viridis.R")

# --- Other custom arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex

# --- Data specific parameters

# --- Model specific parameters