# --- 1. R Packages
# --- 1.1. Tidy environment-related
library(tidyverse)
library(tidymodels)
library(DALEX)
library(DALEXtra)
library(parallel)
library(abind) 

# --- 1.2. Data table opening and storage
library(RSQLite)
library(feather) 
library(vroom)

# --- 1.2. Spatial data and object
library(raster)
library(virtualspecies) 
library(ncdf4)

# --- 1.3. Others
library(RColorBrewer) 
library(pastecs)
library(ecospat)
library(reticulate) 

# --- Seed
set.seed(123)

# --- Input / Output directories
project_wd <- getwd()

# --- Necessary code steps
source(file = "./code/01a_list_bio.R")
source(file = "./code/01b_folder_init.R")
source(file = "./code/02a_query_bio.R")
source(file = "./code/02b_query_env.R")
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
LOCAL_CLUSTERS <- 16
N_CLUSTERS <- 16
