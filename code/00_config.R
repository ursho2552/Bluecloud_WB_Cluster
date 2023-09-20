# --- 1. System arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc from complex
Sys.setenv(RETICULATE_PYTHON = "/UP_home/aschickele/.virtualenvs/r-reticulate/bin/python")

# --- 2. R Packages
# --- 2.1. Tidy environment-related
library(tidyverse)
library(tidymodels)
library(DALEX)
library(DALEXtra)
library(parallel)
library(abind) 
library(caret)
library(xgboost)

# --- 2.2. Data table opening and storage
library(RSQLite)
library(feather) 
library(vroom)

# --- 2.2. Spatial data and object
library(raster)
library(virtualspecies) 
library(ncdf4)
library(usdm)

# --- 2.3. Data access service
library(phyloseq)
library(MGnifyR)
library(robis)
library(rgbif)

# --- 2.4. Others
library(RColorBrewer) 
library(fields)
library(pastecs)
library(ecospat)
library(reticulate) 
library(dendextend)
library(mvrsquared) 

# --- Seed
set.seed(123)

# --- Input / Output directories
project_wd <- getwd()

# --- Necessary code steps
source(file = "./code/01_list_bio_wrapper.R")
source(file = "./code/02_run_init.R")
source(file = "./code/03_query_bio_wrapper.R")
source(file = "./code/04_query_env.R")
source(file = "./code/05_pseudo_abs.R")
source(file = "./code/06_query_check.R")
source(file = "./code/07_folds.R")
source(file = "./code/08_hyperparameters.R")
source(file = "./code/09_model_wrapper.R")
source(file = "./code/10_eval_wrapper.R")
source(file = "./code/11_proj_wrapper.R")
source(file = "./code/12a_standard_maps.R")
source(file = "./code/12b_pdp.R")
source(file = "./code/12c_diversity_maps.R")
source(file = "./code/12d_user_synthesis.R")

# --- Custom functions
source("./function/nc_to_raster.R")
source("./function/sample_raster_NA.R")
source("./function/outlier_iqr_col.R")
source("./function/viridis.R")
source("./function/log_sink.R")
source("./function/QC_recommandations.R")
source("./function/bivar_raster_plot.R")
source("./function/get_cell_neighbors.R")
source("./function/regrid_env.R")

# --- Data specific parameters

# --- Model specific parameters
MAX_CLUSTERS <- 16
