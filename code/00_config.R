# --- 1. System arguments
Sys.setenv(HDF5_USE_FILE_LOCKING="FALSE") # to be able to open .nc
Sys.setenv(RETICULATE_PYTHON = "/UP_home/aschickele/.virtualenvs/r-reticulate/bin/python")

# --- 2. R Packages
# --- 2.1. General use
if(!require("devtools")){install.packages("devtools")}

# --- 2.2. Tidy environment-related
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("tidymodels")){install.packages("tidymodels")}
if(!require("DALEX")){install.packages("DALEX")}
if(!require("DALEXtra")){install.packages("DALEXtra")}
if(!require("parallel")){install.packages("parallel")}
if(!require("caret")){install.packages("caret")}
if(!require("xgboost")){install.packages("xgboost")}

# --- 2.3. Data table opening and storage
if(!require("RSQLite")){install.packages("RSQLite")}
if(!require("feather")){install.packages("feather")}
if(!require("vroom")){install.packages("vroom")}

# --- 2.4. Spatial data and object
if(!require("raster")){install.packages("raster")}
if(!require("virtualspecies")){install.packages("virtualspecies")}
if(!require("ncdf4")){install.packages("ncdf4")}
if(!require("usdm")){install.packages("usdm")}

# --- 2.5. Data access service
if(!require("phyloseq")){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("phyloseq")}
if(!require("MGnifyR")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("MGnifyR")}
if(!require("robis")){install.packages("robis")}
if(!require("rgbif")){install.packages("rgbif")}

# --- 2.6. Others
if(!require("RColorBrewer")){install.packages("RColorBrewer")}
if(!require("fields")){install.packages("fields")}
if(!require("pastecs")){install.packages("pastecs")}
if(!require("ecospat")){install.packages("ecospat")}
if(!require("dendextend")){install.packages("dendextend")}
if(!require("mvrsquared")){install.packages("mvrsquared")}
if(!require("bestNormalize")){install.packages("bestNormalize")}

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
