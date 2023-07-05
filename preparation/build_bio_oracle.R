#' =============================================================================
#' @name build_bio_oracle
#' @description downloads the surface data available within Bio-ORACLE. Builds
#' a raster stack with the data and re-grid to 1x1 degree. The downloaded data
#' are accessible in the /data folder for further use in the pipeline. This 
#' script does not need to be executed at every pipeline run.
#' @param OUT_PATH name of the corresponding output raster (with path)
#' @param RES gridded output resolution (in degree)
#' @return creates a raster file in the /data folder (or else if specified) 
#' corresponding to the latest available Bio-ORACLE data;

build_bio_oracle <- function(OUT_PATH = "/net/meso/work/aschickele/Bluecloud_WB_local/data/bio_oracle",
                             RES = 1){
  
  # --- 1. Initialize minimum libraries
  # To avoid loading the entire pipeline libraries
  library(sdmpredictors)
  library(tidyverse)
  library(raster)
  library(virtualspecies)
  
  # --- 2. List all available layer names in Bio-ORACLE
  # Further selection can be done in the pipeline with ENV_VAR
  all_names <- list_layers() %>% 
    dplyr::filter(dataset_code == "Bio-ORACLE") %>% 
    dplyr::filter(version == "22") %>% 
    dplyr::filter(is_surface == "TRUE") %>% 
    dplyr::select(layer_code)
  
  # --- 3. Download all layers
  # They are stored in a cache somewhere. Not useful for us as we build a raster
  message(paste(Sys.time(), "--- Bio-ORACLE : Downloading the data"))
  env_r <- load_layers(all_names) %>% readAll()
  
  # --- 4. Re-grid to the desired resolution
  message(paste(Sys.time(), "--- Bio-ORACLE : Re-gridding the raster"))
  env_r <- aggregate(env_r, fact = c(RES/res(env_r), 1), fun = mean, na.rm = TRUE)
  env_r <- synchroniseNA(env_r)
  
  # --- 5. Save in the desired output
  writeRaster(env_r, OUT_PATH, overwrite = FALSE)
  
} # END FUNCTION




