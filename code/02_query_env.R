#' =============================================================================
#' @name query_env
#' @description appends the query_bio output with a list of environmental
#' values at the sampling stations and a path to the environmental raster that
#' will be used for projections
#' @param QUERY_BIO output object of 01b_query_bio.R
#' @param ENV_VAR vector of names of environmental variables available within
#' the climatologies available in Blue Cloud. If null all variables are taken.
#' @param ENV_PATH string, path to the .nc or raster of environmental variables
#' @return X: a dataframe of environmetal values at the sampling stations and 
#' @return ENV_PATH: the path to the environmental .nc or raster


query_env <- function(QUERY_BIO = query,
                      ENV_VAR = NULL,
                      ENV_PATH = "/net/meso/work/aschickele/Diversity/data/features_monthly"){
  
  # =============================== ENV QUERY ==================================
  # --- 1. Extract sample info
  # /!\ Depth is not taken into account for now, neither year (1990-2016 = WOA)
  sample <- QUERY_BIO$S %>% 
    dplyr::select(decimallongitude, decimallatitude, month)
  
  # --- 2. Open features gridded data and names
  # /!\ To change later for a direct query in a .nc file : avoid 2 different raster files
  features <- stack(ENV_PATH) %>% readAll()
  features_name <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    names()
  
  # --- 3. Extract the environmental data in the data frame
  X <- NULL
  for(j in 1:nrow(sample)){
    id <- grep(pattern = sample$month[j], names(features))
    xy <- sample[j,] %>% dplyr::select(x = decimallongitude, y = decimallatitude)
    
    tmp <- mclapply(features[[id]]@layers, function(a_layer) sample_raster_NA(a_layer, xy), mc.cores = 10) %>% 
      lapply(mean) %>% # average between two points if a coordinate is EXACTLY at an integer value (i.e., between two cells)
      as.data.frame()
    colnames(tmp) <- features_name
    X <- rbind(X, tmp)
  }
  
  # --- 4. Write X on the disk
  write_feather(X, paste0(project_wd, "/data/X.feather"))
  
  # --- 5. Append query_bio with the environmental values
  return(list(Y = QUERY_BIO$Y,
              X = X,
              S = QUERY_BIO$S,
              CALL = list(DATA_TYPE = QUERY_BIO$CALL$DATA_TYPE,
                          SP_SELECT = QUERY_BIO$CALL$SP_SELECT,
                          SAMPLE_SELECT = QUERY_BIO$CALL$SAMPLE_SELECT,
                          ENV_VAR = ENV_VAR,
                          ENV_PATH = ENV_PATH)))
  
  
} # END FUNCTION
