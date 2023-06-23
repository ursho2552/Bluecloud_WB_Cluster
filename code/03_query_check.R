#' =============================================================================
#' @name query_check
#' @description This function performs a series of last check before modelling,
#' including an outlier check, environmental variable correlation and MESS
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param OUTLIER if TRUE, remove outliers
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

query_check <- function(FOLDER_NAME = NULL,
                        SUBFOLDER_NAME = NULL,
                        OUTLIER = TRUE){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : query_check ********************"))
  # --- 1.2. Load the run metadata and query
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 2. Outlier analysis
  # Outlier check on the query based on z-score (from Nielja code)
  if(OUTLIER == TRUE){
    if(CALL$DATA_TYPE == "pres"){
      message("--- Cannot perform outlier analysis on presence - pseudo absence data ---")
    } else {
      to_remove <- outlier_iqr_col(QUERY$Y, n = 2.5)
      if(length(to_remove > 0)){
        QUERY$Y <- QUERY$Y %>% slice(-to_remove)
        QUERY$X <- QUERY$X %>% slice(-to_remove)
        QUERY$S <- QUERY$S %>% slice(-to_remove)
      } # if to remove !NULL
      
      message(paste("--- OUTLIERS : Removed row number", to_remove, "\n"))
    }
  } # END if outlier TRUE
  
  # --- 3. MESS analysis
  # --- 3.1. Load necessary data
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    readAll() %>% 
    raster::subset(CALL$ENV_VAR)
  
  # --- 3.2. Compute the mess analysis
  tmp <- QUERY$X %>% dplyr::select(all_of(CALL$ENV_VAR))
  r_mess <- dismo::mess(x = features, v = tmp, full = FALSE)
  
  # --- 3.3. Append to query
  QUERY$MESS <- r_mess
  
  # --- 4. Wrap up and save
  # --- 4.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 4.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
