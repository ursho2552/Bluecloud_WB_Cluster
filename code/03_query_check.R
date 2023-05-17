#' =============================================================================
#' @name query_check
#' @description This function performs a series of last check before modelling,
#' including an outlier check, environmental variable correlation and MESS
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param OUTLIER if TRUE, remove outliers
#' @param ENV_COR numeric, removes the correlated environmental values from the
#' query objects and CALL according to the defined threshold. Else NULL.
#' @param MESS if TRUE, performs a MESS analysis that is stored in the query
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

query_check <- function(FOLDER_NAME = NULL,
                        SUBFOLDER_NAME = NULL,
                        OUTLIER = TRUE,
                        ENV_COR = 0.8, 
                        MESS = TRUE){
  
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
      QUERY$Y <- QUERY$Y[-to_remove,]
      QUERY$X <- QUERY$X[-to_remove,]
      QUERY$S <- QUERY$S[-to_remove,]
      message(paste("--- OUTLIERS : Removed row number", to_remove, "\n"))
    }
  } # END if outlier TRUE
  
  # --- 3. Environmental variable correlation
  # Removing correlated environmental variables to avoid correlated model features
  if(is.numeric(ENV_COR) == TRUE){
    # --- 3.1. Open raster as dataframe
    features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
      readAll() %>% 
      rasterToPoints() %>% 
      as.data.frame() %>% 
      dplyr::select(-c(x, y)) %>% 
      tidyr::drop_na()
    
    # --- 3.2. Check correlation/distance between variables
    features_dist <- cor(features) %>% as.dist()
    
    # --- 3.3. Do a clustering and cut
    features_clust <- hclust(-features_dist) # dist = 1/cor
    plot(features_clust)
    abline(h = -ENV_COR, col = "red")
    
    features_group <- cutree(features_clust, h = -ENV_COR)
    
    # --- 3.4. Randomly choose one variable within each inter-correlated clusters
    features_keep <- features_group
    for(i in 1:max(features_group)){
      tmp <- which(features_group == i)
      if(length(tmp) > 1){
        message(paste("--- ENV_COR : Cluster", i, ": Keeping", names(tmp[1])))
        tmp <- tmp[-1]
        message(paste("--- ENV_COR : Cluster", i, ": Removing", names(tmp)))
        features_keep <- features_keep[-tmp]
      }
    }
    
    # --- 3.5. Update CALL$ENV_VAR
    CALL$ENV_VAR <- names(features_keep)
  } # END if env_cor TRUE
  
  # --- 4. MESS analysis
  # TO FIX : probably have to move that to the training set to be more accurate ?
  if(MESS == TRUE){
    # --- 4.1. Load necessary data
    features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
      readAll() %>% 
      raster::subset(CALL$ENV_VAR)
    
    # --- 4.2. Compute the mess analysis
    tmp <- QUERY$X %>% dplyr::select(CALL$ENV_VAR)
    r_mess <- dismo::mess(x = features, v = tmp, full = FALSE)
    
    # --- 4.3. Append to query
    QUERY$MESS <- r_mess
  } # END if mess TRUE
  
  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  save(CALL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
