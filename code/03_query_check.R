#' =============================================================================
#' @name query_check
#' @description This function performs a series of last check before modelling,
#' including an outlier check, environmental variable correlation and MESS
#' @param QUERY the query output from the previous functions
#' @param OUTLIER if TRUE, remove outliers
#' @param ENV_COR numeric, removes the correlated environmental values from the
#' query objects and CALL according to the defined threshold. Else NULL.
#' @param MESS if TRUE, performs a MESS analysis that is stored in the query
#' @return the updated query list

query_check <- function(QUERY = query,
                        OUTLIER = TRUE,
                        ENV_COR = 0.8, 
                        MESS = TRUE){
  
  # =============================== OUTLIER ANALYSIS ===========================
  # --- 1. Outlier check on the query based on z-score (from Nielja code)
  if(OUTLIER == TRUE){
    if(QUERY$CALL$DATA_TYPE == "pres"){
      message("--- Cannot perform outlier analysis on presence - pseudo absence data ---")
    } else {
      to_remove <- outlier_iqr_col(QUERY$Y, n = 2.5)
      QUERY$Y <- QUERY$Y[-to_remove,]
      QUERY$X <- QUERY$X[-to_remove,]
      QUERY$S <- QUERY$S[-to_remove,]
      message(paste("--- OUTLIERS : Removed row number", to_remove, "\n"))
    }
  } # END if outlier TRUE
  
  # ======================== ENV VARIABLE CORRELATION ==========================
  # Deletes correlated variables - SHOULD BE MOVED EARLIER PROBALY
  # TO FIX with the CALL to environmental variable raster

  if(is.numeric(ENV_COR) == TRUE){
    # --- 1. Open raster as dataframe
    features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
      readAll() %>% 
      rasterToPoints() %>% 
      as.data.frame() %>% 
      dplyr::select(-c(x, y)) %>% 
      tidyr::drop_na()
    
    # --- 2. Check correlation/distance between variables
    features_dist <- cor(features) %>% as.dist()
    
    # --- 3. Do a clustering and cut
    features_clust <- hclust(-features_dist) # dist = 1/cor
    plot(features_clust)
    abline(h = -ENV_COR, col = "red")
    
    features_group <- cutree(features_clust, h = -ENV_COR)
    
    # --- 4. Randomly choose one variable within each inter-correlated clusters
    features_keep <- features_group
    for(i in 1:max(features_group)){
      tmp <- which(features_group == i)
      if(length(tmp) > 1){
        message(paste("--- Cluster", i, ": Keeping", names(tmp[1])))
        tmp <- tmp[-1]
        message(paste("--- Cluster", i, ": Removing", names(tmp)))
        features_keep <- features_keep[-tmp]
      }
    }
    
    # --- 5. Update QUERY$CALL$ENV_VAR
    QUERY$CALL$ENV_VAR <- names(features_keep)
    
  } # END if env_cor TRUE
  
  # ============================== MESS ANALYSIS ===============================
  # TO FIX : probably have to move that to the training set to be more accurate
  
  if(MESS == TRUE){
    # --- 1. Load necessary data
    features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
      readAll() %>% 
      raster::subset(QUERY$CALL$ENV_VAR)
    
    # --- 2. Compute the mess analysis
    tmp <- QUERY$X %>% dplyr::select(QUERY$CALL$ENV_VAR)
    r_mess <- dismo::mess(x = features, v = tmp)
    
    # --- 3. Append to query
    QUERY$MESS <- r_mess
    
  } # END if mess TRUE
  
  return(query = QUERY)
  
} # END FUNCTION
