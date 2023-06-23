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
  
  # --- 2. Environmental variable correlation check
  # Removing correlated environmental variables to avoid correlated model features
  if(is.numeric(CALL$ENV_COR) == TRUE){
    # --- 2.1. Opening environmental value at presence points
    features <- QUERY$X
    
    # --- 2.2. Check correlation/distance between variables
    features_dist <- cor(features, method = "pearson")
    features_dist <- as.dist(1-abs(features_dist))
    
    # --- 2.3. Do a clustering and cut
    features_clust <- hclust(features_dist) %>% as.dendrogram()
    features_group <- cutree(features_clust, h = 1-CALL$ENV_COR)
    
    # --- 2.6. Randomly choose one variable within each inter-correlated clusters
    features_keep <- features_group
    for(i in 1:max(features_group)){
      tmp <- which(features_group == i)
      if(length(tmp) > 1){
        message(paste("--- ENV_COR : Cluster", i, ": Keeping", names(tmp[1]), "\n"))
        tmp <- tmp[-1]
        message(paste("--- ENV_COR : Cluster", i, ": Removing", names(tmp), "\n"))
        features_keep <- features_keep[-tmp]
      }
    }
    
    # --- 2.5. Plot the corresponding dentrogram
    pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/env_cor.pdf"))
    pal <- rep("red", length(features))
    pal[get_leaves_attr(features_clust, "label") %in% names(features_keep)] <- "green"
    labels_colors(features_clust) = pal
    plot(features_clust, axes = FALSE, main = "Env. variable Pearson's correlation (r) at the sampling stations")
    axis(side = 2, at = seq(0,1,0.2), labels = seq(1,0,-0.2), las = 1)
    abline(h = seq(0,1,0.2), col = "gray50", lty = "dashed")
    abline(h = 1-CALL$ENV_COR, col = "red")
    dev.off()
    
    # --- 2.6. Update ENV_VAR
    QUERY$SUBFOLDER_INFO$ENV_VAR <- names(features_keep)
  } # END if env_cor TRUE
  
  # --- 3. MESS analysis
  # --- 3.1. Load necessary data
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    readAll() %>% 
    raster::subset(QUERY$SUBFOLDER_INFO$ENV_VAR)
  
  # --- 3.2. Compute the mess analysis
  tmp <- QUERY$X %>% dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
  r_mess <- dismo::mess(x = features, v = tmp, full = FALSE)
  
  # --- 3.3. Append to query
  QUERY$MESS <- r_mess
  
  # --- 4. Wrap up and save
  # --- 4.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 4.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
