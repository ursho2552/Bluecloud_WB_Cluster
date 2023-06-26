#' =============================================================================
#' @name query_check
#' @description This function performs a series of last check before modelling,
#' including an outlier check, environmental variable correlation and MESS
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param OUTLIER if TRUE, remove outliers
#' @param UNIVARIATE if true, performs a univariate predictor pre-selection
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

query_check <- function(FOLDER_NAME = NULL,
                        SUBFOLDER_NAME = NULL,
                        OUTLIER = TRUE,
                        UNIVARIATE = TRUE){
  
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
  
  # --- 2. Univariate variable importance analysis
  # Done with a Random forest using the method developed in the "Caret" library
  if(UNIVARIATE == TRUE){
    if(CALL$DATA_TYPE == "omic"){
      message("A univariate predictor selection is not possible for omic-proportion data,
              please select carefully your predictors")
    } else{
      # --- 2.1. Initialize data and control parameters
      rfe_df <- cbind(QUERY$Y, QUERY$X)
      rfe_control <- rfeControl(functions=rfFuncs, method="cv", number=10)
      # --- 2.2. Run the RFE algorithm
      message(paste(Sys.time(), "--- UNIVARIATE : Fitting the RFE algorithm \n"))
      rfe_fit <- rfe(rfe_df[,-1], rfe_df[,1], sizes = c(1:ncol(rfe_df[,-1])), rfeControl = rfe_control)
      # --- 2.3. Extract the relevant predictors
      ENV_VAR <- caret::predictors(rfe_fit)
      message(paste("--- UNIVARIATE : Selecting", ENV_VAR, "\n"))
      
      # --- 2.4. Produce an information plot
      rfe_vip <- rfe_fit$variables %>% 
        dplyr::select(var, Overall) %>% 
        mutate(var = fct_reorder(as.factor(var), Overall, .desc = TRUE))
      
      pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/univariate_predictor_selection.pdf"))
      par(mfrow = c(2,1))
      # --- Variable importance
      boxplot(rfe_vip$Overall ~ rfe_vip$var, main = "Univariate predictor importance", 
              ylim = c(0,100), xlab = "", ylab = "", axes = FALSE)
      axis(side = 1, at = 1:ncol(QUERY$X), labels = levels(rfe_vip$var), las = 2)
      axis(side = 2, at = seq(0, 100, 20), labels = seq(0, 100, 20), las = 2)
      abline(h = c(0, 20, 40, 60, 80, 100), lty = "longdash", col = "gray50")
      # --- Number of variables
      plot(rfe_fit$results$Variables, rfe_fit$results$RMSE, pch = 20, lwd = 3,
           col = c(rep("green", rfe_fit$optsize), rep("red", ncol(QUERY$X)-rfe_fit$optsize)),
           main = "Optimal predictor number", xlab = "Nb. of predictors", ylab = "RMSE")
      grid(col = "gray50")
      dev.off()
    }
  } # END if univariate TRUE

  # --- 3. Environmental variable correlation check
  # Removing correlated environmental variables to avoid correlated model features
  if(is.numeric(CALL$ENV_COR) == TRUE){
    # --- 3.1. Opening environmental value at presence points
    if(UNIVARIATE == TRUE){features <- QUERY$X[,ENV_VAR]
    } else {features <- QUERY$X}
    
    # --- 3.2. Check correlation/distance between variables
    features_dist <- cor(features, method = "pearson")
    features_dist <- as.dist(1-abs(features_dist))
    
    # --- 3.3. Do a clustering and cut
    features_clust <- hclust(features_dist) %>% as.dendrogram()
    features_group <- cutree(features_clust, h = 1-CALL$ENV_COR)
    
    # --- 3.4. Randomly choose one variable within each inter-correlated clusters
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
    
    # --- 3.5. Plot the corresponding dentrogram
    pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/env_cor.pdf"))
    pal <- rep("red", length(features))
    pal[get_leaves_attr(features_clust, "label") %in% names(features_keep)] <- "green"
    labels_colors(features_clust) = pal
    plot(features_clust, axes = FALSE, main = "Env. variable Pearson's correlation (r) at the sampling stations")
    axis(side = 2, at = seq(0,1,0.2), labels = seq(1,0,-0.2), las = 1)
    abline(h = seq(0,1,0.2), col = "gray50", lty = "dashed")
    abline(h = 1-CALL$ENV_COR, col = "red")
    dev.off()
    
    # --- 3.6. Update ENV_VAR
    QUERY$SUBFOLDER_INFO$ENV_VAR <- names(features_keep)
  } # END if env_cor TRUE

  # --- 4. MESS analysis
  # --- 4.1. Load necessary data
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
    readAll() %>% 
    raster::subset(QUERY$SUBFOLDER_INFO$ENV_VAR)
  
  # --- 4.2. Compute the mess analysis
  tmp <- QUERY$X %>% dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
  r_mess <- dismo::mess(x = features, v = tmp, full = FALSE)
  
  # --- 4.3. Append to query
  QUERY$MESS <- r_mess
  
  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
