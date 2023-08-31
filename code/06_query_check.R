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
  set.seed(123)
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : query_check ********************"))
  # --- 1.2. Load the run metadata and query
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 1.3. Moving average function
  # Short and only used here, thats why it is not in the function folder
  ma <- function(x, n = 10){stats::filter(x, rep(1 / n, n), sides = 2)}
  
  # --- 2. Outlier analysis
  # Outlier check on the query based on z-score (from Nielja code)
  if(OUTLIER == TRUE){
    if(CALL$DATA_TYPE == "binary"){
      message("--- Cannot perform outlier analysis on presence - pseudo absence data ---")
    } else {
      to_remove <- outlier_iqr_col(QUERY$Y, n = 2.5) %>% as.vector()
      if(length(to_remove > 0)){
        # QUERY$Y <- QUERY$Y[-to_remove,]
        QUERY$Y <- dplyr::slice(QUERY$Y, -to_remove)
        QUERY$X <- dplyr::slice(QUERY$X, -to_remove)
        QUERY$S <- dplyr::slice(QUERY$S, -to_remove)
      } # if to remove !NULL
      
      message(paste("--- OUTLIERS : Removed row number", to_remove, "\n"))
    }
  } # END if outlier TRUE
  
  # --- 2. Environmental variable correlation check
  # Removing correlated environmental variables to avoid correlated model features
  if(is.numeric(CALL$ENV_COR) == TRUE){
    # --- 2.1. Opening environmental value at presence points
    features <- QUERY$X
    if(CALL$DATA_TYPE == "binary"){
      features <- features[which(QUERY$Y$measurementvalue != 0),]
    }
    
    # --- 2.2. Re-order features by Variation Inflation Factor
    # Later, we will keep the variable with the lower VIF among correlated clusters
    features_vif <- vif(features) %>% 
      arrange(VIF) %>% 
      dplyr::filter(!is.infinite(VIF) & !is.na(VIF))
    features <- features[, features_vif$Variables]
    
    # --- 2.3. Check correlation/distance between variables
    features_dist <- cor(features, method = "pearson")
    features_dist <- as.dist(1-abs(features_dist))
    
    # --- 2.4. Do a clustering and cut
    features_clust <- hclust(features_dist) %>% as.dendrogram()
    features_group <- cutree(features_clust, h = 1-CALL$ENV_COR)
    
    # --- 2.5. Randomly choose one variable within each inter-correlated clusters
    features_keep <- NULL
    for(i in 1:max(features_group)){
      tmp <- which(features_group == i)
      message(paste("--- ENV_COR : Cluster", i, ": Keeping", names(tmp[1]), "\n"))
      features_keep <- c(features_keep, features_group[tmp[1]])
      if(length(tmp) > 1){
        out <- tmp[-1]
        message(paste("--- ENV_COR : Cluster", i, ": Removing", names(out), "\n"))
      }
    }
    
    # --- 2.6. Plot the corresponding dentrogram
    pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/02_env_cor.pdf"))
    par(mfrow = c(2,1), mar = c(8,5,5,3), cex = 0.6)
    pal <- rep("#B64A60", length(features))
    pal[get_leaves_attr(features_clust, "label") %in% names(features_keep)] <- "#1F867B"
    labels_colors(features_clust) = pal
    plot(features_clust, axes = FALSE, main = "Env. variable Pearson's correlation (r) at the sampling stations", cex = 0.5)
    axis(side = 2, at = seq(0,1,0.2), labels = seq(1,0,-0.2), las = 1, cex.axis = 0.6)
    abline(h = seq(0,1,0.2), col = "gray50", lty = "dashed")
    abline(h = 1-CALL$ENV_COR, col = "#B64A60")
    dev.off()
    
    # --- 2.7. Update ENV_VAR
    QUERY$SUBFOLDER_INFO$ENV_VAR <- names(features_keep)
  } # END if env_cor TRUE
  
  
  # --- 3. Univariate variable importance analysis
  # Done with a Random forest using the method developed in the "Caret" library
  if(UNIVARIATE == TRUE){
    if(CALL$DATA_TYPE == "proportions"){
      message("A univariate predictor selection is not possible for proportion data,
              please select carefully your predictors")
    } else {
      # --- 3.1. Initialize data and control parameters
      features <- QUERY$X[,QUERY$SUBFOLDER_INFO$ENV_VAR]
      rfe_df <- cbind(QUERY$Y, features)
      # Modify rfFuncs : we modify the selectSize code to use pickSizeBest instead of a 1.5 tolerance around the best
      rfFuncs$selectSize <- function(...) pickSizeBest(...) # weird functioning :-)
      rfe_control <- rfeControl(functions=rfFuncs, method="cv", number=5, rerank = FALSE)
      # --- 3.2. Run the RFE algorithm
      message(paste(Sys.time(), "--- UNIVARIATE : Fitting the RFE algorithm \n"))
      rfe_fit <- rfe(rfe_df[,-1], rfe_df[,1], sizes = c(1:ncol(rfe_df[,-1])), rfeControl = rfe_control)
      
      # --- 3.3. Extract the relevant predictors
      # --- 3.3.1. Compute the moving average loss
      loss_ma <- ma(rfe_fit$results$RMSE, n = 10)
      # --- 3.3.2. Compute the percentage loss
      loss_ma_pct <- (loss_ma[-1] - loss_ma[-length(loss_ma)])/loss_ma[-length(loss_ma)]*100
      # --- 3.3.3. Find the first minimum or <1% loss percentage
      # Consider all variables if "id" is NA, due to moving average not working for low variable number
      id <- which(loss_ma_pct > -1)[1]
      if(is.na(id) == TRUE){
        message("--- The moving window could not find a minimum loss, please considering 
                removing this option from your run due to insufficient number of environmental variables")
        id <- ncol(features)
      }
      
      # --- 3.4. Compute variable importance
      rfe_vip <- rfe_fit$variables %>% 
        dplyr::select(var, Overall) %>% 
        mutate(var = fct_reorder(as.factor(var), Overall, .desc = TRUE))
      
      # --- 3.5. Extract the environmental variables
      ENV_VAR <- rfe_vip$var[1:id] %>% as.character()
      message(paste("--- UNIVARIATE : Selecting", ENV_VAR, "\n"))
      
      # --- 3.6. Produce an information plot
      pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/03_univariate_predictor_selection.pdf"))
      par(mfrow = c(2,1))
      # --- Variable importance
      boxplot(rfe_vip$Overall ~ rfe_vip$var, main = paste("Univariate predictor importance for:", SUBFOLDER_NAME), 
              xlab = "", ylab = "", axes = FALSE, outline = FALSE,
              col = c(rep("#1F867B", id), rep("#B64A60", ncol(features)-id)))
      axis(side = 1, at = 1:ncol(features), labels = levels(rfe_vip$var), las = 2, cex.axis = 0.3)
      axis(side = 2, at = c(seq(0, 15, 5), seq(0, 100, 20)), labels = c(seq(0, 15, 5), seq(0, 100, 20)), las = 2)
      abline(h = c(seq(0, 15, 5), seq(0, 100, 20)), lty = "longdash", col = "gray50")
      # --- Number of variables
      plot(rfe_fit$results$Variables, rfe_fit$results$RMSE, pch = 20, lwd = 3,
           col = c(rep("#1F867B", id), rep("#B64A60", ncol(features)-id)),
           main = paste("Optimal predictor number for:", SUBFOLDER_NAME), xlab = "Nb. of predictors", ylab = "RMSE")
      grid(col = "gray50")
      dev.off()
      
      # --- 3.7. Update ENV_VAR
      QUERY$SUBFOLDER_INFO$ENV_VAR <- ENV_VAR
      
    }
  } # END if univariate TRUE

  # --- 4. MESS analysis
  # --- 4.1. Load necessary data
  features <- stack(CALL$ENV_PATH) %>% 
    readAll() %>% 
    raster::subset(QUERY$SUBFOLDER_INFO$ENV_VAR)
  names(features) <- QUERY$SUBFOLDER_INFO$ENV_VAR
  
  # --- 4.2. Compute the mess analysis
  # --- 4.2.1. Load environmental samples data
  tmp <- QUERY$X %>% dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
  
  # --- 4.2.2. Remove the pseudo-absences from the analysis
  # We should not consider them as true input data as they are user defined
  if(CALL$DATA_TYPE == "binary"){
    tmp <- tmp[which(QUERY$Y != 0),]
  } # if pres remove pseudo abs
  
  # --- 4.2.3. Analysis
  r_mess <- dismo::mess(x = features, v = as.data.frame(tmp), full = FALSE)
  
  # --- 4.3. Append to query
  QUERY$MESS <- r_mess
  
  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
