#' =============================================================================
#' @name query_check
#' @description This function performs a series of last check before modelling,
#' including an outlier check, environmental variable correlation and MESS
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

query_check <- function(FOLDER_NAME = NULL,
                        SUBFOLDER_NAME = NULL){
  
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
  ma <- function(x, n = 10){
                 if(length(x) > n){ma = stats::filter(x, rep(1 / n, n), sides = 2)
                 } else {ma = NA}
                 return(ma)
    } # end function
  
  # --- 2. Outlier analysis
  # Outlier check on the query based on z-score (from Nielja code)
  if(CALL$OUTLIER == TRUE){
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
    # Remove variables with less than X unique values (to be able to fit the lm within VIF)
    tmp <- which(apply(features, 2, function(x)(length(unique(x)))) <= 10) %>% as.numeric()
    # Later, we will keep the variable with the lower VIF among correlated clusters
    if(length(tmp) != 0){
      features_vif <- usdm::vif(features[,-tmp]) %>%
        arrange(VIF)
      features <- features[, features_vif$Variables]
    }
    
    # --- 2.3. Check correlation/distance between variables
    features_dist <- cor(features, method = "pearson")
    features_dist <- as.dist(1-abs(features_dist))
    
    # --- 2.4. Do a clustering and cut
    features_clust <- hclust(features_dist) %>% as.dendrogram()
    features_group <- cutree(features_clust, h = 1-CALL$ENV_COR)
    
    # --- 2.5. Choose one variable within each inter-correlated clusters
    # Lowest VIF if possible to compute it
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
    
    # --- 2.6. Plot the corresponding dendrogram
    pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/02_env_cor.pdf"))
    par(mfrow = c(1,1), mar = c(3,2,4,30), cex = 0.6)
    pal <- rep("#B64A60", length(features))
    pal[get_leaves_attr(features_clust, "label") %in% names(features_keep)] <- "#1F867B"
    labels_colors(features_clust) = pal
    plot(features_clust, horiz = TRUE, axes = FALSE,
         main = "ENVIRONMENTAL PREDICTORS \n Pearson's correlation (r) at the sampling stations", cex = 0.5)
    axis(side = 1, at = seq(0,1,0.2), labels = seq(1,0,-0.2), las = 1, cex.axis = 1)
    abline(v = seq(0,1,0.2), col = "gray50", lty = "dashed")
    abline(v = 1-CALL$ENV_COR, col = "#B64A60")
    dev.off()
    
    # --- 2.7. Update ENV_VAR
    QUERY$SUBFOLDER_INFO$ENV_VAR <- names(features_keep)
  } # END if env_cor TRUE
  
  # # --- 3. Univariate variable importance analysis
  # Done with a Random forest using the method developed in the "Caret" library
  if(CALL$UNIVARIATE == TRUE){
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
      pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/03_feature_pre_selection.pdf"))
      par(mfrow = c(2,1), mar = c(4,4,3,20))
      # --- Variable importance
      boxplot(rfe_vip$Overall ~ rfe_vip$var,
              main = paste("ENVIRONMENTAL PREDICTORS \n A-priori importance for ID:", SUBFOLDER_NAME),
              xlab = "Estimated importance (%)", ylab = "", axes = FALSE, outline = FALSE, horizontal = TRUE,
              col = c(rep("#1F867B", id), rep("#B64A60", ncol(features)-id)), cex.main = 0.7, cex.lab = 0.7)
      axis(side = 4, at = 1:ncol(features), labels = levels(rfe_vip$var), las = 2, cex.axis = 0.8)
      axis(side = 1, at = c(seq(0, 15, 5), seq(0, 100, 20)), labels = c(seq(0, 15, 5), seq(0, 100, 20)))
      abline(v = c(seq(0, 15, 5), seq(0, 100, 20)), lty = "longdash", col = "gray50")
      box()
      box("figure", col="black", lwd = 1)
      # --- Number of variables
      par(mar = c(4,4,3,20))
      plot(rfe_fit$results$RMSE, rfe_fit$results$Variables, pch = 18, cex = 2,
           col = c(rep("#1F867B", id), rep("#B64A60", ncol(features)-id)),
           main = paste("ENVIRONMENTAL PREDICTORS \n Optimal number for ID:", SUBFOLDER_NAME),
           ylab = "Nb. of considered predictors", xlab = "Loss metric (RMSE)", cex.main = 0.7, cex.lab = 0.7)
      axis(side = 4, at = 1:ncol(features), labels = levels(rfe_vip$var), las = 2, cex.axis = 0.8)
      grid(col = "gray50")
      dev.off()

      # --- 3.7. Update ENV_VAR
      QUERY$SUBFOLDER_INFO$ENV_VAR <- ENV_VAR

    }
  } # END if univariate TRUE

  # --- 4. MESS analysis
  r_mess <- NULL
  for(m in 1:length(CALL$ENV_DATA)){
    # --- 4.1. Load necessary data
    features <- CALL$ENV_DATA[[m]] %>%
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
    r_mess[[m]] <- dismo::mess(x = features, v = as.data.frame(tmp), full = FALSE)

  } # for m month

  # --- 4.3. Append to query
  QUERY$MESS <- stack(r_mess)
  
  # --- 5. Verification of feature pre-selection
  if(CALL$DATA_TYPE != "proportions"){
    # --- 5.1. Initialize plot
    pal <- rep("#B64A60", ncol(QUERY$X))
    pal[which(names(QUERY$X) %in% names(features_keep) == TRUE)] <- "orange" # passed correlation
    if(CALL$UNIVARIATE == TRUE){
      pal[which(names(QUERY$X) %in% QUERY$SUBFOLDER_INFO$ENV_VAR == TRUE)] <- "#1F867B" # passed RFE
    } else {
      pal[which(names(QUERY$X) %in% names(features_keep) == TRUE)] <- "#1F867B" # passed correlation
    }
    
    # --- 5.2. Iterate data preparation and save in list
    plot_data <- lapply(1:ncol(QUERY$X), function(i){
      
      x0 <- as.numeric(QUERY$X[,i]) %>% unlist()
      q <- quantile(x0, probs = seq(0, 1, length.out = 51)) %>% unique()
      q_id <- cut(x0, q, include.lowest = TRUE, labels = FALSE)
      x <- q[q_id]
      
      df <- data.frame(x = x, y = QUERY$Y$measurementvalue) %>%
        dplyr::group_by(x) %>%
        dplyr::summarise(y = mean(y, na.rm = TRUE))
      p <- cor(df$x, df$y, method = "spearman", use = "pairwise.complete.obs") %>% round(2)
      return(list(df = df, p = p))
    }) # end function
    
    # --- 5.3. Compute the QC
    # We calculate the average absolute spearman correlation between obs. and env. (i.e., meaningful feature set)
    tmp <- plot_data[which(names(QUERY$X) %in% QUERY$SUBFOLDER_INFO$ENV_VAR == TRUE)]
    PRE_VIP <- lapply(1:length(tmp), function(x){
      out <- tmp[[x]]$p %>% abs()
    }) %>% unlist() %>% mean() %>% round(2)
    
    # --- 5.4. Do the plot
    pdf(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/04_feature_verification.pdf"))
    par(mfrow = c(4,4), mar = c(4,1,3,1))
    # --- 5.4.1. Plot QC and legend
    plot.new()
    text(x = 0.5, y = 0.5, paste("FEATURE PRE-SELECTION \n", QUERY$annotations$scientificname, "\n ID:", QUERY$annotations$worms_id), cex = 1)
    par(xpd = TRUE, mar = c(4,4,3,4))
    plot(x = c(0.1,0.1,0.1), y = c(0.8, 0.5, 0.2), cex = 3, pch = 15, col = c("#B64A60", "orange", "#1F867B"), xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "", main = "Quantile average")
    text(x = 0.1, y = 0.8, "Excluded due to intercollinearity", pos = 4, offset = 1, cex = 0.8)
    text(x = 0.1, y = 0.5, "Do not bring additional information", pos = 4, offset = 1, cex = 0.8)
    text(x = 0.1, y = 0.2, "Selected feature" , pos = 4, offset = 1, cex = 0.8)
    text(x = 0.1, y = 0, paste("Avg. Spearman:", PRE_VIP) , pos = 4, offset = 1, cex = 0.8)
    plot(x = 0.1, y = 0.8, cex = 2, pch = 20, col = "black", xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "", main = "Observations")
    text(x = 0.1, y = 0.8, "Raw values", pos = 4, offset = 1, cex = 0.8)
    plot.new()
    
    # --- 5.4.2. Plot all single obs. vs env. with quantiles average and raw data
    par(mar = c(5,4,3,1))
    for(i in 1:ncol(QUERY$X)){
      plot(QUERY$X[,i], QUERY$Y$measurementvalue, xlab = "feature value", ylab = "target observations",
           main = names(QUERY$X)[i], cex.main = 0.8,
           sub = paste("Spearman:", plot_data[[i]]$p),
           pch = 20, col = alpha("black", 0.1), cex = 1)
      points(plot_data[[i]]$df$x, plot_data[[i]]$df$y, pch = 20, col = alpha(pal[i], 0.5), cex = 2)
      box("figure", col="black", lwd = 1)
    } # end for i
    
    # --- 5.4.3. Stop pdf and add QC to the query
    dev.off()
    QUERY[["eval"]][["PRE_VIP"]] <- PRE_VIP
  } # end if !proportions
  
  # --- 6. Wrap up and save
  # --- 6.1. Save file(s)
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  # --- 6.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 6.3. Pretty return
  # Updates the list of species to model depending on the PRE_VIP (> 0.25) to avoid non meaningful feature ensembles
  if(QUERY$eval$PRE_VIP >= 0.25 | CALL$FAST == FALSE){
    return(SUBFOLDER_NAME)
  } else {
    message("The selected features do not present significant trends to the observations; please work on the predictors and data. \n
            The species is discarded from further analysis")
    return(NA)
  }
  
} # END FUNCTION
