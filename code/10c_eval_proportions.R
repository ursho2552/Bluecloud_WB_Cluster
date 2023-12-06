#' =============================================================================
#' @name eval_proportions
#' @description sub-pipeline for model evaluation corresponding to presence data
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @return the MODEL object updated with evaluation metric values (model performance
#' metric and variable importance metric)
#' @return variable importance plots as PDF file

eval_proportions <- function(CALL,
                             QUERY,
                             MODEL){
  
  # --- 1. Model performance assessment
  # --- 1.1. Source the MBTR functions
  library(reticulate)
  source_python(paste0(project_wd,"/function/mbtr_function.py"))
  
  # --- 1.2. Load final model data and object
  final_fit <- py_load_object(MODEL$MBTR$final_fit, pickle = "pickle")[[1]]
  
  # --- 1.3. Extract observations and predictions
  x <- QUERY$FOLDS$test %>% 
    dplyr::select(QUERY$SUBFOLDER_INFO$ENV_VAR)
  
  y <- QUERY$FOLDS$test %>% 
    dplyr::select(as.character(CALL$SP_SELECT))
  
  y_hat <- mbtr_predict(final_fit, x, MODEL$MBTR$final_wf$NBOOST)
  
  # --- 1.4. Compute R2 into MODEL object
  MODEL[["MBTR"]][["eval"]][["R2"]] <- calc_rsquared(as.matrix(y), y_hat) %>% round(3)
  
  # --- 1.5. Compute an over-fitting rate
  # Calculated as the deviation during training to deviation in testing ratio
  resample_dev <- MODEL[["MBTR"]][["best_fit"]]
  eval_dev <- sqrt(mean(as.matrix((y-y_hat)^2), na.rm=TRUE))
  MODEL[["MBTR"]][["eval"]][["OVF"]] <- ((resample_dev/eval_dev)-1)*100
  
  # --- 2. Variable importance
  # Performed according to the loss evolution corresponding to each tree split
  # Re-fitting a model for every variable would be too costly for MBTR
  # --- 2.1. Initialize function
  var_imp <- matrix(0, 1, length(QUERY$SUBFOLDER_INFO$ENV_VAR)) %>% as.data.frame()
  
  n_tree <- length(final_fit$trees)
  
  get_split_loss <- function(t, nvar){
    tmp <- matrix(0, ncol = nvar, nrow = 1)
    n_node <-length(final_fit$trees[[t]]$g$nodes$`_nodes`)
    for(n in 1:n_node){
      var_nb <- final_fit$trees[[t]]$g$nodes$`_nodes`[[n]]$variable+1 #py index start at 0, R at 1
      if(!is.null(var_nb)){
        dloss <- final_fit$trees[[t]]$g$nodes$`_nodes`[[n]]$loss*(-1)
        tmp[var_nb] <- tmp[var_nb]+dloss
      }
      return(tmp)}
  } # end function
  
  # --- 2.2. Iteratively extract the loss
  # Per tree split x environmental variable
  var_imp0 <- mcmapply(FUN = get_split_loss,
                       t = 1:n_tree,
                       nvar = length(var_imp),
                       SIMPLIFY = FALSE,
                       USE.NAMES = FALSE,
                       mc.cores = 10)
  var_imp <- Reduce(`+`, var_imp0)
  colnames(var_imp) <- QUERY$SUBFOLDER_INFO$ENV_VAR
  
  # --- 2.3. Re-ordering and re_scaling
  var_imp_order <- order(var_imp, decreasing = TRUE)
  var_imp <- var_imp[,var_imp_order]
  var_imp <- var_imp/sum(var_imp)*100
  
  # --- 2.4. Compute the cumulative variable importance
  MODEL[["MBTR"]][["eval"]][["CUM_VIP"]] <- sum(var_imp[1:3])
  
  # --- 3. Removing low quality algorithm
  # In this case the only algorithm available...
  for(i in MODEL$MODEL_LIST){
    # --- 3.1. Based on model performance
    # Fixed at 0.1 for R2
    # /!\ /!\ /!\ /!\ treshold at 0.1 to be able to prototype something. TO CHANGE LATER OBVIOUSLY /!\ /!\ /!\ /!\
    if(MODEL[[i]][["eval"]][["R2"]] < 0.25 | is.na(MODEL[[i]][["eval"]][["R2"]])){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to R2 =", MODEL[[i]][["eval"]][["R2"]], "< 0.25 \n"))
    }
    
    # --- 3.2. Based on cumulative variable importance
    # Fixed at 30% cumulative importance for the top three predictors
    if(MODEL[[i]][["eval"]][["CUM_VIP"]] < 50 | is.na(MODEL[[i]][["eval"]][["CUM_VIP"]])){
      MODEL$MODEL_LIST <- MODEL$MODEL_LIST[MODEL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to CUM_VIP =", MODEL[[i]][["eval"]][["CUM_VIP"]], "< 50% \n"))
    }
  } # for each model loop
  
  # --- 4. Variable importance - Plot
  if(CALL$FAST == FALSE | (length(MODEL$MODEL_LIST) == 1)){
    # Define the color (green = QC passed; red = no)
    if(length(MODEL$MODEL_LIST) == 1){pal <- "#1F867B"
    } else {pal <- "#B64A60"}
    
    # Simple barplot of the MBTR vip if the model passed QC
    par(mar = c(20,20,3,2))
    barplot(var_imp, rep(1,length(var_imp)), axes = FALSE, horiz = TRUE, cex.names = 0.5,
            main = "PREDICTOR IMPORTANCE ( MBTR )",
            sub = paste("Predictive performance (R2) =", round(MODEL[["MBTR"]]$eval$R2, 2), "; Cumulated var. importance (%; top 3) =", round(MODEL[["MBTR"]]$eval$CUM_VIP, 0)), 
            col = pal, ylab = "", xlab = "Variable importance (%)", las = 2)
    axis(side = 1, at = seq(0, 100, 10), labels = seq(0, 100, 10), las = 2)
    abline(v = seq(0, 100, 10), lty = "dotted")
    box()
    box("figure", col = "black", lwd = 1)
  }

  # --- 5. Wrap up and save
  return(MODEL)
  
} # END FUNCTION
