#' =============================================================================
#' @name eval_proportions
#' @description sub-pipeline for model evaluation corresponding to presence data
#' @param CALL the call object from the master pipeline
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline
#' @param PDF_PATH preset path for the pdf save
#' @return the MODEL object updated with evaluation metric values (model performance
#' metric and variable importance metric)
#' @return variable importance plots as PDF file

eval_proportions <- function(CALL,
                             QUERY,
                             MODEL,
                             PDF_PATH){

  # --- 1. Model performance assessment
  # --- 1.1. Source the MBTR functions
  library(reticulate)
  source_python(paste0(project_wd,"/function/mbtr_function.py"))

  # --- 1.2. Loop over the cross validation folds
  y_all <- NULL
  y_hat_all <- NULL
  
  for(cv in 1:CALL$NFOLD){
    # --- 1.2.1. Load final model data and object
    final_fit <- py_load_object(paste0(MODEL$MBTR$final_fit[cv],"_m"), pickle = "pickle")[[1]]
    
    # --- 1.2.1. Extract observations and predictions
    x <- read_feather(paste0(project_wd, "/data/MBTR_cache/", cv, "_X_val.feather"))
    y <- read_feather(paste0(project_wd, "/data/MBTR_cache/", cv, "_Y_val.feather"))
    y_hat <- mbtr_predict(final_fit, x, MODEL$MBTR$final_wf$NBOOST)
    
    # --- 1.2.2. Concatenate
    y_all <- rbind(y_all, as.matrix(y))
    y_hat_all <- rbind(y_hat_all, y_hat)
    
  } # end CV loop
  
  # --- 1.3. Compute R2 into MODEL object
  MODEL[["MBTR"]][["eval"]][["R2"]] <- calc_rsquared(y_all, y_hat_all) %>% round(3)

  # --- 2. Variable importance
  # Performed according to the loss evolution corresponding to each tree split
  # Re-fitting a model for every variable would be too costly for MBTR
  
  # Loop over cross validation folds again
  var_imp <- list()
  for(cv in 1:CALL$NFOLD){
    
    # --- 2.1. Reload final fit
    final_fit <- py_load_object(paste0(MODEL$MBTR$final_fit[cv],"_m"), pickle = "pickle")[[1]]
    
    # --- 2.2. Initialize function
    # var_imp <- matrix(0, 1, length(QUERY$SUBFOLDER_INFO$ENV_VAR)) %>% as.data.frame()
    n_tree <- min(MODEL$MBTR$final_wf$NBOOST, final_fit$trees %>% length())
    
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
    
    # --- 2.3. Iteratively extract the loss
    # Per tree split x environmental variable
    var_imp0 <- mcmapply(FUN = get_split_loss,
                         t = 1:n_tree,
                         nvar = length(QUERY$SUBFOLDER_INFO$ENV_VAR),
                         SIMPLIFY = FALSE,
                         USE.NAMES = FALSE,
                         mc.cores = 10)
    var_imp[[cv]] <- Reduce(`+`, var_imp0)
    colnames(var_imp[[cv]]) <- QUERY$SUBFOLDER_INFO$ENV_VAR
    
    # --- 2.4. Re-scaling
    var_imp[[cv]] <- var_imp[[cv]]/sum(var_imp[[cv]])*100
    
  } # end CV loop
  
  # --- 2.3. Transform to a matrix and re-order
  var_imp <- Reduce(`rbind`, var_imp)
  var_imp_order <- apply(var_imp, 2, median) %>% order() %>% rev()
  var_imp <- var_imp[,var_imp_order]

  # --- 2.4. Compute the cumulative variable importance
  MODEL[["MBTR"]][["eval"]][["CUM_VIP"]] <- apply(var_imp[,1:3], 2, mean) %>% sum()

  # --- 2.5. Save raw vip
  MODEL[["MBTR"]][["vip"]] <- var_imp

  # --- 3. Removing low quality algorithm
  # In this case the only algorithm available...
  for(i in MODEL$MODEL_LIST){
    # --- 3.1. Based on model performance
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
    pdf(PDF_PATH)
    # Define the color (green = QC passed; red = no)
    if(length(MODEL$MODEL_LIST) == 1){pal <- "#1F867B"
    } else {pal <- "#B64A60"}

    # Simple barplot of the MBTR vip if the model passed QC
    par(mar = c(8,8,3,12))
    boxplot(var_imp, axes = FALSE, horizontal = TRUE, cex.names = 0.5,
            main = "PREDICTOR IMPORTANCE ( MBTR )",
            sub = paste("Predictive performance (R2) =", round(MODEL[["MBTR"]]$eval$R2, 2), "; Cumulated var. importance (%; top 3) =", round(MODEL[["MBTR"]]$eval$CUM_VIP, 0)),
            col = pal, ylab = "", xlab = "Variable importance (%)", las = 2, cex.axis = 0.6)
    axis(side = 4, at = 1:ncol(var_imp), labels = colnames(var_imp), las = 2, cex.axis = 0.6)
    axis(side = 1, at = seq(0, 100, 10), labels = seq(0, 100, 10), las = 1)
    abline(v = seq(0, 100, 10), lty = "dotted")
    box()
    box("figure", col = "black", lwd = 1)
    dev.off()
  }

  # --- 5. Wrap up and save
  return(MODEL)

} # END FUNCTION
