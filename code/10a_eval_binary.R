#' =============================================================================
#' @name eval_binary
#' @description sub-pipeline for model evaluation corresponding to binary data
#' @param QUERY the query object from the master pipeline
#' @param MODEL the models object from the master pipeline#'
#' @param ENSEMBLE if TRUE, computes variable importance metrics for the ensemble 
#' model as well
#' @return the MODEL object updated with evaluation metric values (model performance
#' metric and variable importance metric)
#' @return variable importance plots as PDF file

eval_binary <- function(QUERY,
                        MODEL,
                        ENSEMBLE = TRUE){

  # --- 1. Model performance assessment
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 1.1. Load final model data 
    final_fit <- MODEL[[i]][["final_fit"]] %>% 
      collect_predictions()
    
    # --- 1.2. Extract observations and predictions
    y <- final_fit$measurementvalue
    y_hat <- final_fit$.pred
    
    # --- 1.3. Compute Continuous Boyce Index into MODELS object
    MODEL[[i]][["eval"]][["CBI"]] <- ecospat.boyce(fit = y_hat,
                                                            obs = y_hat[which(y == 1)],
                                                            PEplot = FALSE,
                                                            method = "spearman",
                                                            rm.duplicate = TRUE,
                                                            res = 100) %>%
      .$cor
    
    # --- 1.4. Compute an over-fitting rate
    # Calculated as the deviation during training to deviation in testing ratio
    resample_dev <- MODEL[[i]][["best_fit"]]$mean
    eval_dev <- rmse(data = data.frame(truth = y, estimate = y_hat), 1, 2)$.estimate
    MODEL[[i]][["eval"]][["overfit_rate"]] <- ((resample_dev/eval_dev)-1)*100
    
  } # for each model loop
  
  # --- 2. Variable importance - algorithm level
  # --- 2.1. Initialize function
  # --- 2.1.1. Storage and graphical specification
  var_imp <- NULL
  par(mfrow = c(3,2), mar = c(5,5,5,5))
  
  # --- 2.1.2. Model related data
  features <- QUERY$FOLDS$train %>% 
    dplyr::select(all_of(QUERY$SUBFOLDER_INFO$ENV_VAR))
  target <- QUERY$FOLDS$train %>% 
    dplyr::select(measurementvalue)
 
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 2.2. Extract final model fit
    m <- extract_fit_parsnip(MODEL[[i]][["final_fit"]])
    
    # --- 2.3. Build model explainer
    explainer <- explain_tidymodels(model = m,
                                    data = features,
                                    y = target)
    
    # --- 2.4. Compute variable importance
    # --- 2.4.1. First in terms of RMSE, i.e., raw var importance for later ensemble computing
    message(paste("--- VAR IMPORTANCE : compute for", i))
    var_imp[[i]][["Raw"]] <- model_parts(explainer = explainer,
                                loss_function = loss_root_mean_square) %>% 
      dplyr::filter(permutation != 0) %>% 
      dplyr::filter(variable != "_baseline_") %>% 
      group_by(permutation) %>% 
      mutate(value = dropout_loss - dropout_loss[variable == "_full_model_"]) %>% 
      dplyr::filter(variable != "_full_model_")
    
    # --- 2.4.2. Further compute it as percentage for model-level plot
    # Security if a model did not fit, hence var_imp = 0 for all predictors
    # Avoids an error leading to a function stop, while other models could be OK
    if(sum(var_imp[[i]][["Raw"]][["value"]]) > 0){
      var_imp[[i]][["Percent"]] <- var_imp[[i]][["Raw"]] %>% 
        mutate(value = value / sum(value) * 100)  %>% 
        ungroup() %>% 
        dplyr::select(variable, value) %>% 
        mutate(variable = fct_reorder(variable, value, .desc = TRUE, .na_rm = TRUE))
    } else {
      var_imp[[i]][["Percent"]] <- var_imp[[i]][["Raw"]]
    }
    
    # --- 2.4.3. Compute cumulative variable importance
    MODEL[[i]][["eval"]][["CUM_VIP"]] <- var_imp[[i]][["Percent"]] %>% 
      group_by(variable) %>% 
      summarise(average = mean(value)) %>% 
      dplyr::slice(1:3) %>% 
      dplyr::select(average) %>% 
      sum()
    
  } # for each model loop
  
  # --- 3. Removing low quality algorithms
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 3.1. Based on model performance
    # Fixed at 0.3 for CBI value or NA (in case of a 0 & 1 binary model prediction)
    if(MODEL[[i]][["eval"]][["CBI"]] < 0.5 | is.na(MODEL[[i]][["eval"]][["CBI"]])){
      MODEL$CALL$MODEL_LIST <- MODEL$CALL$MODEL_LIST[MODEL$CALL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to CBI =", MODEL[[i]][["eval"]][["CBI"]], "< 0.5 \n"))
    }
    
    # --- 3.2. Based on cumulative variable importance
    # Fixed at 30% cumulative importance for the top three predictors
    if(MODEL[[i]][["eval"]][["CUM_VIP"]] < 50 | is.na(MODEL[[i]][["eval"]][["CUM_VIP"]])){
      MODEL$CALL$MODEL_LIST <- MODEL$CALL$MODEL_LIST[MODEL$CALL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to CUM_VIP =", MODEL[[i]][["eval"]][["CUM_VIP"]], "< 50% \n"))
    }
  } # for each model loop
  
  # --- 4. Variable importance - Ensemble level
  # Variable importance values also scaled by the CBI metric value
  if(ENSEMBLE == TRUE & length(MODEL$CALL$MODEL_LIST > 1)){
    # --- 4.1. Aggregate and weight raw data
    ens_imp <- NULL
    message("--- VAR IMPORTANCE : compute ensemble")
    
    for(i in MODEL$CALL$MODEL_LIST){
      # Concatenate eval-weighted raw variable importance
      ens_imp <- rbind(ens_imp, var_imp[[i]][["Raw"]])
    } # End i model loop
    
    # --- 4.2. Further compute it as percentage
    var_imp[["ENSEMBLE"]][["Percent"]] <- ens_imp %>% 
      group_by(permutation) %>% 
      mutate(value = value / sum(value) * 100 * length(MODEL$CALL$MODEL_LIST))  %>% 
      ungroup() %>% 
      dplyr::select(variable, value) %>% 
      mutate(variable = fct_reorder(variable, value, .desc = TRUE))
  } # End if ensemble TRUE
  
  # --- 5. Variable importance - Plot
  # Plot algorithm level and ensemble variable importance
  if(length(MODEL$CALL$MODEL_LIST > 1)){
    for(i in MODEL$CALL$MODEL_LIST){
      tmp <- var_imp[[i]][["Percent"]]
      boxplot(tmp$value ~ tmp$variable, axes = FALSE, 
              main = paste("Model-level for :", i), col = "gray50",
              xlab = "", ylab = "Variable importance (%)")
      axis(side = 1, at = 1:ncol(features), labels = levels(tmp$variable), las = 2)
      axis(side = 2, at = seq(0, 100, 10), labels = seq(0, 100, 10), las = 2)
      abline(h = seq(0, 100, 10), lty = "dotted")
      box()
    } # End i model loop
    if(ENSEMBLE == TRUE){
      tmp <- var_imp[["ENSEMBLE"]][["Percent"]] 
      boxplot(tmp$value ~ tmp$variable, axes = FALSE, 
              main = paste("Ensemble"), col = "gray50",
              xlab = "", ylab = "Variable importance (%)")
      axis(side = 1, at = 1:ncol(features), labels = levels(tmp$variable), las = 2)
      axis(side = 2, at = seq(0, 100, 10), labels = seq(0, 100, 10), las = 2)
      abline(h = seq(0, 100, 10), lty = "dotted")
      box()
    } # End if ensemble TRUE
  } # End if model list > 1
  
  # --- 6. Wrap up and save
  return(MODEL)
  
} # END FUNCTION