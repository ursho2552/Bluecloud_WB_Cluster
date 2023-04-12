#' =============================================================================
#' @name eval_cont
#' @description sub-pipeline for model evaluation corresponding to continuous data
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @return the MODELS object updated with evaluation metric values

eval_cont <- function(QUERY,
                      MODELS){
  
  for(i in MODELS$CALL$MODEL_LIST){
    # --- 1. Load final model data 
    final_fit <- MODELS[[i]][["final_fit"]] %>% 
      collect_predictions()
    
    # --- 2. Extract observations and predictions
    y <- final_fit$measurementvalue
    y_hat <- final_fit$.pred
    
    # --- 3. Compute R-squared into MODELS object
    df <- data.frame(truth = y, estimate = y_hat)
    MODELS[[i]][["eval"]][["R2"]] <- yardstick::rsq(data = df, truth, estimate ) %>% 
      .$.estimate %>% 
      round(3)
    
    # --- 4. Removing model from list if low quality fit
    # Fixed at 0.3 for R2 value
    if(MODELS[[i]][["eval"]][["R2"]] < 0.3){
      MODELS$CALL$MODEL_LIST <- MODELS$CALL$MODEL_LIST[MODELS$CALL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to R2 =", MODELS[[i]][["eval"]][["R2"]], "< 0.3 \n"))
    }
    
  } # for each model loop
  
  return(MODELS)
  
} # END FUNCTION