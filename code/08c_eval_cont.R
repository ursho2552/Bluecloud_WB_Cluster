#' =============================================================================
#' @name eval_cont
#' @description sub-pipeline for model evaluation corresponding to continuous data
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @return the MODELS object updated with evaluation metric values

eval_cont <- function(QUERY,
                      MODEL){
  
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 1. Load final model data 
    final_fit <- MODEL[[i]][["final_fit"]] %>% 
      collect_predictions()
    
    # --- 2. Extract observations and predictions
    y <- final_fit$measurementvalue
    y_hat <- final_fit$.pred
    
    # --- 3. Compute R-squared into MODELS object
    df <- data.frame(truth = y, estimate = y_hat)
    MODEL[[i]][["eval"]][["R2"]] <- yardstick::rsq(data = df, truth, estimate ) %>% 
      .$.estimate %>% 
      round(3)
    
    # --- 4. Removing model from list if low quality fit
    # Fixed at 0.3 for R2 value
    if(MODEL[[i]][["eval"]][["R2"]] < 0.3){
      MODEL$CALL$MODEL_LIST <- MODEL$CALL$MODEL_LIST[MODEL$CALL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to R2 =", MODEL[[i]][["eval"]][["R2"]], "< 0.3 \n"))
    }
    
  } # for each model loop
  
  return(MODEL)
  
} # END FUNCTION