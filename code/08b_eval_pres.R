#' =============================================================================
#' @name eval_pres
#' @description sub-pipeline for model evaluation corresponding to presence data
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @return the MODELS object updated with evaluation metric values

eval_pres <- function(QUERY,
                      MODEL){
  
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 1. Load final model data 
    final_fit <- MODEL[[i]][["final_fit"]] %>% 
      collect_predictions()
    
    # --- 2. Extract observations and predictions
    y <- final_fit$measurementvalue
    y_hat <- final_fit$.pred
    
    # --- 3. Compute Continuous Boyce Index into MODELS object
    MODEL[[i]][["eval"]][["CBI"]] <- ecospat.boyce(fit = y_hat,
                                                    obs = y_hat[which(y == 1)],
                                                    PEplot = FALSE) %>% 
      .$cor
    
    # --- 4. Removing model from list if low quality fit
    # Fixed at 0.3 for CBI value or NA (in case of a 0 & 1 binary model prediction)
    if(MODEL[[i]][["eval"]][["CBI"]] < 0.3 | is.na(MODEL[[i]][["eval"]][["CBI"]])){
      MODEL$CALL$MODEL_LIST <- MODEL$CALL$MODEL_LIST[MODELS$CALL$MODEL_LIST != i]
      message(paste("--- EVAL : discarded", i, "due to CBI =", MODEL[[i]][["eval"]][["CBI"]], "< 0.3 \n"))
    }

  } # for each model loop
  
  return(MODEL)
  
} # END FUNCTION