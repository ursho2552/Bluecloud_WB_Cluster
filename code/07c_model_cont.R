#' =============================================================================
#' @name model_cont
#' @description sub-pipeline corresponding to the model fitting procedure for
#' continuous data. Called via the model_wrapper.R function, thus no default
#' parameter values are defined.
#' @param QUERY the query object from the master pipeline
#' @param HP the hyperparameter object from the master pipeline. In case of no
#' hyperparameter selection, please pass an empty list as it is returned from
#' the hyperparameter function.
#' @param MODEL_LIST vector of string model names. It can be different from the
#' list passed to the hyperparameter function in the previous step
#' @return returns a list object containing the best model, associated hyper
#' parameters and predicted values per re sampling folds

model_cont <- function(QUERY,
                       HP,
                       MODEL_LIST){
  
  # --- 1. Define formula common to the model workflows
  tmp <- QUERY$CALL$ENV_VAR %>% paste(collapse = " + ")
  formula <- paste0("measurementvalue ~ ", tmp) %>% as.formula()
  
  # --- 2. Loop with all selected models
  for(i in 1:length(MODEL_LIST)){
    # --- Display information
    message(paste(Sys.time(), "--- Start hyper parameter tuning for", MODEL_LIST[i], "---"))
    
    # --- 2.1. Define workflow adapted to hyper parameter tuning
    # /!\ THE FORMULA ARGUMENT IS NECESSARY FOR GAM: DOUBLE CHECK CONSEQUENCES
    # Needed as a pre-processor for GAMs.
    model_wf <- workflow() %>% 
      add_model(HP[[MODEL_LIST[i]]][["model_spec"]], formula = formula) %>% 
      add_formula(formula)
    
    # --- 2.2. Register parallel
    # Only if the number of species is less then the number of available clusters
    # Otherwise, the parallel computing is done by species
    if(length(QUERY$CALL$SP_SELECT) < LOCAL_CLUSTERS){
      cl <- makePSOCKcluster(LOCAL_CLUSTERS)
      doParallel::registerDoParallel(cl)
      message(paste(Sys.time(), "--- Parallel grid tuning for", MODEL_LIST[[i]], ": START"))
    }
    
    # --- 2.3. Run the model for each fold x (hyper parameter grid rows)
    # Runs hyper parameter tuning if a grid is present in the HP (e.g. no GLM tune)
    if(!is.null(HP[[MODEL_LIST[i]]][["model_grid"]])){
      model_res <- model_wf %>% 
        tune_grid(resamples = QUERY$FOLDS$resample_split,
                  grid = HP[[MODEL_LIST[i]]][["model_grid"]],
                  metrics = yardstick::metric_set(rmse),
                  control = control_grid(allow_par = TRUE,
                                         parallel_over = "everything"))
    } else {
      model_res <- model_wf %>% 
        fit_resamples(resamples = QUERY$FOLDS$resample_split)
    }
    
    # --- 2.4. Stop parallel backend
    stopCluster(cl)
    message(paste(Sys.time(), "--- Parallel grid tuning for", MODEL_LIST[[i]], ": DONE"))
    
    # --- 2.5. Select best hyper parameter set
    # Based on RMSE values per model run (rsq does not work with 0's)
    model_best <- model_res %>% 
      select_best("rmse")
    
    # --- 2.6. Define final workflow
    final_wf <- model_wf %>% 
      finalize_workflow(model_best)
    
    # --- 2.7. Run the model on the initial split
    # Save the workflow and model object in a list to be passed to further steps
    final_fit <- final_wf %>%
      last_fit(QUERY$FOLDS$init_split) 
    
    HP[[MODEL_LIST[i]]][["final_wf"]] <- final_wf
    HP[[MODEL_LIST[i]]][["final_fit"]] <- final_fit
    
    # --- Display information
    message(paste(Sys.time(), "--- DONE ---"))
  }
  
  return(HP)
  
} # END FUNCTION



