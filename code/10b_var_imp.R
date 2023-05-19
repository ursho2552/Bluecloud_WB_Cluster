#' =============================================================================
#' @name var_imp
#' @description Function to compute and extract the variable importance related
#' to model fitting. The proportion data variable importance is already embedded
#' in the MBTR model. Thus, this function is only related to pres or cont data.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param ENSEMBLE if TRUE, computes the variable importance as an ensemble ? 
#' @return plots mean and uncertainty maps per model or ensemble
#' @return the cumulative variable importance saved in the MODEL object. Used as
#' a quality assessment of the model later

var_imp <- function(FOLDER_NAME = NULL,
                    SUBFOLDER_NAME = NULL,
                    ENSEMBLE = FALSE){
  
  # --- 1. Initialize function
  # --- 1.1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 1.2. Create PDf saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/variable_importance.pdf"))
  
  # --- 1.3. Set initial parameters
  # --- 1.3.1. Model related data
  features <- QUERY$FOLDS$train %>% 
    dplyr::select(all_of(names(QUERY$X)))
  target <- QUERY$FOLDS$train %>% 
    dplyr::select(measurementvalue)
  
  # --- 1.3.2. Storage and graphical specification
  var_imp <- NULL
  par(mfrow = c(3,2), mar = c(5,5,5,5))
  
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 2. Extract final model fit
    m <- extract_fit_parsnip(MODEL[[i]][["final_fit"]])
    
    # --- 3. Build model explainer
    explainer <- explain_tidymodels(model = m,
                                    data = features,
                                    y = target)
    
    # --- 4. Compute variable importance
    # --- 4.1. First in terms of RMSE, i.e., raw var importance for later ensemble computing
    message(paste("--- VAR IMPORTANCE : compute for", i))
    var_imp[[i]] <- model_parts(explainer = explainer,
                                loss_function = loss_root_mean_square) %>% 
      dplyr::filter(permutation != 0) %>% 
      dplyr::filter(variable != "_baseline_") %>% 
      group_by(permutation) %>% 
      mutate(value = dropout_loss - dropout_loss[variable == "_full_model_"]) %>% 
      dplyr::filter(variable != "_full_model_")
    
    # --- 4.2. Further compute it as percentage for model-level plot
    tmp <- var_imp[[i]] %>% 
      mutate(value = value / sum(value) * 100)  %>% 
      ungroup() %>% 
      dplyr::select(variable, value) %>% 
      mutate(variable = fct_reorder(variable, value, .desc = TRUE))
    
    # --- 4.3 Plot the variable importance
    boxplot(tmp$value ~ tmp$variable, axes = FALSE, 
            main = paste("Model-level for :", i), col = "gray50",
            xlab = "", ylab = "Variable importance (%)")
    axis(side = 1, at = 1:ncol(features), labels = levels(tmp$variable), las = 2)
    axis(side = 2, at = seq(0, 100, 10), labels = seq(0, 100, 10), las = 2)
    abline(h = seq(0, 100, 10), lty = "dotted")
    box()
    
  } # End i model loop
  
  
  if(ENSEMBLE == TRUE & length(MODEL$CALL$MODEL_LIST > 1)){
    # --- 5. Compute an ensemble plot
    # Variable importance values also scaled by the metric value
    
    # --- 5.1. Aggregate and weight raw data
    ens_imp <- NULL
    message("--- VAR IMPORTANCE : compute ensemble")
    
    for(i in MODEL$CALL$MODEL_LIST){
      # Concatenate eval-weighted raw variable importance
      tmp <- var_imp[[i]] %>% 
        mutate(value = value * MODEL[[i]][["eval"]][[1]])
      ens_imp <- rbind(ens_imp, tmp)
    } # End i model loop
    
    # --- 5.2. Further compute it as percentage
    tmp <- ens_imp %>% 
      group_by(permutation) %>% 
      mutate(value = value / sum(value) * 100 * length(MODEL$CALL$MODEL_LIST))  %>% 
      ungroup() %>% 
      dplyr::select(variable, value) %>% 
      mutate(variable = fct_reorder(variable, value, .desc = TRUE))
      
    # --- 5.3 Plot the variable importance
    boxplot(tmp$value ~ tmp$variable, axes = FALSE, 
            main = paste("Ensemble"), col = "gray50",
            xlab = "", ylab = "Variable importance (%)")
    axis(side = 1, at = 1:ncol(features), labels = levels(tmp$variable), las = 2)
    axis(side = 2, at = seq(0, 100, 10), labels = seq(0, 100, 10), las = 2)
    abline(h = seq(0, 100, 10), lty = "dotted")
    box()
    
  } # End if ensemble TRUE

  # --- 6. Compute cumulative variable importance
  # --- 6.1. Compute
  CUM_VIP <- tmp %>% 
    group_by(variable) %>% 
    summarise(average = mean(value)) %>% 
    dplyr::slice(1:3) %>% 
    dplyr::select(average) %>% 
    sum()
  # --- 6.2. Append
  MODEL$ENSEMBLE$EVAL$CUM_VIP <- CUM_VIP
  
  # --- 7. Wrap up and save
  # --- 7.1. Stop PDF saving
  dev.off()
  
  # --- 7.2. Save file(s)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"),
        compress = "gzip", compression_level = 6)
  
} # END FUNCTION
  
  