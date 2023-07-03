#' =============================================================================
#' @name pdp
#' @description Compute and plot partial dependency plots at the model and
#' ensemble level. The proportion data have their own PDP function embedded in
#' MBTR. (FOR NOW)
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param ENSEMBLE if TRUE, computes the ensemble pdp ?

pdp <- function(FOLDER_NAME = NULL,
                SUBFOLDER_NAME = NULL,
                ENSEMBLE = FALSE,
                N_BOOTSTRAP = 10){
  
  # --- 1. Initialize function
  # --- 1.1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 1.2. Create PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/pdp.pdf"))
  
  # --- 2. Define bootstraps
  # --- 2.1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X, QUERY$S)
  
  # --- 2.2. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = N_BOOTSTRAP,
                           strata = measurementvalue)
  
  # --- 2.3. Specify target and features names
  features <- QUERY$FOLDS$train %>% 
    dplyr::select(all_of(names(QUERY$X)))
  target <- QUERY$FOLDS$train %>% 
    dplyr::select(measurementvalue)
  
  # --- 2.4. Initialize global PDP storage
  pdp_all <- NULL
  
  # ----------------------------------------------------------------------------
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 3. Fit the bootstrap
    # --- 3.1. Fit model on bootstrap and save fit object
    # We use control_resample > extract_fit_parsnip to save the fit object
    boot_fit <- MODEL[[i]][["final_wf"]] %>% 
      fit_resamples(resamples = boot_split,
                    control = control_resamples(extract = function (x) extract_fit_parsnip(x))) %>% 
      unnest(.extracts)

    # --- 3.2. Prepare the pdp object
    pdp_m <- NULL
    message(paste(Sys.time(), "--- PDP : computing for", i))
    
    for(j in 1:N_BOOTSTRAP){
      # --- 4. Build model explainer
      explainer <- explain_tidymodels(model = boot_fit$.extracts[[j]],
                                      data = features,
                                      y = target,
                                      verbose = FALSE)
      
      # --- 5. Compute the raw partial dependency plot
      # --- 5.1. Compute
      tmp <- model_profile(explainer = explainer,
                           variables = names(features),
                           N = 100,
                           type = "partial") 
      tmp <- tmp[["agr_profiles"]] %>% 
        dplyr::select(-"_label_", -"_ids_")
      
      # --- 5.2. Save in the model-level object
      names(tmp) <- c("var","x","yhat") # Clean naming
      if(j == 1){pdp_m <- tmp} 
      else {pdp_m <- cbind(pdp_m, tmp$yhat)}
    } # End j bootstrap loop
    
    # --- 6. Compute mean and coefficient of variation
    # --- 6.1. Compute
    tmp_m <- apply(pdp_m[,-c(1,2)], 1, mean)
    tmp_cv <- apply(pdp_m[,-c(1,2)], 1, cv)
    pdp_m <- cbind(pdp_m[,c(1,2)], tmp_m, tmp_cv)
    names(pdp_m) <- c("var","x","y_hat_m","y_hat_cv")
    
    # --- 6.2. Save in an all_pdp object
    pdp_all[[i]] <- pdp_m
  } # End model loop
  # ----------------------------------------------------------------------------
  
  # --- 7. Build ensemble pdp
  # TO BE IMPLEMENTED


  # --- 8. Plotting PDPs
  # --- 8.1. Initial plot settings
  par(mfrow = c(4,4), mar = c(2,2,4,5))
  pal <- c(brewer.pal(length(MODEL$CALL$MODEL_LIST), "Spectral"), "black")
  
  # --- 8.2. Iteratively compute the plots itself
  for(i in names(features)){
    for(j in 1:length(MODEL$CALL$MODEL_LIST)){
      # --- 8.2.1. Prepare data table
      tmp <- pdp_all[[MODEL$CALL$MODEL_LIST[j]]] %>% 
        dplyr::filter(var == i)
      
      # --- 8.2.2. Plot
      if(j == 1){
        plot(tmp$x, tmp$y_hat_m, type = 'l', ylim = c(0,1), lwd = 1, col = pal[j],
             xlab = "", ylab = "", main = i)
        polygon(x = c(tmp$x, rev(tmp$x)),
                y = c(tmp$y_hat_m-tmp$y_hat_m*tmp$y_hat_cv/100, rev(tmp$y_hat_m+tmp$y_hat_m*tmp$y_hat_cv/100)),
                col = scales::alpha(pal[j], 0.3), border = NA)
        mtext(side = 4, at = tail(tmp$y_hat_m, 1), text = MODEL$CALL$MODEL_LIST[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
        grid(col = "gray20")
      } else {
        lines(tmp$x, tmp$y_hat_m, type = 'l', ylim = c(0,1), lwd = 1, col = pal[j],
             xlab = "", ylab = "", main = i)
        polygon(x = c(tmp$x, rev(tmp$x)),
                y = c(tmp$y_hat_m-tmp$y_hat_m*tmp$y_hat_cv/100, rev(tmp$y_hat_m+tmp$y_hat_m*tmp$y_hat_cv/100)),
                col = scales::alpha(pal[j], 0.3), border = NA)
        mtext(side = 4, at = tail(tmp$y_hat_m, 1), text = MODEL$CALL$MODEL_LIST[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
      } # End if
    } # End j model loop
  } # End i plot loop
  
  # --- 9. Wrap up and save
  # --- 9.1. Stop PDF saving
  dev.off()
  
} # END FUNCTION

