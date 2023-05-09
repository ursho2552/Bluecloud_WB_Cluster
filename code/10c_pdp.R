#' =============================================================================
#' @name pdp
#' @description Compute and plot partial dependency plots at the model and
#' ensemble level. The proportion data have their own PDP function embedded in
#' MBTR. (FOR NOW)
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @param ENSEMBLE if TRUE, computes the ensemble pdp ?

pdp <- function(SP_SELECT = NULL,
                FOLDER_NAME = NULL,
                ENSEMBLE = FALSE,
                N_BOOTSTRAP = 10){
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"))
  
  # ============================= BUILDING PDPs ================================
  # With PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SP_SELECT,"/pdp.pdf"))
  
  # --- 1. Define bootstraps
  # --- 1.1. Re-assemble all query tables
  tmp <- cbind(QUERY$Y, QUERY$X, QUERY$S)
  
  # --- 1.2. Run the bootstrap generation from tidy models
  boot_split <- bootstraps(tmp, times = N_BOOTSTRAP,
                           strata = measurementvalue)
  
  # --- 1.3. Specify target and features names
  features <- QUERY$FOLDS$train %>% 
    dplyr::select(all_of(names(QUERY$X)))
  target <- QUERY$FOLDS$train %>% 
    dplyr::select(measurementvalue)
  
  # --- 1.4. Initialize global PDP storage
  pdp_all <- NULL
  
  # ======================= MODEL/BOOTSTRAP LOOP SECTION =======================
  for(i in MODEL$CALL$MODEL_LIST){
    # --- 2. Fit the bootstrap
    # --- 2.1. Register parallel
    # Only if the number of species is less then the number of available clusters
    # Otherwise, the parallel computing is done by species
    if(length(CALL$SP_SELECT) < LOCAL_CLUSTERS){
      cl <- makePSOCKcluster(LOCAL_CLUSTERS)
      doParallel::registerDoParallel(cl)
      message(paste(Sys.time(), "--- Parallel bootstrap for", i, ": START"))
    }
    # --- 2.2. Fit model on bootstrap and save fit object
    # We use control_resample > extract_fit_parsnip to save the fit object
    boot_fit <- MODEL[[i]][["final_wf"]] %>% 
      fit_resamples(resamples = boot_split,
                    control = control_resamples(extract = function (x) extract_fit_parsnip(x)),
                    allow_par = TRUE,
                    parallel_over = "everything") %>% 
      unnest(.extracts)
    
    # --- 2.3. Stop parallel backend
    stopCluster(cl)
    message(paste(Sys.time(), "--- Parallel grid tuning for", i, ": DONE"))
    
    # --- 2.4. Prepare the pdp object
    pdp_m <- NULL
    message(paste(Sys.time(), "--- PDP : computing for", i))
    
    for(j in 1:N_BOOTSTRAP){
      # --- 3. Build model explainer
      explainer <- explain_tidymodels(model = boot_fit$.extracts[[j]],
                                      data = features,
                                      y = target,
                                      verbose = FALSE)
      
      # --- 4. Compute the raw partial dependency plot
      # --- 4.1. Compute
      tmp <- model_profile(explainer = explainer,
                           variables = names(features),
                           N = NULL,
                           type = "partial") 
      tmp <- tmp[["agr_profiles"]] %>% 
        dplyr::select(-"_label_", -"_ids_")
      
      # --- 4.2. Save in the model-level object
      names(tmp) <- c("var","x","yhat") # Clean naming
      pdp_m <- rbind(pdp_m, tmp)
    } # End j bootstrap loop
    
    # --- 5. Compute mean and coefficient of variation
    # --- 5.1. Compute
    pdp_m <- pdp_m %>% 
      group_by(var, x) %>% 
      summarise(y_hat_m = mean(yhat),
                y_hat_cv = cv(yhat))
    
    # --- 5.2. Save in an all_pdp object
    pdp_all[[i]] <- pdp_m
  } # End model loop
  # ============================== END LOOP SECTION ============================
  
  # --- 6. Build ensemble pdp
  # pdp_all[["ENSEMBLE"]] <- reduce(rbind, pdp_all) %>%
  #   group_by()


  # --- 7. Plotting PDPs
  # --- 7.1. Initial plot settings
  par(mfrow = c(4,4), mar = c(2,2,4,5))
  pal <- c(brewer.pal(length(MODEL$CALL$MODEL_LIST), "Spectral"), "black")
  
  # --- 7.2. Iteratively compute the plots itself
  for(i in 1:length(features)){
    for(j in 1:length(MODEL$CALL$MODEL_LIST)){
      # --- 7.2.1. Prepare data table
      tmp <- pdp_all[[MODEL$CALL$MODEL_LIST[j]]] %>% 
        dplyr::filter(var == names(features)[i])
      
      # --- 7.2.2. Plot
      if(j == 1){
        plot(tmp$x, tmp$y_hat_m, type = 'l', ylim = c(0,1), lwd = 1, col = pal[j],
             xlab = "", ylab = "", main = names(features)[i])
        polygon(x = c(tmp$x, rev(tmp$x)),
                y = c(tmp$y_hat_m-tmp$y_hat_m*tmp$y_hat_cv/100, rev(tmp$y_hat_m+tmp$y_hat_m*tmp$y_hat_cv/100)),
                col = alpha(pal[j], 0.3), border = NA)
        mtext(side = 4, at = tail(tmp$y_hat_m, 1), text = MODEL$CALL$MODEL_LIST[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
        grid(col = "gray20")
      } else {
        lines(tmp$x, tmp$y_hat_m, type = 'l', ylim = c(0,1), lwd = 1, col = pal[j],
             xlab = "", ylab = "", main = names(features)[i])
        polygon(x = c(tmp$x, rev(tmp$x)),
                y = c(tmp$y_hat_m-tmp$y_hat_m*tmp$y_hat_cv/100, rev(tmp$y_hat_m+tmp$y_hat_m*tmp$y_hat_cv/100)),
                col = alpha(pal[j], 0.3), border = NA)
        mtext(side = 4, at = tail(tmp$y_hat_m, 1), text = MODEL$CALL$MODEL_LIST[j], col = pal[j], padj = 0.5, las = 1, cex = 0.6)
      } # End if
    } # End j model loop
  } # End i plot loop
  
  # Stop pdf saving
  dev.off()
  
} # END FUNCTION

