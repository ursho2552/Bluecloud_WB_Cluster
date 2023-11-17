#' =============================================================================
#' @name diversity_maps
#' @description Function extracting all projections by species and bootstrap.
#' Then computes a set of diversity metrics and plots
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param N_BOOTSTRAP the number of bootstrap in the projection step
#' @param BETA the list of beta diversity dissimilarity index. Among those
#' available in vegan::vegdist()
#' @param BUFFER the radius of neighboring cells to consider for beta diversity
#' @return plots mean and uncertainty maps per diversity metric
#' @return a diversity and mess object per projection, species and bootstrap.
#' Saved in FOLDERNAME.

diversity_maps <- function(FOLDER_NAME = NULL,
                           SUBFOLDER_NAME = NULL,
                           MONTH = list(c(10,11,12,1,2,3),
                                        4:9)){
  
  # --- 1. Initialize function
  # --- 1.1. Global parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 1.2. Baseline raster
  r0 <- CALL$ENV_DATA[[1]][[1]]
  
  # --- 1.3. Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 1.4. Beta diversity function
  beta_div <- function(ID, NX, NY, VALUE, BUFFER, BETA){
    # --- Get cell identifiers
    cells <- get_cell_neighbors(NX = NX,
                                NY = NY, 
                                ID = ID,
                                BUFFER = BUFFER)
    # --- Compute average pairwise dissimilarity
    diss <- vegan::vegdist(VALUE[c(cells$center_id, cells$neighbors_id),,,], method = BETA, na.rm = TRUE) %>% 
      as.matrix() %>% 
      .[2:nrow(.), 1] %>% 
      mean(na.rm = TRUE)
    return(diss)
  }
  
  # --- 2. Extract ensembles by species
  # --- 2.1. For binary or continuous data
  if(CALL$DATA_TYPE != "proportions"){
    all_ens <- mclapply(SUBFOLDER_NAME, FUN = function(s){
      # --- 2.1.1. Load subfolder information
      # We have a file exist security here in case the diversity is run later
      f <- paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/MODEL.RData")
      if(file.exists(f)){load(f)} else {MODEL$MODEL_LIST <- NULL}
      
      # --- 2.1.2. Early return if there is not enough algorithms passing QC
      if(length(MODEL$MODEL_LIST) == 0){
        s_ens <- NULL
      } # if early return
      
      # --- 2.1.3. Compute ensemble across algorithms
      if(length(MODEL$MODEL_LIST) != 0){
        s_ens <- lapply(MODEL$MODEL_LIST, FUN = function(m){
          MODEL[[m]]$proj$y_hat
        }) %>% abind(along = 4) %>% apply(c(1,2,3), mean)
      } # compute ensemble
    },
    mc.cores = MAX_CLUSTERS) %>% abind(along = 4) %>% aperm(c(1,4,2,3))
  } # if binary or continuous
  
  # --- 2.2. For proportions
  if(CALL$DATA_TYPE == "proportions"){
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/MODEL.RData"))
    s_ens <- MODEL$MBTR$proj$y_hat %>% aperm(c(1,3,2,4))
  } # if proportions
  
  # --- 2.3. Return if not enough ensembles
  if(dim(all_ens)[[2]] < 5){
    message("--- DIVERSITY : not enough ensembles to compute diversity from")
  }
  
  # --- 3. Compute alpha diversity indices
  a_shannon <- apply(all_ens, c(1,3,4), function(x)(x = vegan::diversity(x, "shannon")))
  a_richness <- apply(all_ens, c(1,3,4), function(x)(x = sum(x, na.rm = TRUE)))
  a_evenness <- a_shannon/(a_richness)
  a_invsimpson <- apply(all_ens, c(1,3,4), function(x)(x = vegan::diversity(x, "invsimpson")))
  
  # --- 4. Compute beta diversity indices
  # To add
  
  # --- 5. Stack diversities
  div_all <- abind(a_shannon, a_richness, a_evenness, a_invsimpson, along = 4)
  div_m <- apply(div_all, c(1,3,4), mean)
  div_sd <- apply(div_all, c(1,3,4), sd)
  div_names <- c("a_shannon", "a_richness", "a_evenness", "a_invsimpson")
  
  # --- 6. Extract MESS analysis
  message(paste(Sys.time(), "--- Extract MESS"))
  mess_all <- lapply(SUBFOLDER_NAME, FUN = function(s){
    # --- 6.1. Load corresponding objects
    # We have a file exist security here in case the diversity is run later
    f <- paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/MODEL.RData")
    if(file.exists(f)){
      load(paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/QUERY.RData"))
      load(paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/MODEL.RData"))
      
      # --- 6.2. Extract mess
      if(length(MODEL$MODEL_LIST) != 0){
        mess_s <- getValues(QUERY$MESS*-1)
      } # if ensemble
    } # if exist
  }) %>% abind(along = 3) %>% apply(c(1,2), mean)
  
  # --- 6.3. Final adjustments
  mess_all[mess_all > 100] <- 100
  
  # --- 7. Diversity plots
  # --- 7.1. Initialize pdf and plot
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/diversity_maps.pdf"))
  par(mfrow = c(3,2), mar = c(2,3,3,3))
  
  # --- 7.2. Plot the legends
  # --- 7.2.1. Diversity legend
  plot.new()
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.1,
                col = inferno_pal(100), border = "black")
  axis(side = 1)
  text(x = 0.5, y = 0.3, "Diversity value [0 - 1]", adj = 0.5)
  
  # --- 7.2.2. MESS x CV 2D color scale
  par(mar = c(5,5,1,2))
  bivar_pal <- colmat(nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Standard deviation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = round(seq(0, 0.25, length.out = 6), 2))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(5,3,3,1))
  
  # --- 7.3. Plot diversity maps
  # Loop over diversity maps and projections
  for(i in 1:dim(div_m)[[3]]){
    for(m in MONTH){
      
      # --- 7.3.1. Prepare the data
      tmp_m <- div_m[,m,] %>% apply(c(1,3), mean) %>% .[,i]
      tmp_sd <- div_sd[,m,] %>% apply(c(1,3), mean) %>% .[,i]
      plot_scale <- quantile(tmp_m, 0.95, na.rm = TRUE)
      
      # --- 7.3.2. Assign to raster
      # Average raster is cut at plot scale (Q95)
      r_m <- r0 %>% setValues(tmp_m)
      r_m[r_m > plot_scale] <- plot_scale
      
      r_sd <- r0 %>% setValues(tmp_sd)
      
      r_mess <- r0 %>% setValues(apply(mess_all[, m], 1, mean))
      r_mess[r_mess<0] <- 1e-10 # temporary fix : modify bivarmap so that 0 is included
      r_mess[r_mess>100] <- 100
      
      # --- 7.3.3. Average projection
      plot(r_m, col = inferno_pal(100), legend=FALSE,
           main = paste("DIVERSITY (",div_names[i], ") \n Month:", paste(m, collapse = ",")), cex.main = 1)
      mtext(text = paste("Projection rescaling factor: ", round(plot_scale,2), 
                         "\n Number of species ensembles:", dim(all_ens)[[2]]),
            side = 1, line = 3, cex = 0.7)
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      plot(raster::rasterToContour(r_m, nlevels = 4), add = TRUE)
      box("figure", col="black", lwd = 1)
      
      # --- 7.3.4. Uncertainties projection
      # First we draw a land mask
      plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties", cex.main = 1)
      # Then compute and plot the 2-dimensional color scale for SD x MESS
      r <- bivar_map(rasterx = r_sd, rastery = r_mess, colormatrix = bivar_pal,
                     cutx = seq(0,max(getValues(r_sd), na.rm = TRUE), length.out = 101), cuty = 0:100)
      plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
      plot(raster::rasterToContour(r_sd, nlevels = 4), add = TRUE)
      box("figure", col="black", lwd = 1)
      
    } # m month
  } # i diversity
  
  # --- 8. Stop PDF
  dev.off()
  
} # END FUNCTION


