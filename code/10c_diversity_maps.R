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
                           N_BOOTSTRAP = 10,
                           BETA = c("hellinger","bray","altGower"),
                           BUFFER = 1){
  
  # --- 1. Initialize function
  # --- 1.1. Global parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 1.2. Baseline raster
  r0 <- raster(paste0(project_wd, "/data/features_mean_from_monthly"))
  
  # --- 2. Extract and compute diversity indices
  # --- 2.1. Create global object
  div_all <- NULL
  div_names <- c(paste("Alpha [", c("shannon","richness","evenness","invsimpson"), "]"), 
                 paste("Beta [", BETA, "]"))
  
  # --- 2.2. Loop over BOOTSTRAP
  # For memory efficiency - computing time trade off
  for(b in 1:N_BOOTSTRAP){
    message(paste(Sys.time(), "--- Computing diversity indices for bootstrap n°", b))
    ens_by_b <- NULL
    
    # --- 2.2.1. Loop over SUBFOLDER_NAME
    message(paste(Sys.time(), "--- Extract ensembles"))
    for(s in SUBFOLDER_NAME){
      # --- 2.2.1.1. Load MODEL object
      load(paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/MODEL.RData"))
      
      # --- 2.2.1.2. Compute ensemble across each algorithm
      ens_by_s <- NULL
      if(length(MODEL$CALL$MODEL_LIST != 0)){
        for(m in MODEL$CALL$MODEL_LIST){
          tmp <- MODEL[[m]][["proj"]][["y_hat"]][,b]*MODEL[[m]][["eval"]][[1]] # Take care for abundances (sum at 1)
          ens_by_s <- cbind(ens_by_s, tmp)
        } # m model_list loop
        ens_by_s <- ens_by_s/max(ens_by_s, na.rm = TRUE) # Divide by the sum of metrics, not by max (for abundances ?)
        ens_by_s <- apply(ens_by_s, 1, function(x)(x = mean(x, na.rm = TRUE)))
      } # if
      ens_by_b <- cbind(ens_by_b, ens_by_s)
    } # s sub folder name loop
    ens_by_b[ens_by_b < 0] <- 0 # Security
    
    # --- 2.2.2. Compute diversity indices
    message(paste(Sys.time(), "--- Compute diversity"))
    # --- 2.2.2.1. Alpha diversity
    a_shannon <- apply(ens_by_b, 1, function(x)(x = vegan::diversity(x, "shannon")))
    a_richness <- apply(ens_by_b, 1, function(x)(x = sum(x, na.rm = TRUE)))
    a_evenness <- a_shannon/(a_richness)
    a_invsimpson <- apply(ens_by_b, 1, function(x)(x = vegan::diversity(x, "invsimpson")))
    
    div_by_b <- cbind(a_shannon, a_richness, a_evenness, a_invsimpson)
    
    # --- 2.2.2.2. Beta diversity
    # --- 2.2.2.2.1. Define function
    beta_div <- function(ID, NX, NY, VALUE, BUFFER, BETA){
      # --- Get cell identifiers
      cells <- get_cell_neighbors(NX = NX,
                                  NY = NY, 
                                  ID = ID,
                                  BUFFER = BUFFER)
      # --- Compute average pairwise dissimilarity
      diss <- vegan::vegdist(VALUE[c(cells$center_id, cells$neighbors_id),], method = BETA, na.rm = TRUE) %>% 
        as.matrix() %>% 
        .[2:nrow(.), 1] %>% 
        mean(na.rm = TRUE)
      return(diss)
    }
    
    # --- 2.2.2.2.2. List of cells
    id <- 1:dim(ens_by_b)[[1]]
    # --- 2.2.2.2.3. Compute indices
    for(i in BETA){
      tmp <- mclapply(id, function(x) beta_div(ID = x, NX = ncol(r0), NY = nrow(r0), BUFFER = BUFFER, 
                                               VALUE = ens_by_b, BETA = i),
                      mc.cores = MAX_CLUSTERS) %>% 
        unlist()
      div_by_b <- cbind(div_by_b, tmp)
    } # i beta diversity loop
    
    # --- 2.2.3. Store diversity in global object
    div_all <- abind(div_all, div_by_b, along = 3)
    
  } # b bootstrap loop
  
  # --- 2.4. Save the diversity object
  save(div_all, file = paste0(project_wd, "/output/", FOLDER_NAME, "/DIVERSITY.RData"),
       compress = "gzip", compression_level = 6)
  
  # --- 3. Extract MESS analysis
  # --- 3.1. Create global object
  mess_all <- NULL
  
  # --- 3.2. Loop over SUBFOLDER_NAME
  message(paste(Sys.time(), "--- Extract MESS"))
  for(s in SUBFOLDER_NAME){
    # --- 3.2.1. Load QUERY and MODEL object
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/QUERY.RData"))
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", s, "/MODEL.RData"))
    
    # --- 3.2.2. Extract MESS
    if(length(MODEL$CALL$MODEL_LIST != 0)){
      mess_all <- cbind(mess_all, getValues(QUERY$MESS*-1))
    } # if
  } # s sub folder name loop
  
  # --- 4. Diversity plots
  # --- 4.1. Compute necessary objects
  # --- 4.1.1. Land mask
  land <- r0
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 4.1.2. Average MESS
  mess_all[mess_all > 100] <- 100
  r_mess <- r0 %>% 
    setValues(apply(mess_all, 1, function(x)(x = mean(x, na.rm = TRUE))))
  r_mess[r_mess<0] <- NA
  
  # --- 4.1.3. Average diversity
  div_m <- apply(div_all, c(1,2), function(x)(x = mean(x, na.rm = TRUE)))
  
  # --- 4.1.3. Average diversity
  div_cv <- apply(div_all, c(1,2), function(x)(x = cv(x, na.rm = TRUE)))
  div_cv[div_cv > 100] <- 100
  
  # --- 4.2. Start PDF saving
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/diversity_maps.pdf"))
  
  # --- 4.3. Plot the legends
  par(mfrow = c(4,2), mar = c(2,3,3,3))
  # --- 4.3.1. Diversity legend
  hsi_pal <- inferno_pal(100)
  plot.new()
  colorbar.plot(x = 0.5, y = 0, strip = seq(0,1,length.out = 100),
                strip.width = 0.3, strip.length = 2.1,
                col = hsi_pal, border = "black")
  axis(side = 1)
  text(x = 0.5, y = 0.3, "Diversity value [0 - 1]", adj = 0.5)
  # --- 4.3.2. MESS x CV 2D color scale
  par(mar = c(5,5,1,2))
  bivar_pal <- colmat(xmax = "deepskyblue4", ymax = "darkgoldenrod2", nbreaks = 100)
  colmat_plot(bivar_pal, xlab = "Coefficient of variation", ylab = "MESS value")
  axis(side = 1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
  axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, -20, -40, -60, -80, -100), las = 2)
  par(mar = c(1,3,3,1))
  
  # --- 4.4. Plot the diversity maps
  for(i in 1:dim(div_all)[[2]]){
    # --- 4.4.1. Average map
    r_m <- r0 %>% setValues(div_m[,i])
    plot(r_m, col = hsi_pal, legend = FALSE,  main = div_names[i])
    plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
    
    # --- 4.4.2. Uncertainties map
    r_cv <- r0 %>% setValues(div_cv[,i])
    plot(land, col = "antiquewhite4", legend=FALSE, main = "Uncertainties")
    r <- bivar_map(rasterx = r_cv, rastery = r_mess, colormatrix = bivar_pal,
                   cutx = 0:100, cuty = 0:100)
    plot(r[[1]], col = r[[2]], legend=FALSE, add = TRUE)
  } # i diversity loop
  
  # --- 4.5. Stop PDF
  dev.off()
  
} # END FUNCTION


