#' =============================================================================
#' @name generate_virtual
#' @description generates a set of virtual species for a given data type
#' @param ENV_VAR a list of .nc files to extract the main variable from, located in ENV_PATH
#' @param ENV_PATH string or vector of path to the root where the .nc are.
#' @param DATA_TYPE the output type of the data, which can influence the sub-folder
#' architecture. See details in run_init function
#' @param MONTH integer of the month on which to perform the initial species. An explicit
#' monthly model is not necessary as the test is carried out in the environmental space
#' @param MEANS a list of vectors of length 2, giving the mean coordinate of the species
#' in the environmental PCA space
#' @param SDS a list of vectors of length 2, giving the coordinate sd of the species
#' in the environmental PCA space
#' @return a CSV file containing the sampled data. To be sourced as custom DATA_SOURCE
#' in list_bio function
#' @return a VIRTUAL object in the memory, containing all metadata and virtual
#' species distribution. Can be copied in the run folder after list_bio.

generate_virtual <- function(ENV_VAR = c("!dist2coast_allmonths"),
                             ENV_PATH = "/net/meso/work/nknecht/Masterarbeit/General_Pipeline/Data/environmental_climatologies",
                             MONTH = 6,
                             SP_NB = 20,
                             N = c(50, 200, 1000),
                             NOISE_SD = c(1,2,3),
                             DATA_TYPE = "binary"){

  # --- 1. Prepare the raster stack
  # It has to be the same as the one used in the master script later (i.e; variables and month)
  # For the sake of simplicity, we only produce virtual species for a specific month
  # --- 1.1. Initialize
  set.seed(123)
  ENV_DATA <- list()
  
  # --- 1.2. Get the list of files
  message("--- VIRTUAL : preparing the environmental data")
  list_nc <- list.files(ENV_PATH) %>% 
    .[grep(pattern = ".nc", x = .)]
  
  # --- 1.3. Get the list of variables
  var_out <- ENV_VAR %>% .[grep("!", . , invert = FALSE)] %>% gsub("!", "", .)
  var_in <- ENV_VAR %>% .[grep("!", . , invert = TRUE)]
  
  var_all <- str_sub(list_nc, 1, -4)
  if(length(var_in) != 0){var_all <- var_all[var_all %in% var_in]}
  if(length(var_out) != 0){var_all <- var_all[!c(var_all %in% var_out)]}
  ENV_VAR <- var_all
  
  # --- 2. Extract the data from the .nc files
  stack_month <- lapply(paste0(ENV_PATH, "/", ENV_VAR, ".nc"), 
                        FUN = function(x){x = nc_to_raster(MONTH = MONTH,
                                                           NC = x,
                                                           MIN_DEPTH = 0,
                                                           MAX_DEPTH = 50)
                        }) 
  stack_month <- lapply(stack_month, FUN = function(x)(x = x[[1]])) %>% 
    raster::stack() %>% synchroniseNA()
  names(stack_month) <- ENV_VAR
  
  # --- 3. Generate species from PC
  message("--- VIRTUAL : generate distribution from environmental PCA")
  pca_sp <- list()
  for(i in 1:SP_NB){
    # --- 3.1. Generate species
    tmp <- generateSpFromPCA(raster.stack = stack_month, 
                             niche.breadth = "any", plot = FALSE)
    # --- 3.2. Convert to presence raster
    pca_sp[[i]] <- convertToPA(x = tmp, beta = 0.5)
  }
  
  # --- 4. Generate background samples
  message("--- VIRTUAL : generate bias background from which to sample")
  # --- 4.1. No bias
  background_null <- raster(pca_sp[[1]][["suitab.raster"]]/pca_sp[[1]][["suitab.raster"]])
  
  # --- 4.2. Coast bias
  dist <- background_null
  dist[is.na(dist)] <- 9999
  dist[dist<9999] <- NA
  background_coast <- raster::distance(dist)
  background_coast[background_coast==9999] <- NA
  background_coast <- (background_coast/max(getValues(background_coast), na.rm = TRUE))*(-1) + 1
  background_coast[is.na(background_null)] <- NA
  
  # --- 4.3. North Atlantic bias
  bias_coord <- data.frame(x = -50:0, y = rep(50,51))
  background_natl <- distanceFromPoints(background_null, xy = bias_coord)
  background_natl <- (background_natl/max(getValues(background_natl), na.rm = TRUE))*(-1) + 1
  background_natl[is.na(background_null)] <- NA
  
  # --- 4.4. Stack them to better call them
  bias <- stack(background_null, background_coast, background_natl)
  names(bias) <- c("null", "coast", "natl")
  
  # --- 5. Sample the data
  # --- 5.1. Define parameters
  message(paste("--- VIRTUAL : generate sample for", DATA_TYPE, "data type"))
  sampling_par <- expand.grid(SP = 1:SP_NB, N = N, BIAS = names(bias))
  
  # --- 5.2. Sample binary data
  # Weights have to use terra::rast format... we  **5 to increase the bias a bit
  if(DATA_TYPE == "binary"){
    all_target <- lapply(1:nrow(sampling_par),
                         function(x){
                           tmp <- sampling_par[x,]
                           target <- sampleOccurrences(x = pca_sp[[tmp$SP]], n = tmp$N, plot = TRUE, replacement = TRUE,
                                                       bias = "manual", weights = rast(bias[[tmp$BIAS]]**5)) %>% 
                             .$sample.points %>% 
                             dplyr::select(c(x,y, Real))
                           names(target) <- c("decimallongitude", "decimallatitude", "measurementvalue")
                           
                           target <- target %>% 
                             mutate(worms_id = paste(tmp$SP, tmp$N, names(bias[[tmp$BIAS]]), collapse = "-"),
                                    scientificname = tmp$SP,
                                    depth = 1,
                                    year = "2010",
                                    month = MONTH,
                                    measurementunit = "virtualspecies",
                                    taxonrank = "virtualspecies") %>% 
                             dplyr::select(c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","month","measurementvalue","measurementunit", "taxonrank"))
                           return(target)
                         }) %>% 
      bind_rows()
    
    # Save as CSV file
    write.csv(all_target, file = "/net/meso/work/aschickele/Bluecloud_WB_local/data/virtual_binary.csv", row.names = FALSE)
  } # if binary
  
  # --- 5.3. Sample continuous data
  # Weights have to use terra::rast format... we  **5 to increase the bias a bit
  # --- 5.3.1. Generate samples
  if(DATA_TYPE == "continuous"){
    all_target <- lapply(1:nrow(sampling_par),
                         function(x){
                           tmp <- sampling_par[x,]
                           target <- sampleOccurrences(x = pca_sp[[tmp$SP]], n = tmp$N, plot = TRUE, type = "presence-absence", extract.probability = TRUE, replacement = TRUE,
                                                       bias = "manual", weights = rast(bias[[tmp$BIAS]]**5)) %>% 
                             .$sample.points %>% 
                             dplyr::select(c(x,y, true.probability))
                           names(target) <- c("decimallongitude", "decimallatitude", "measurementvalue")
                           
                           target <- target %>% 
                             mutate(worms_id = paste(tmp$SP, tmp$N, names(bias[[tmp$BIAS]]), collapse = "-"),
                                    scientificname = tmp$SP,
                                    depth = 1,
                                    year = "2010",
                                    month = MONTH,
                                    measurementunit = "virtualspecies",
                                    taxonrank = "virtualspecies") %>% 
                             dplyr::select(c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","month","measurementvalue","measurementunit", "taxonrank"))
                           return(target)
                         }) %>% 
      bind_rows()    
    
    # --- 5.3.2. Add noise
    if(length(NOISE_SD) > 0){
      new_target <- NULL
      new_sampling_par <- NULL
      
      for(i in 1:length(NOISE_SD)){
        noise <- rlnorm(n = nrow(all_target), mean = log(1), sd = log(NOISE_SD[i]))
        tmp <- all_target %>% 
          mutate(measurementvalue = measurementvalue*noise,
                 worms_id = paste(worms_id, NOISE_SD[i]))
        new_target <- bind_rows(new_target, tmp)
        
        tmp <- sampling_par %>% 
          mutate(NOISE_SD = NOISE_SD[i])
        new_sampling_par <- bind_rows(new_sampling_par, tmp)
      } # for i in length noise
    } # end if noise
    
    all_target <- new_target
    sampling_par <- new_sampling_par
    
    # --- 5.3.3. Save as CSV file
    write.csv(all_target, file = "/net/meso/work/aschickele/Bluecloud_WB_local/data/virtual_continuous.csv", row.names = FALSE)
  } # if continuous
  
  # --- 5.3. Sample proportion data
  # TO IMPLEMENT
  
  # --- 6. Wrap up and save
  VIRTUAL <- list()
  VIRTUAL[["par"]] <- list(MONTH = MONTH, SP_NB = SP_NB, N = N, NOISE_SD = NOISE_SD)
  VIRTUAL[["distribution"]] <- pca_sp
  VIRTUAL[["bias"]] <- bias
  VIRTUAL[["sample_parameters"]] <- sampling_par
  
  return(VIRTUAL)
} # END FUNCTION
  



#' =============================================================================
#' @name evaluate_virtual
#' @description evaluate the quality of the pipeline and QC according to the
#' true distributions generated from the virtual species
#' @param FOLDER_NAME the name of the run folder
#' @return a set of plots including the initial virtual species, samples 
#' correlation and taylor diagrams

evaluate_virtual <- function(FOLDER_NAME){
  
  # --- 1. Initialize
  # --- 1.1. Load necessary objects
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/VIRTUAL.RData"))
  SUBFOLDER_NAME <- CALL$SP_SELECT
  
  # --- 1.2. Get virtual projections
  message("--- VIRTUAL : Get virtual projections")
  virtual_proj <- lapply(1:length(VIRTUAL$distribution), function(x){
    out <- VIRTUAL$distribution[[x]]$suitab.raster %>% raster()
    return(out)
  }) %>% stack()
  
  # --- 1.3. Get parameter grid
  message("--- VIRTUAL : Get parameter grid")
  tmp <- VIRTUAL$sample_parameters %>% 
    mutate(SUBFOLDER_NAME = paste(SP, N, BIAS, NOISE_SD))
  virtual_grid <- merge(CALL$MODEL_LIST, tmp)
  colnames(virtual_grid)[1] <- "algorithm"
  
  # --- 1.4. Extract estimated projections
  message("--- VIRTUAL : extract estimated projections")
  estimated_proj <- mclapply(SUBFOLDER_NAME, function(x){
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", x, "/MODEL.RData"))
    
    tmp <- lapply(CALL$MODEL_LIST, function(z){
      out <- MODEL[[z]][["proj"]][["y_hat"]][,,VIRTUAL$par$MONTH] %>%  # take June again
        apply(1, mean)
    })
    tmp <- do.call(rbind, tmp)
    return(tmp)
  }, mc.cores = 20) %>% 
    abind(along = 3) %>% 
    aperm(c(2,1,3))
  dimnames(estimated_proj)[[2]] <- CALL$MODEL_LIST
  dimnames(estimated_proj)[[3]] <- SUBFOLDER_NAME
  
  # --- 2. Build metric data frame
  metric_comparison <- mclapply(1:nrow(virtual_grid), function(x){
    ID <- virtual_grid[x,]
    load(paste0(project_wd, "/output/", FOLDER_NAME,"/", ID$SUBFOLDER_NAME, "/MODEL.RData"))
    EST_FIT <- MODEL[[ID$algorithm]]$eval[[1]] %>% as.numeric()
    CUM_VIP <- MODEL[[ID$algorithm]]$eval$CUM_VIP %>% as.numeric()
    NSD <- MODEL[[ID$algorithm]]$eval$NSD %>% as.numeric()
    true_proj <- getValues(virtual_proj[[ID$SP]])
    est_proj <- estimated_proj[,ID$algorithm,ID$SUBFOLDER_NAME]
    TRUE_FIT <- cor(true_proj, est_proj, method = "pearson", use = "pairwise.complete.obs") %>% 
      as.numeric()
    if(ID$algorithm %in% MODEL$MODEL_LIST){IN_ENSEMBLE <- "Y"} else {IN_ENSEMBLE <- "N"}
    
    out_df <- data.frame(ID,
                         EST_FIT = EST_FIT,
                         CUM_VIP = CUM_VIP,
                         NSD = NSD,
                         TRUE_FIT = TRUE_FIT,
                         IN_ENSEMBLE = IN_ENSEMBLE) 
    return(out_df)
    
  }, mc.cores = 20) %>% do.call(rbind, .) %>% as.data.frame()
  
  # --- 3. Compute Shannon estimate
  # --- 3.1. For the true distributions
  virtual_shannon <- apply(getValues(virtual_proj), 1, function(x)(x = vegan::diversity(x, "shannon")))
  virtual_richness <- apply(getValues(virtual_proj), 1, function(x)(x = sum(x, na.rm = TRUE)))
  
  # --- 3.2. For the estimated distributions
  # --- 3.2.1. Group by ensembles
  estimated_ensemble_grid <- metric_comparison %>% 
    dplyr::select(algorithm, N, BIAS, SUBFOLDER_NAME, IN_ENSEMBLE) %>% 
    dplyr::group_by(SUBFOLDER_NAME) %>% 
    dplyr::filter(IN_ENSEMBLE == "Y") %>% 
    mutate(ID = cur_group_id())
  
  estimated_ensemble <- mclapply(unique(estimated_ensemble_grid$ID), function(x){
    ID <- estimated_ensemble_grid %>% 
      dplyr::filter(ID == x)
    tmp <- estimated_proj[, unique(ID$algorithm), unique(ID$SUBFOLDER_NAME)]
    if(length(dim(tmp)) > 1){tmp <- apply(tmp, 1, mean)}
    return(tmp)
  }, mc.cores = 20) %>% do.call(cbind, .) %>% as.data.frame()
  colnames(estimated_ensemble) <- unique(estimated_ensemble_grid$SUBFOLDER_NAME)
  
  # --- 3.2.2. Estimates across N x BIAS
  estimated_ensemble_grid <- estimated_ensemble_grid %>% 
    dplyr::select(SUBFOLDER_NAME, N, BIAS) %>% 
    dplyr::group_by(N, BIAS) %>% 
    mutate(ID = cur_group_id())
  
  estimated_shannon <- mclapply(unique(estimated_ensemble_grid$ID), function(x){
    ID <- estimated_ensemble_grid %>% 
      dplyr::filter(ID == x)
    tmp <- estimated_ensemble[, unique(ID$SUBFOLDER_NAME)]
    if(length(dim(tmp)) > 1){tmp <- apply(tmp, 1, function(x)(x = vegan::diversity(x, "shannon")))
    } else {tmp <- NULL}
    return(tmp)
  }, mc.cores = 20) %>% do.call(cbind, .) %>% as.data.frame()
  colnames(estimated_shannon) <- unique(estimated_ensemble_grid$ID)
  
  estimated_richness <- mclapply(unique(estimated_ensemble_grid$ID), function(x){
    ID <- estimated_ensemble_grid %>% 
      dplyr::filter(ID == x)
    tmp <- estimated_ensemble[, unique(ID$SUBFOLDER_NAME)]
    if(length(dim(tmp)) > 1){tmp <- apply(tmp, 1, function(x)(x = sum(x, na.rm = TRUE)))
    } else {tmp <- NULL}
    return(tmp)
  }, mc.cores = 20) %>% do.call(cbind, .) %>% as.data.frame()
  colnames(estimated_richness) <- unique(estimated_ensemble_grid$ID)
  
  # --- 4. Do the plots
  # --- 4.1. Initialize
  pdf(paste0(project_wd, "/output/", FOLDER_NAME,"/virtual.pdf"))
  par(xpd = FALSE, mar = c(2,2,2,2))
  
  # --- 4.2. Create land
  land <- virtual_proj[[1]]
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 4.3. Virtual species example
  plot(virtual_proj[[3]], col = inferno_pal(100), legend = FALSE)
  plot(land, add = TRUE, legend = FALSE, col = "antiquewhite4")
  plot(raster::rasterToContour(virtual_proj[[3]], nlevels = 4), add = TRUE)
  
  # --- 4.4. Sample bias
  plot(VIRTUAL$bias$coast, col = viridis_pal(100), legend = FALSE)
  plot(land, add = TRUE, legend = FALSE, col = "antiquewhite4")
  plot(raster::rasterToContour(VIRTUAL$bias$coast, nlevels = 4), add = TRUE)
  
  plot(VIRTUAL$bias$natl, col = viridis_pal(100), legend = FALSE)
  plot(land, add = TRUE, legend = FALSE, col = "antiquewhite4")
  plot(raster::rasterToContour(VIRTUAL$bias$natl, nlevels = 4), add = TRUE)
  
  # --- 4.5. Taylor diagram
  # --- 4.5.1. Initialize - to calibrate Y and X axis with fake distribution
  library(plotrix)
  pal <- alpha(viridis_pal(101), 0.5)
  taylor.diagram(ref = getValues(virtual_proj[[1]]), model = getValues(virtual_proj[[1]])*1.34,
                 normalize = TRUE, col = "black", pcex = 0.01, ref.sd = TRUE)
  # --- 4.5.2. Add all other points
  for(i in 1:nrow(metric_comparison)){
    ID <- metric_comparison[i,]
    REF <- virtual_proj[[ID$SP]] %>% getValues()
    EST <- estimated_proj[, ID$algorithm, ID$SUBFOLDER_NAME]
    PREVALENCE <- VIRTUAL$distribution[[ID$SP]]$PA.conversion["species.prevalence"] %>% 
      as.numeric() %>% round(2) %>% .[]*100
    # COL <- pal[PREVALENCE]
    COL <- pal[ID$NOISE_SD]
    
    if(ID$IN_ENSEMBLE == "Y"){
      taylor.diagram(ref = REF, model = EST, normalize = TRUE, col = COL, add = TRUE)
      taylor.diagram(ref = REF, model = EST, normalize = TRUE, pch = 1, col = alpha("black", 0.5), add = TRUE)
    } else {
      taylor.diagram(ref = REF, model = EST, normalize = TRUE, pch = 1, col = COL, add = TRUE)
    }
  } # end i loop
  
  # --- 4.5.3. Add richness metrics
  for(i in 1:ncol(estimated_shannon)){
    taylor.diagram(ref = virtual_shannon, model = estimated_shannon[,i], normalize = TRUE, col = alpha("orange", 0.5), pch = 18, pcex = 2, add = TRUE)
    taylor.diagram(ref = virtual_shannon, model = estimated_shannon[,i], normalize = TRUE, col = alpha("black", 0.5), pch = 5, pcex = 2, add = TRUE)
  }
  for(i in 1:ncol(estimated_richness)){
    taylor.diagram(ref = virtual_richness, model = estimated_richness[,i], normalize = TRUE, col = alpha("red", 0.5), pch = 18, pcex = 2, add = TRUE)
    taylor.diagram(ref = virtual_richness, model = estimated_richness[,i], normalize = TRUE, col = alpha("black", 0.5), pch = 5, pcex = 2, add = TRUE)
  }
  
  # --- 4.6. Estimated QC vs True QC
  # --- 4.6.1. Extract corresponding metric - for future decision making
  decision <- metric_comparison %>% 
    dplyr::filter(IN_ENSEMBLE == "Y") %>% 
    dplyr::select(EST_FIT, TRUE_FIT)
  
  # --- 4.6.1. Initialize - to calibrate Y and X axis with fake distribution
  par(xpd = FALSE)
  plot(x = 1, y = 1, type = 'n', xlim = c(0, 1), ylim = c(0,1), xlab = "Pearson R", ylab = "Fit estimate",
       main = paste("R = ", round(cor(decision$EST_FIT, decision$TRUE_FIT),2), "-",
                    "RMSE =", round(cor(decision$EST_FIT, decision$TRUE_FIT),2), "-",
                    "Median =", round(median(decision$EST_FIT), 2)
       ))
  abline(coef = c(0,1))
  abline(h = 0.5, v = 0.5, col = "darkred")
  grid(col = "gray50")
  # --- 4.6.2. Add all other points
  for(i in 1:nrow(metric_comparison)){
    ID <- metric_comparison[i,]
    PREVALENCE <- VIRTUAL$distribution[[ID$SP]]$PA.conversion["species.prevalence"] %>% 
      as.numeric() %>% round(2) %>% .[]*100
    # COL <- pal[PREVALENCE]
    COL <- pal[ID$NOISE_SD]
    
    if(ID$IN_ENSEMBLE == "Y"){
      points(x = ID$TRUE_FIT, y = ID$EST_FIT, bg = COL, pch = 21)
    } else {
      points(x = ID$TRUE_FIT, y = ID$EST_FIT, col = COL, pch = 1)
    }
  } # end i loop
  
  # --- 4.7. Wrap up plots
  dev.off()

} # End function
# --- END --- 