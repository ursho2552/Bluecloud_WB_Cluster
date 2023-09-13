#' =============================================================================
#' @name pseudo_abs
#' @description computes pseudo absences added to the X and Y matrices from
#' query_env and bio according to various background selection methods.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return X updated with the pseudo-absence values (= 0)
#' @return Y updated with the environmental values corresponding
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

pseudo_abs <- function(FOLDER_NAME = NULL,
                       SUBFOLDER_NAME = NULL){
  
  # --- 1. Initialize function
  set.seed(123)
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : pseudo_absences ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  if(is.null(CALL$NB_PA)){CALL$NB_PA = nrow(QUERY$S)}
  
  # --- 1.3. Base raster and land
  r <- CALL$ENV_DATA[[1]][[1]]
  r[!is.na(r)] <- 0
  
  land <- r
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  
  # --- 2. Double check data type - plot observation locations anyway
  if(CALL$DATA_TYPE != "binary"){
    message("No Pseudo-absence generation necessary for this data type")
    pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_observations.pdf"))
    par(mfrow = c(2,2))
    
    # --- 2.1.1. Geographical plot for continuous data
    # We plot an artificial land with the scale of the observation first, for the legend
    if(CALL$DATA_TYPE == "continuous"){
      tmp <- (land-9998)*max(QUERY$Y$measurementvalue)
      tmp[1] <- 0
      plot(tmp, col = inferno_pal(100), main = paste("Observation locations for:", SUBFOLDER_NAME), 
           sub = paste("NB_OBS :",  nrow(QUERY$S)))
      plot(land, col = "antiquewhite4", legend=FALSE, add = TRUE)
      points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, 
             col = col_numeric("inferno", domain = range(QUERY$Y$measurementvalue, na.rm = TRUE))(QUERY$Y$measurementvalue),
             pch = 20)
    }

    # --- 2.1.2. Geographical plot for proportion type
    if(CALL$DATA_TYPE == "proportions"){
      plot(land, col = "antiquewhite4", legend=FALSE, main = paste("Observation locations for:", SUBFOLDER_NAME), 
           sub = paste("NB_OBS :",  nrow(QUERY$S)))
      points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, 
             col = "black", pch = 20)
    }
    
    # --- 2.2. Histogram of values
    hist(unlist(QUERY$Y), breaks = 25, col = alpha("black", 0.5), main = "Histogram of values", xlab = CALL$DATA_TYPE)
    
    # --- 2.3. Longitude and latitude profile
    hist(QUERY$S$decimallongitude, breaks = seq(-180,180, 10), col = alpha("black", 0.5), 
         main = "Longitudinal spectrum", xlab = "Longitude")
    
    hist(QUERY$S$decimallatitude, breaks = seq(-180,180, 10), col = alpha("black", 0.5), 
         main = "Latitudinal spectrum", xlab = "Latitude")
    
    
    dev.off()
    
    log_sink(FILE = sinkfile, START = FALSE)
    return(NULL)
  } # if datatype !binary
  
  # --- 3. Background definition
  # --- 3.1. Based on a geographical distance
  # Defining the background as cells distant from more than n-km from presence
  if(CALL$METHOD_PA == "mindist"){
    # --- 3.1.1. Extract presence points
    val <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude)
    
    # --- 3.1.2. Calculate distance to presences in raster object
    background <- rasterize(val, r, update=TRUE)
    background[background < 1] <- NA
    background <- raster::distance(background)
    
    background <- synchroniseNA(stack(background, r))[[1]] %>% 
      rasterToPoints() %>% 
      as.data.frame() %>% 
      dplyr::filter(layer < CALL$DIST_PA & !is.na(layer)) %>% 
      dplyr::select(x, y)
  } # End if geo
  
  # --- 3.2. Random but biased by cumulative-distance to presence
  if(CALL$METHOD_PA == "cumdist"){
    # --- 3.2.1. Extract presence points
    presence <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude)
    
    # --- 3.2.2. Compute cumulative distance
    background <- r %>% rasterToPoints() %>% .[,1:2]
    cumdist <- pointDistance(presence, background, lonlat = TRUE) %>% 
      apply(2, sum)
    
    # --- 3.2.3. Define weighted background raster
    background <- r
    background[!is.na(background)] <- cumdist
    background <- (background/max(getValues(background), na.rm = TRUE)-1)*-1
    
    background <- synchroniseNA(stack(background, r))[[1]] %>% 
      rasterToPoints() %>% 
      as.data.frame()
  } # End if bias_env
  
  # --- 3.3. Density of presence within a buffer
  if(CALL$METHOD_PA == "density"){
    # --- 3.3.1. Extract presence points
    presence <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude) %>% 
      rasterize(r)
    presence[!is.na(presence)] <- 1
    
    # --- 3.3.2. Compute density
    focal_w <- focalWeight(r, 20, type = "Gauss")
    dens <- focal(presence, focal_w, fun = function(x){sum(x, na.rm = TRUE)}, pad = TRUE)
    dens <- (dens/max(getValues(dens), na.rm = TRUE)*(1-CALL$PER_RANDOM))+CALL$PER_RANDOM # try to get a fix random PA generation
    dens[!is.na(presence)] <- NA
    
    # --- 3.3.3. Define the weighted background raster
    background <- synchroniseNA(stack(dens, r))[[1]] %>% 
      rasterToPoints() %>% 
      as.data.frame()
  } # End if density
  
  # --- 4. Additional background filter
  if(!is.null(CALL$BACKGROUND_FILTER)){
    if(is.data.frame(CALL$BACKGROUND_FILTER)){
      background <- CALL$BACKGROUND_FILTER
    } else {
      background <- raster(CALL$BACKGROUND_FILTER) %>% 
        rasterToPoints() %>% 
        as.data.frame()
    }
  } # END if background_filter TRUE

  # --- 5. Sample within the background data
  # --- 5.1. Conditional sampling locations
  # Add a resample option if there is not enough background available
  if(ncol(background == 3)){
    if(nrow(background) < CALL$NB_PA){
      message(" PSEUDO-ABS : background too small, selection with replacement !")
      tmp <- sample(x = 1:nrow(background), size = CALL$NB_PA, replace = TRUE, prob = background[,3]**2) # sqr the distance to bias more
    } else {
      tmp <- sample(x = 1:nrow(background), size = CALL$NB_PA, prob = background[,3]**2) # sqr the distance to bias more
    } 
  } else {
    if(nrow(background) < CALL$NB_PA){
      message(" PSEUDO-ABS : background too small, selection with replacement !")
      tmp <- sample(x = 1:nrow(background), size = CALL$NB_PA, replace = TRUE)
    } else {
      tmp <- sample(x = 1:nrow(background), size = CALL$NB_PA)
    } 
  }

  # --- 5.2. Sampling month
  # According to the same distribution as the one in the data
  month_freq <- summary(as.factor(QUERY$S$month))
  month_freq <- month_freq/sum(month_freq)
  month_prob <- sample(x = 1:length(month_freq), size = CALL$NB_PA, replace = TRUE, prob = month_freq)
  
  # --- 5.2. Subset the background coordinates
  xym <- background[tmp,1:2] %>% cbind(month_prob)
  
  # --- 5.4. Fast PDF to check the absences location
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_pseudo_abs.pdf"))
  par(mfrow = c(2,2))
  
  # --- 5.4.1. Geographical distribution
  plot(land, col = "antiquewhite4", legend=FALSE, main = paste("Obs. presence for:", SUBFOLDER_NAME), 
       sub = paste("NB_OBS :",  nrow(QUERY$S)))
  points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, col = alpha("black", 0.2), pch = 20)
  
  plot(land, col = "antiquewhite4", legend=FALSE, main = paste("Pseudo Abs for:", SUBFOLDER_NAME), 
       sub = paste("NB_PA :", CALL$NB_PA, "// METHOD_PA :", CALL$METHOD_PA))
  points(xym[,1:2], col = alpha("red", 0.2), pch = 20)
  
  # --- 5.4.2. Longitude and latitude profile
  hist(QUERY$S$decimallongitude, breaks = seq(-180,180, 10), col = alpha("black", 0.5), 
       main = "Longitudinal spectrum", xlab = "Longitude")
  hist(xym$x, breaks = seq(-180,180, 10), col = alpha("red", 0.5), add = TRUE)
  
  hist(QUERY$S$decimallatitude, breaks = seq(-180,180, 10), col = alpha("black", 0.5), 
       main = "Latitudinal spectrum", xlab = "Latitude")
  hist(xym$y, breaks = seq(-180,180, 10), col = alpha("red", 0.5), add = TRUE)
  dev.off()
  
  # --- 6. Append the query
  # --- 6.1. Feature table
  # Extract from the corresponding monthly raster
  X <- NULL
  for(i in 1:nrow(xym)){
    tmp <- raster::extract(CALL$ENV_DATA[[xym[i,3]]], xym[i,1:2]) %>% 
      as.data.frame()
    X <- rbind(X, tmp)
  }
  QUERY$X <- rbind(QUERY$X, X)
  
  # --- 6.2. Target table - replace by 0 and 1's
  QUERY$Y <- data.frame(measurementvalue = c(rep(1, nrow(QUERY$Y)), 
                                             rep(0, nrow(xym))))
  
  # --- 6.3. Sample table
  S <- data.frame(decimallongitude = xym$x,
                  decimallatitude = xym$y,
                  month = xym$month_prob,
                  measurementtype = "Pseudo-absence") %>% 
    mutate(ID = row_number()+nrow(QUERY$S))
  QUERY$S <- QUERY$S %>% 
    bind_rows(S)

  # --- 7. Wrap up and save
  # --- 7.1. Save QUERY object
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
