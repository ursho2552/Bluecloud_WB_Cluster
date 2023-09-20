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
    
    hist(QUERY$S$decimallatitude, breaks = seq(-90,90, 10), col = alpha("black", 0.5), 
         main = "Latitudinal spectrum", xlab = "Latitude")
    
    
    dev.off()
    
    log_sink(FILE = sinkfile, START = FALSE)
    return(NULL)
  } # if datatype !binary
  
  # --- 3. Initialize background definition
  # --- 3.1. Extract the monthly bias
  # To know how many pseudo-absences to draw from each month background
  month_freq <- summary(as.factor(QUERY$S$month))
  month_freq <- trunc(month_freq/(sum(month_freq)/CALL$NB_PA))
  
  # --- 3.2. Initialize lon x lat x time df
  xym <- NULL
  
  # --- 4. Extract pseudo-absence
  for(m in names(month_freq)){
    # --- 4.1. Extract presence points
    presence <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude, month) %>% 
      dplyr::filter(month == m) %>% 
      dplyr::select(decimallongitude, decimallatitude)
    
    # --- 4.2. Background definition
    # Need to redefine base raster here in case of different NA between months
    r <- CALL$ENV_DATA[[as.numeric(m)]][[1]]
    r[!is.na(r)] <- 0
    
    # --- 4.2.1. Based on a geographical distance
    # Defining the background as cells distant from more than n-km from presence
    if(CALL$METHOD_PA == "mindist"){
      # Calculate distance to presences in raster object
      background <- raster::rasterize(presence, r)
      background[background < 1] <- NA
      background <- raster::distance(background)
      
      # Filter by minimum distance to presence
      background <- synchroniseNA(stack(background, r))[[1]] %>% 
        rasterToPoints() %>% 
        as.data.frame() %>% 
        dplyr::filter(layer < CALL$DIST_PA & !is.na(layer)) %>% 
        dplyr::select(x, y)
    } # End if mindist
    
    # --- 4.2.2. Biased by cumulative-distance to presence
    if(CALL$METHOD_PA == "cumdist"){
      
      # Compute cumulative distance
      background <- r %>% rasterToPoints() %>% .[,1:2]
      cumdist <- pointDistance(presence, background, lonlat = TRUE) %>% 
        apply(2, sum)
      
      # Define weighted background raster
      background <- r
      background[!is.na(background)] <- cumdist
      background <- (background/max(getValues(background), na.rm = TRUE)-1)*-1
      
      background <- synchroniseNA(stack(background, r))[[1]] %>% 
        rasterToPoints() %>% 
        as.data.frame()
    } # End if cumdist
    
    # --- 4.2.3. Density of presence within a buffer
    if(CALL$METHOD_PA == "density"){
      # Prepare base raster
      presence <- raster::rasterize(presence, r)
      presence[!is.na(presence)] <- 1
      
      # Compute density
      focal_w <- focalWeight(r, 20, type = "Gauss")
      dens <- focal(presence, focal_w, fun = function(x){sum(x, na.rm = TRUE)}, pad = TRUE)
      dens <- (dens/max(getValues(dens), na.rm = TRUE)*(1-CALL$PER_RANDOM))+CALL$PER_RANDOM # try to get a fix random PA generation
      dens[!is.na(presence)] <- NA
      
      # Define the weighted background raster
      background <- synchroniseNA(stack(dens, r))[[1]] %>% 
        rasterToPoints() %>% 
        as.data.frame()
    } # End if density
    
    # --- 4.2.4. Additional background filter
    if(!is.null(CALL$BACKGROUND_FILTER)){
      if(is.data.frame(CALL$BACKGROUND_FILTER)){
        background <- CALL$BACKGROUND_FILTER
      } else {
        background <- raster(CALL$BACKGROUND_FILTER) %>% 
          rasterToPoints() %>% 
          as.data.frame()
      }
    } # END if background_filter TRUE
    
    # --- 4.3. Sample within the background data
    # Add a resample option if there is not enough background available
    if(ncol(background == 3)){
      if(nrow(background) < month_freq[m]){
        message(" PSEUDO-ABS : background too small, selection with replacement !")
        tmp <- sample(x = 1:nrow(background), size = month_freq[m], replace = TRUE, prob = background[,3]) # sqr the distance to bias more
      } else {
        tmp <- sample(x = 1:nrow(background), size = month_freq[m], prob = background[,3]) # sqr the distance to bias more
      } 
    } else {
      if(nrow(background) < month_freq[m]){
        message(" PSEUDO-ABS : background too small, selection with replacement !")
        tmp <- sample(x = 1:nrow(background), size = month_freq[m], replace = TRUE)
      } else {
        tmp <- sample(x = 1:nrow(background), size = month_freq[m])
      } 
    }
    
    # --- 4.4. Concatenate over month
    xym0 <- background[tmp,1:2] %>% cbind(rep(as.numeric(m), month_freq[m]))
    xym <- rbind(xym, xym0)
    
  } # m month loop
  
  # --- 5.Fast PDF to check the absences location
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_pseudo_abs.pdf"))
  par(mfrow = c(2,2))
  
  # --- 5.1. Geographical distribution
  plot(land, col = "antiquewhite4", legend=FALSE, main = paste("Obs. presence for:", SUBFOLDER_NAME), 
       sub = paste("NB_OBS :",  nrow(QUERY$S)))
  points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, col = alpha("black", 0.2), pch = 20)
  
  plot(land, col = "antiquewhite4", legend=FALSE, main = paste("Pseudo Abs for:", SUBFOLDER_NAME), 
       sub = paste("NB_PA :", CALL$NB_PA, "// METHOD_PA :", CALL$METHOD_PA))
  points(xym[,1:2], col = alpha("red", 0.2), pch = 20)
  
  # --- 5.2. Longitude and latitude profile
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
                  month = xym[,3],
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
