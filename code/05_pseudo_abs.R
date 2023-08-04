#' =============================================================================
#' @name pseudo_abs
#' @description computes pseudo absences added to the X and Y matrices from
#' query_env and bio according to various background selection methods.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param METHOD_PA method of pseudo-absence, either "mindist" or "cumdist" or "density"
#' @param NB_PA number of pseudo-absences to generate
#' @param DIST_PA if METHOD_PA = "mindist", distance from presences (in meters),
#'  from which to define the background data. Expert use only.
#' @param BACKGROUND_FILTER additional background filter for finer tuning, such
#' as selecting pseudo-absences within the sampled background of a given campaign
#' or instrument deployment. Passed by the user in the form of a 2 column 
#' data frame, x = longitude and y = latitude where the pseudo-absences
#' can be sampled. Or a path to a raster object where pseudo-absences are sampled in
#' non NA cells, weighted by the cell values.
#' @return X updated with the pseudo-absence values (= 0)
#' @return Y updated with the environmental values corresponding
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

pseudo_abs <- function(FOLDER_NAME = NULL,
                       SUBFOLDER_NAME = NULL,
                       METHOD_PA = "cumdist",
                       NB_PA = NULL,
                       DIST_PA = 1000e3,
                       PER_RANDOM = 0.25,
                       BACKGROUND_FILTER = NULL){
  
  # --- 1. Initialize function
  set.seed(123)
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : pseudo_absences ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  if(is.null(NB_PA)){NB_PA = nrow(QUERY$S)}
  
  # --- 1.3. Double check data type - plot observation locations anyway
  if(CALL$DATA_TYPE != "binary"){
    message("No Pseudo-absence generation necessary for this data type")
    r <- raster(CALL$ENV_PATH)
    r[!is.na(r)] <- 0
    
    land <- r
    land[is.na(land)] <- 9999
    land[land != 9999] <- NA
    
    pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_observations.pdf"))
    plot(land, col = "antiquewhite4", legend=FALSE, main = paste("Observation locations for:", SUBFOLDER_NAME), 
         sub = paste("NB_OBS :",  nrow(QUERY$S)))
    points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, col = "black", pch = 3)
    dev.off()
    
    log_sink(FILE = sinkfile, START = FALSE)
    return(NULL)
  } 

  # --- 2. Load features and create base raster
  features <- stack(CALL$ENV_PATH) %>% readAll()
  r <- features[[1]]
  r[!is.na(r)] <- 0
  
  # --- 3. Background definition
  # --- 3.1. Based on a geographical distance
  # Defining the background as cells distant from more than n-km from presence
  if(METHOD_PA == "mindist"){
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
      dplyr::filter(layer < DIST_PA & !is.na(layer)) %>% 
      dplyr::select(x, y)
  } # End if geo
  
  # --- 3.2. Random but biased by cumulative-distance to presence
  if(METHOD_PA == "cumdist"){
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
  if(METHOD_PA == "density"){
    # --- 3.3.1. Extract presence points
    presence <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude) %>% 
      rasterize(r)
    presence[!is.na(presence)] <- 1
    
    # --- 3.3.2. Compute density
    focal_w <- focalWeight(r, 20, type = "Gauss")
    dens <- focal(presence, focal_w, fun = function(x){sum(x, na.rm = TRUE)}, pad = TRUE)
    dens <- (dens/max(getValues(dens), na.rm = TRUE)*(1-PER_RANDOM))+PER_RANDOM # try to get a fix random PA generation
    dens[!is.na(presence)] <- NA
    
    # --- 3.3.3. Define the weighted background raster
    background <- synchroniseNA(stack(dens, r))[[1]] %>% 
      rasterToPoints() %>% 
      as.data.frame()
  } # End if density
  
  # --- 4. Additional background filter
  if(!is.null(BACKGROUND_FILTER)){
    if(is.data.frame(BACKGROUND_FILTER)){
      background <- BACKGROUND_FILTER
    } else {
      background <- raster(BACKGROUND_FILTER) %>% 
        rasterToPoints() %>% 
        as.data.frame()
    }
  } # END if background_filter TRUE

  # --- 5. Sample within the background data
  # --- 5.1. Conditional sampling
  # Add a resample option if there is not enough background available
  if(ncol(background == 3)){
    if(nrow(background) < NB_PA){
      message(" PSEUDO-ABS : background too small, selection with replacement !")
      tmp <- sample(x = 1:nrow(background), size = NB_PA, replace = TRUE, prob = background[,3]**2) # sqr the distance to bias more
    } else {
      tmp <- sample(x = 1:nrow(background), size = NB_PA, prob = background[,3]**2) # sqr the distance to bias more
    } 
  } else {
    if(nrow(background) < NB_PA){
      message(" PSEUDO-ABS : background too small, selection with replacement !")
      tmp <- sample(x = 1:nrow(background), size = NB_PA, replace = TRUE)
    } else {
      tmp <- sample(x = 1:nrow(background), size = NB_PA)
    } 
  }

  # --- 5.2. Subset the background coordinates
  xy <- background[tmp,1:2]
  
  # --- 5.3. Fast PDF to check the absences location
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/01_pseudo_abs.pdf"))
  land <- r
  land[is.na(land)] <- 9999
  land[land != 9999] <- NA
  plot(land, col = "antiquewhite4", legend=FALSE, main = paste("Presence - Pseudo Abs for:", SUBFOLDER_NAME), 
       sub = paste("NB_OBS :",  nrow(QUERY$S),"// NB_PA :", NB_PA, "// METHOD_PA :", METHOD_PA))
  points(xy, col = "red", pch = 3)
  points(QUERY$S$decimallongitude, QUERY$S$decimallatitude, col = "black", pch = 3)
  dev.off()
  
  # --- 6. Append the query
  # --- 6.1. Feature table
  X <- raster::extract(features, xy) %>% 
    as.data.frame()
  QUERY$X <- rbind(QUERY$X, X)
  
  # --- 6.2. Target table - replace by 0 and 1's
  QUERY$Y <- data.frame(measurementvalue = c(rep(1, nrow(QUERY$Y)), 
                                             rep(0, nrow(xy))))
  
  # --- 6.3. Sample table
  S <- data.frame(decimallongitude = xy$x,
                  decimallatitude = xy$y,
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
