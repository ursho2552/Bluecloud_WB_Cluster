#' =============================================================================
#' @name pseudo_abs
#' @description computes pseudo absences added to the X and Y matrices from
#' query_env and bio according to various background selection methods.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param METHOD_PA method of pseudo-absence, either "env" or "geo"
#' @param NB_PA number of pseudo-absences to generate
#' @param DIST_PA if METHOD_PA = "geo", distance from presences (in meters),
#'  from which to define the background data. Expert use only.
#' @param BACKGROUND_FILTER additional background filter for finer tuning, such
#' as selecting pseudo-absences within the sampled background of a given campaign
#' or instrument deployment. Passed by the user in the form of a 2 column 
#' data frame, x = longitude and y = latitude where the pseudo-absences
#' can be sampled. Expert use only. 
#' @return X updated with the pseudo-absence values (= 0)
#' @return Y updated with the environmental values corresponding
#' @return Updates the output in a QUERY.RData and CALL.Rdata files

pseudo_abs <- function(FOLDER_NAME = NULL,
                       SUBFOLDER_NAME = NULL,
                       METHOD_PA = "env",
                       NB_PA = NULL,
                       DIST_PA = 100e3,
                       BACKGROUND_FILTER = NULL){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : pseudo_absences ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  if(is.null(NB_PA)){NB_PA = nrow(QUERY$S)}
  
  # --- 1.3. Double check data type
  if(CALL$DATA_TYPE != "pres"){
    stop("No Pseudo-absence generation necessary for this data type")
  } 
  
  # --- 2. Open environmental data
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>%
    readAll()
  
  # --- 3. Create base raster
  r <- features[[1]]
  r[!is.na(r)] <- 0
  
  # --- 4. Environmental background computation
  # --- 4.1. Based on an environmental envelope
  # Here, we are using the MESS analysis to approximate environmental background
  # outside the environmental space of presence
  if(METHOD_PA == "env"){
    if(is.null(QUERY$MESS)){
      stop("env Pseudo-Absence generation requires a MESS analysis in the previous step \n")
    } else {
      background <- QUERY$MESS
      background <- synchroniseNA(stack(background, r))[[1]] %>% 
        rasterToPoints() %>% 
        as.data.frame() %>% 
        dplyr::filter(mess < 0 & !is.na(mess)) %>% 
        dplyr::select(x, y)
    }
  } # End if env

  # --- 4.2. Based on a geographical distance
  # Defining the background as cells distant from more than n-km from presence
  if(METHOD_PA == "geo"){
    # --- 4.2.1. Extract presence points
    val <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude)
    
    # --- 4.2.2. Calculate distance to presences in raster object
    background <- rasterize(val, r, update=TRUE)
    background[background < 1] <- NA
    background <- distance(background)
    
    background <- synchroniseNA(stack(background, r))[[1]] %>% 
      rasterToPoints() %>% 
      as.data.frame() %>% 
      dplyr::filter(layer > DIST_PA & !is.na(layer)) %>% 
      dplyr::select(x, y)
  } # End if geo
  
  # --- 5. Additional background filter
  # TO UPDATE LATER

  # --- 6. Sample within the background data
  # Add a resample option if there is not enough background available
  if(nrow(background) < NB_PA){
    message(" PSEUDO-ABS : background too small, selection with replacement !")
    tmp <- sample(x = 1:nrow(background), size = NB_PA, replace = TRUE)
  } else {
    tmp <- sample(x = 1:nrow(background), size = NB_PA)
  }
  # Subset the background coordinates
  xy <- background[tmp,]
  
  # --- 7. Append the query
  # --- 7.1. Feature table
  X <- raster::extract(features, xy) %>% 
    as.data.frame()
  QUERY$X <- rbind(QUERY$X, X)
  
  # --- 7.2. Target table - replace by 0 and 1's
  QUERY$Y <- data.frame(measurementvalue = c(rep(1, nrow(QUERY$Y)), 
                                             rep(0, nrow(xy))))
  
  # --- 7.3. Sample table
  S <- data.frame(decimallongitude = xy$x,
                  decimallatitude = xy$y,
                  measurementtype = "Pseudo-absence") %>% 
    mutate(ID = row_number()+nrow(QUERY$S))
  QUERY$S <- QUERY$S %>% 
    bind_rows(S)

  # --- 8. Wrap up and save
  # --- 8.1. Save QUERY object
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  # --- 8.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
