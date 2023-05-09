#' =============================================================================
#' @name pseudo_abs
#' @description computes pseudo absences added to the X and Y matrices from
#' query_env and bio according to various background selection methods.
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
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

pseudo_abs <- function(SP_SELECT = NULL,
                       FOLDER_NAME = NULL,
                       METHOD_PA = "env",
                       NB_PA = NULL,
                       DIST_PA = 100e3,
                       BACKGROUND_FILTER = NULL){
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  if(is.null(NB_PA)){NB_PA = nrow(QUERY$S)}
  
  # =========================== DATA TYPE CHECK ================================
  if(CALL$DATA_TYPE != "pres"){
    stop("No Pseudo-absence generation necessary for this data type")
  } 
  
  # --- 1. Open environmental data
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>%
    readAll()
  
  # --- 2. Create base raster
  r <- features[[1]]
  r[!is.na(r)] <- 0
  
  # --- 3. Environmental background computation
  # --- 3.1. Based on an environmental envelope
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

  # --- 3.2. Based on a geographical distance
  # Defining the background as cells distant from more than n-km from presence
  if(METHOD_PA == "geo"){
    # --- 3.2.1. Extract presence points
    val <- QUERY$S %>% 
      dplyr::select(decimallongitude, decimallatitude)
    
    # --- 3.2.2. Calculate distance to presences in raster object
    background <- rasterize(val, r, update=TRUE)
    background[background < 1] <- NA
    background <- distance(background)
    
    background <- synchroniseNA(stack(background, r))[[1]] %>% 
      rasterToPoints() %>% 
      as.data.frame() %>% 
      dplyr::filter(layer > DIST_PA & !is.na(layer)) %>% 
      dplyr::select(x, y)
  } # End if geo
  
  # --- 4. Additional background filter
  # TO UPDATE LATER

  # --- 5. Sample within the background data
  # Add a resample option if there is not enough background available
  if(nrow(background) < NB_PA){
    message(" PSEUDO-ABS : background too small, selection with replacement !")
    tmp <- sample(x = 1:nrow(background), size = NB_PA, replace = TRUE)
  } else {
    tmp <- sample(x = 1:nrow(background), size = NB_PA)
  }
  # Subset the background coordinates
  xy <- background[tmp,]
  
  # --- 6. Append the query
  # --- 6a. Feature table
  X <- raster::extract(features, xy) %>% 
    as.data.frame()
  QUERY$X <- rbind(QUERY$X, X)
  
  # --- 6b. Target table - replace by 0 and 1's
  QUERY$Y <- data.frame(measurementvalue = c(rep(1, nrow(QUERY$Y)), 
                                             rep(0, nrow(xy))))
  
  # --- 6c. Sample table
  S <- data.frame(decimallongitude = xy$x,
                  decimallatitude = xy$y,
                  measurementtype = "Pseudo-absence") %>% 
    mutate(ID = row_number()+nrow(QUERY$S))
  QUERY$S <- QUERY$S %>% 
    bind_rows(S)

  # --- 7. Save QUERY object
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  
} # END FUNCTION
