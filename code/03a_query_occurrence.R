#' =============================================================================
#' @name query_occurrence
#' @description extracts biological from the OBIS database, according to a
#' user provided list of species, time and depth range. The extracted data is
#' formatted to be directly usable by the models available in this workbench.
#' @param FOLDER_NAME name of the corresponding folder
#' @param QUERY the QUERY.RData object containing the list of species
#' (depth range later - TO IMPLEMENT)
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth
#' @return Annotations: a data frame of taxonomic (maybe functional) annotation
#' for each target species.
#' @return the output in a QUERY object

query_occurrence <- function(FOLDER_NAME = NULL,
                             QUERY = NULL){
  
  # --- 1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 2. Query OBIS occurrences
  # --- 2.1. Default univariate target
  target <- occurrence(taxonid = QUERY$SUBFOLDER_INFO$SP_SELECT) %>% 
    dplyr::filter(aphiaID == QUERY$SUBFOLDER_INFO$SP_SELECT) %>% # triple check !
    dplyr::filter(basisOfRecord == "Occurrence") %>% 
    dplyr::filter(decimalLatitude != 0 | decimalLatitude < -90 | decimalLatitude > 90) %>% 
    dplyr::filter(decimalLongitude != 0 | decimalLongitude < -180 | decimalLongitude > 180) %>% 
    dplyr::filter(occurrenceStatus == "present") %>% 
    dplyr::select(any_of(c("scientificName", "aphiaID", "decimalLatitude", "decimalLongitude", "depth", "date_year", "occurrenceStatus", "basisOfRecord", "taxonRank"))) %>% 
    dplyr::filter(date_year >= CALL$SAMPLE_SELECT$START_YEAR & date_year <= CALL$SAMPLE_SELECT$STOP_YEAR) %>% 
    dplyr::filter(depth >= CALL$SAMPLE_SELECT$MIN_DEPTH & depth <= CALL$SAMPLE_SELECT$MAX_DEPTH) %>% 
    distinct()
  
  # --- 2.2. Column repair
  if("taxonRank" %in% names(target) == FALSE){
    target <- mutate(target, taxonrank = NA)
  }
  
  # --- 2.3. Add number of occurrence
  target <- mutate(target, nb_occ = n())
  
  # --- 2.3. Nice names
  colnames(target) <- c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","measurementvalue","measurementunit", "taxonrank","nb_occ")
  
  # --- 3. Create Y target table
  Y <- target %>% 
    dplyr::select(measurementvalue)
  
  # --- 4. Create S sample table
  # Add the row ID in S to keep track of rows in sub-sample, train, test etc...
  S <- target %>% 
    dplyr::select(-any_of(c("measurementvalue", "worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           ID = row_number())
  
  # --- 5. Create an Annotation table
  annotations <- target %>% 
    dplyr::select(any_of(c("worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    distinct()
  
  # --- 6. Save in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)
  
} # END FUNCTION
