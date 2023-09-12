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
    dplyr::filter(basisOfRecord == "Occurrence" | basisOfRecord == "HumanObservation" | basisOfRecord == "LivingSpecimen") %>% 
    dplyr::filter(decimalLatitude != 0 | decimalLatitude < -90 | decimalLatitude > 90) %>% 
    dplyr::filter(decimalLongitude != 0 | decimalLongitude < -180 | decimalLongitude > 180) %>% 
    dplyr::filter(occurrenceStatus == "present") %>% 
    dplyr::select(any_of(c("scientificName", "aphiaID", "decimalLatitude", "decimalLongitude", "depth", "date_year", "month", "occurrenceStatus", "basisOfRecord", "taxonRank"))) %>% 
    dplyr::filter(date_year >= CALL$SAMPLE_SELECT$START_YEAR & date_year <= CALL$SAMPLE_SELECT$STOP_YEAR) %>% 
    dplyr::filter(depth >= CALL$SAMPLE_SELECT$MIN_DEPTH & depth <= CALL$SAMPLE_SELECT$MAX_DEPTH) %>% 
    distinct()
  
  # --- 2.2. Column repair
  if("taxonRank" %in% names(target) == FALSE){
    target <- mutate(target, taxonrank = "undefined")
  }
  
  # --- 2.3. Nice names
  colnames(target) <- c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year", "month","measurementvalue","measurementunit", "taxonrank")
  
  # --- 3. Query GBIF occurrences
  # --- 3.1. Prepare the scientific name - Aphia ID does not work for GBIF
  SNAME <- target$scientificname[1] %>% as.character()
  
  # --- 3.2. Default univariate target
  target_gbif <- lapply(CALL$SAMPLE_SELECT$START_YEAR:CALL$SAMPLE_SELECT$STOP_YEAR,
                 FUN = function(YEAR){
                   occ_data(scientificName = SNAME,
                            year = YEAR,
                            depth = paste0(CALL$SAMPLE_SELECT$MIN_DEPTH, ",", CALL$SAMPLE_SELECT$MAX_DEPTH),
                            occurrenceStatus = 'PRESENT',
                            limit = 99000)$data
                 }) %>% 
    bind_rows() 
  
  # --- 3.3. Re-format data frame
  # Security if there is not GBIF data available for the species
  if(nrow(target_gbif) != 0){
    # --- 3.3.1. Select columns
    # And re-select the scientificname as a double check
    target_gbif <- target_gbif %>% 
      dplyr::filter(grepl(SNAME, scientificName)) %>% 
      dplyr::select(any_of(c("decimalLatitude","decimalLongitude","depth","year","month"))) %>% 
      dplyr::filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))
    colnames(target_gbif) <- c("decimallatitude","decimallongitude","depth","year","month")
    
    # --- 3.3.2. Add missing columns
    target_gbif <- target_gbif %>% 
      mutate(scientificname = target$scientificname[1],
             worms_id = QUERY$SUBFOLDER_INFO$SP_SELECT,
             measurementvalue = "present",
             measurementunit = "Occurrence",
             taxonrank = target$taxonrank[1]) %>% 
      dplyr::select(any_of(c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","month","measurementvalue","measurementunit", "taxonrank")))
    
  }
  
  # --- 4. Bind both datasets and add nb_occ
  # Again, security if there is not GBIF data available for the species
  if(nrow(target_gbif) != 0){
     target <- rbind(target, target_gbif) %>% 
       mutate(nb_occ = n())
  } else {
    target <- target %>% 
      mutate(nb_occ = n())
  }
  
  # --- 5. Create Y target table
  Y <- target %>% 
    dplyr::select(measurementvalue)
  
  # --- 6. Create S sample table
  # Add the row ID in S to keep track of rows in sub-sample, train, test etc...
  S <- target %>% 
    dplyr::select(-any_of(c("measurementvalue", "worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           month = as.numeric(month),
           ID = row_number())
  
  # --- 7. Create an Annotation table
  annotations <- target %>% 
    dplyr::select(any_of(c("worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    distinct()
  
  # --- 8. Save in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)
  
} # END FUNCTION
