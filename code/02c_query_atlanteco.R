#' =============================================================================
#' @name query_atlanteco
#' @description extracts biological from the Atlanteco database, according to a
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

query_atlanteco <- function(FOLDER_NAME = NULL,
                            QUERY = NULL){
  
  # --- 1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT
  
  # --- 2. Connect to database
  # For now only the presence and abundance data from AtlantECO are available
  # in a temporary SQLite database because working on the data access service
  db <- dbConnect(RSQLite::SQLite(), paste0(project_wd,"/data/DB_clean.sqlite"))
  
  # --- 3. Extract data from our Aphia_ID of interest
  target <- tbl(db, paste0(CALL$DATA_TYPE, "_data")) %>% 
    dplyr::filter(worms_id %in% !!CALL$SP_SELECT) %>% 
    collect() %>% 
    mutate(month = str_pad(month, 2, pad = "0"))
  
  # --- 4. Filter target according to SAMPLE_SELECT requirements
  target <- target %>% 
    dplyr::filter(depth >= !!CALL$SAMPLE_SELECT$MIN_DEPTH & 
                    depth <= !!CALL$SAMPLE_SELECT$MAX_DEPTH & 
                    year >= !!CALL$SAMPLE_SELECT$START_YEAR & 
                    year <= !!CALL$SAMPLE_SELECT$STOP_YEAR &
                    measurementvalue != "Absence") %>% 
    group_by(worms_id) %>% 
    mutate(nb_occ = n()) %>% 
    ungroup()
  
  # --- 5. Create Y target table
  Y <- target %>% 
    dplyr::select(measurementvalue)
  
  # --- 6. Create S sample table
  # Add the row ID in S to keep track of rows in sub-sample, train, test etc...
  S <- target %>% 
    dplyr::select(-measurementvalue, -worms_id, -taxonrank, -scientificname, -nb_occ) %>% 
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           ID = row_number())
  
  # --- 7. Create an Annotation table
  annotations <- target %>% 
    dplyr::select(worms_id, taxonrank, scientificname, nb_occ) %>% 
    distinct()
  
  # --- 8. Disconnect from database
  dbDisconnect(db)
  
  # --- 9. Save in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)
  
} # END FUNCTION
