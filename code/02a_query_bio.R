#' =============================================================================
#' @name query_bio
#' @description extracts biological data from a user defined dataset among those
#' available within the Blue Cloud data access service. The extracted data is
#' formatted to be directly usable by the models available in this workbench.
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth
#' @return Saves the output in a QUERY.RData file

query_bio <- function(SP_SELECT = NULL,
                      FOLDER_NAME = NULL){
  
  # ========================== CONNECT TO DATABASE =============================
  # For now only the presence and abundance data from AtlantECO are available
  # in a temporary SQLite database because working on the data access service
  
  db <- dbConnect(RSQLite::SQLite(), paste0(project_wd,"/data/DB_clean.sqlite"))
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  DATA_TYPE <- CALL$DATA_TYPE
  SAMPLE_SELECT <- CALL$SAMPLE_SELECT
  
  # ============================== DATA QUERY ==================================
  # --- 1. Extract data from our Aphia_ID of interest
  target <- tbl(db, paste0(DATA_TYPE, "_data")) %>% 
    dplyr::filter(worms_id %in% !!SP_SELECT) %>% 
    collect() %>% 
    mutate(month = str_pad(month, 2, pad = "0"))
  
  # --- 2. Filter target according to SAMPLE_SELECT requirements
  target <- target %>% 
    dplyr::filter(depth >= !!SAMPLE_SELECT$MIN_DEPTH & 
                    depth <= !!SAMPLE_SELECT$MAX_DEPTH & 
                    year >= !!SAMPLE_SELECT$START_YEAR & 
                    year <= !!SAMPLE_SELECT$STOP_YEAR &
                    measurementvalue != "Absence") %>% 
    group_by(worms_id) %>% 
    mutate(nb_occ = n()) %>% 
    ungroup()
  
  # --- 4. Create Y target table
  Y <- target %>% 
    dplyr::select(measurementvalue)
  
  # --- 5. Create S sample table
  # Add the row ID in S to keep track of rows in sub-sample, train, test etc...
  S <- target %>% 
    dplyr::select(-measurementvalue) %>% 
    mutate(ID = row_number())

  # --- 6. Write all objects in the species folder
  QUERY <- list(Y = Y, S = S, CALL = list(DATA_TYPE = DATA_TYPE,
                                          SP_SELECT = SP_SELECT,
                                          SAMPLE_SELECT = SAMPLE_SELECT))
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME, "/", SP_SELECT,"/QUERY.RData"))
  
} # END FUNCTION
