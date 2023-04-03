#' =============================================================================
#' @name query_bio
#' @description extracts biological data from a user defined dataset among those
#' available within the Blue Cloud data access service. The extracted data is
#' formatted to be directly usable by the models available in this workbench.
#' @param DATA_TYPE string among "cont", "pres" or "prop" corresponding to
#' continuous (e.g. Ecotaxa), presence-only (e.g. OBIS) or proportions (e.g. 
#' omics) data.
#' @param SP_SELECT integer, corresponding to the worms Aphia ID to extract
#' @param SAMPLE_SELECT list of sample selection criteria, including :
#' MIN_SAMPLE : minimum number of geographical points to consider the selection
#' MIN_DEPTH and MAX_DEPTH: minimum and maximum depth levels in m
#' START_YEAR and STOP_YEAR: start and stop years for the considered period
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth

query_bio <- function(DATA_TYPE = "cont",
                      SP_SELECT = 1,
                      SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016)){
  
  # ========================== CONNECT TO DATABASE =============================
  # For now only the presence and abundance data from AtlantECO are available
  # in a temporary SQLite database because working on the data access service
  
  db <- dbConnect(RSQLite::SQLite(), paste0(project_wd,"/data/DB_clean.sqlite"))
  
  # =========================== PARAMETER CHECK ================================
  if(DATA_TYPE != "cont" & DATA_TYPE != "pres" & DATA_TYPE != "omic"){
    stop("The specified data type should be 'cont', 'pres' or 'omic'")
  } else if(is.numeric(SP_SELECT) == FALSE){
    stop("SP_SELECT has to be a numeric value of WoRMS_ID or a NULL object")
  }

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
    ungroup() %>%  
    dplyr::filter(nb_occ >= !!SAMPLE_SELECT$MIN_SAMPLE)
  
  # --- 4. Create Y target table
  Y <- target %>% 
    dplyr::select(measurementvalue)
  
  # --- 5. Create S sample table
  # Add the row ID in S to keep track of rows in sub-sample, train, test etc...
  S <- target %>% 
    dplyr::select(-measurementvalue) %>% 
    mutate(ID = row_number())
  
  # --- 6. Write Y and S on the disk
  write_feather(Y, paste0(project_wd, "/data/Y.feather"))
  write_feather(S, paste0(project_wd, "/data/S.feather"))
  
  return(list(Y = Y, S = S, CALL = list(DATA_TYPE = DATA_TYPE,
                                        SP_SELECT = SP_SELECT,
                                        SAMPLE_SELECT = SAMPLE_SELECT)))

} # END FUNCTION