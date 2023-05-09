#' =============================================================================
#' @name list_bio
#' @description extracts available Aphia_ID data from a user defined dataset
#' among those available within the Blue Cloud data access service.
#' @param DATA_TYPE string among "cont", "pres" or "prop" corresponding to
#' continuous (e.g. Ecotaxa), presence-only (e.g. OBIS) or proportions (e.g. 
#' omics) data
#' @param SAMPLE_SELECT list of sample selection criteria, including :
#' MIN_SAMPLE : minimum number of geographical points to consider the selection (i.e. non-zero records)
#' MIN_DEPTH and MAX_DEPTH: minimum and maximum depth levels in m
#' START_YEAR and STOP_YEAR: start and stop years for the considered period
#' @return worms_ID: a data frame of Aphia ID available

list_bio <- function(DATA_TYPE = "cont",
                     SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016)){
  
  # ========================== CONNECT TO DATABASE =============================
  # For now only the presence and abundance data from AtlantECO are available
  # in a temporary SQLite database because working on the data access service
  
  db <- dbConnect(RSQLite::SQLite(), paste0(project_wd,"/data/DB_clean.sqlite"))
  
  # =========================== PARAMETER CHECK ================================
  if(DATA_TYPE != "cont" & DATA_TYPE != "pres" & DATA_TYPE != "omic"){
    stop("The specified data type should be 'cont', 'pres' or 'omic'")
  }
  
  # ========================= LAZY DATA QUERY ==================================
  # If no Aphia ID is specific, returns the ID available 
  # within the data type and minimum sample criteria
  worms_list <- tbl(db, paste0(DATA_TYPE,"_data")) %>% 
    dplyr::select(worms_id, decimallongitude, decimallatitude, depth, year, measurementvalue) %>% 
    dplyr::filter(depth >= !!SAMPLE_SELECT$MIN_DEPTH & 
                    depth <= !!SAMPLE_SELECT$MAX_DEPTH & 
                    year >= !!SAMPLE_SELECT$START_YEAR & 
                    year <= !!SAMPLE_SELECT$STOP_YEAR &
                    measurementvalue != "Absence" &
                    measurementvalue != 0) %>% 
    group_by(worms_id) %>% 
    summarise(nb_occ = n()) %>% 
    dplyr::filter(nb_occ >= !!SAMPLE_SELECT$MIN_SAMPLE) %>% 
    collect()
  return(list(worms_list = worms_list, CALL = list(DATA_TYPE = DATA_TYPE,
                                                   SAMPLE_SELECT = SAMPLE_SELECT)))


} # END FUNCTION
