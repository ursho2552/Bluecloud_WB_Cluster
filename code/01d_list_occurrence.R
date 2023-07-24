#' =============================================================================
#' @name list_occurrence
#' @description extracts available Aphia_ID corresponding to user defined criteria
#' among the available data within the OBIS database
#' @param DATA_SOURCE parameter passed from the wrapper function
#' @param SAMPLE_SELECT parameter passed from the wrapper function
#' @return a list of available Worms ID or Aphia ID and number of occurrences
#' within the data type and sample criteria

list_occurrence <- function(DATA_SOURCE,
                            SAMPLE_SELECT){
  
  # --- 1. List the available taxa
  # Returns the Worms ID, or Aphia ID and number of occurrences from OBIS
  message("--- LIST BIO : retrieving species occurrences available in OBIS")
  
  # --- 1.1. Raw query
  data_list <- checklist(startdate = as.Date(as.character(SAMPLE_SELECT$START_YEAR), format = "%Y"),
                        enddate = as.Date(as.character(SAMPLE_SELECT$STOP_YEAR), format = "%Y"),
                        startdepth = SAMPLE_SELECT$MIN_DEPTH,
                        enddepth = SAMPLE_SELECT$MAX_DEPTH) %>% 
    dplyr::filter(records >= SAMPLE_SELECT$MIN_SAMPLE) %>% 
    dplyr::select(taxonID, taxonRank, scientificName, records)
  
  # --- 1.2. Nice names
  colnames(data_list) <- c("worms_id","taxonrank","scientificname","nb_occ")
  
  # --- 1.3. Wrap up
  return(data_list)
  
} # END FUNCTION
