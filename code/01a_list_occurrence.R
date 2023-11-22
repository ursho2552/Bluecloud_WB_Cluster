#' =============================================================================
#' @name list_occurrence
#' @description extracts available Aphia_ID corresponding to user defined criteria
#' among the available data within the OBIS and GBIF database
#' @param DATA_SOURCE parameter passed from the wrapper function
#' @param SAMPLE_SELECT parameter passed from the wrapper function
#' @return a list of available Worms ID or Aphia ID and number of occurrences
#' within the data type and sample criteria

list_occurrence <- function(DATA_SOURCE,
                            SAMPLE_SELECT){
  
  # --- 1. List the available taxa
  # Returns the Worms ID, or Aphia ID and number of occurrences from OBIS
  message("--- LIST BIO : retrieving species occurrences available in OBIS and GBIF")
  
  # --- 1.1. Raw query on OBIS
  # Keeping species with more than 10 records (to speed up the process). 
  # It will be filtered to MIN_SAMPLE after GBIF query
  data_list <- checklist(startdate = as.Date(as.character(SAMPLE_SELECT$START_YEAR), format = "%Y"),
                        enddate = as.Date(as.character(SAMPLE_SELECT$STOP_YEAR), format = "%Y"),
                        startdepth = SAMPLE_SELECT$TARGET_MIN_DEPTH,
                        enddepth = SAMPLE_SELECT$TARGET_MAX_DEPTH) %>% 
    dplyr::filter(records >= 2) %>% 
    dplyr::select(acceptedNameUsageID, taxonRank, acceptedNameUsage, records) %>% 
    distinct()
  
  # --- 1.2. Nice names
  colnames(data_list) <- c("worms_id","taxonrank","scientificname","obis_occ")
  
  # --- 1.3. Remove spurious names from OBIS (e.g., including NA, brackets etc...)
  # There has been a OBIS update and now some names contain brackets...
  data_list <- data_list %>% 
    dplyr::filter(!grepl("\\[", scientificname)) # double backslash to escape the meaning of "["
  
  # --- 2. Complete with GBIF
  # --- 2.1. Retrieve occurrence number
  # We cannot apply year and depth filter here. It will be done later in the query check step
  message(paste(Sys.time(), "--- Retriving GBIF data : START"))
  gbif_occ <- mclapply(data_list$scientificname,
                     FUN = function(NAME){
                       occ_count(scientificName = NAME,
                                 year = paste0(SAMPLE_SELECT$START_YEAR, ",", SAMPLE_SELECT$STOP_YEAR),
                                 depth = paste0(SAMPLE_SELECT$TARGET_MIN_DEPTH, ",", SAMPLE_SELECT$TARGET_MAX_DEPTH),
                                 occurrenceStatus = 'PRESENT')
                     },
                     mc.cores = MAX_CLUSTERS) %>% unlist()
  message(paste(Sys.time(), "--- Retriving GBIF data : DONE"))
  
  # --- 2.2. Adding it to the initial occurrence number
  data_list <- data.frame(data_list,
                          gbif_occ = gbif_occ,
                          nb_occ = data_list$obis_occ + gbif_occ) %>% 
    dplyr::filter(nb_occ >= SAMPLE_SELECT$MIN_SAMPLE)

  # --- 3. Wrap up
  return(data_list)
  
} # END FUNCTION
