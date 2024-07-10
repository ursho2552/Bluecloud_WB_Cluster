#' =============================================================================
#' @name list_abundance_biomass
#' @description extracts available Aphia_ID corresponding to user defined criteria
#' among the available data within the Atlanteco database, build locally.
#' @param DATA_SOURCE parameter passed from the wrapper function
#' @param SAMPLE_SELECT parameter passed from the wrapper function
#' @return a list of available Worms ID or Aphia ID and number of occurrences
#' within the data type and sample criteria

list_abundance_biomass <- function(DATA_SOURCE,
                                   SAMPLE_SELECT){

  # --- 1. Connect to database
  # This database is a light online copy of the AtlantECO base v1
  # Further details are available in the source material referenced in the database
  db <- dbConnect(
    drv=PostgreSQL(),
    host="postgresql-srv.d4science.org",
    dbname="bc2026_wb3_db",
    user="bluecloud_wb3_reader",
    password="1a6414f89d8a265c8bdd",
    port=5432
  )

  # --- 2. Query the database - abundance or biomass depending on the source
  # Returns the Worms ID, or Aphia ID and number of occurrences within the data
  # type and sample criteria
  message("--- LIST BIO : retrieving species available in ATLANTECO")
  data_list <- tbl(db, paste0(DATA_SOURCE,"_data")) %>%
    dplyr::filter(depth >= !!SAMPLE_SELECT$TARGET_MIN_DEPTH &
                  depth <= !!SAMPLE_SELECT$TARGET_MAX_DEPTH &
                  year >= !!SAMPLE_SELECT$START_YEAR &
                  year <= !!SAMPLE_SELECT$STOP_YEAR &
                  measurementvalue != 0 & !is.na(measurementvalue),
                  worms_id != "Not found") %>%
    group_by(worms_id) %>%
    mutate(nb_occ = n()) %>%
    dplyr::select(-c(decimallongitude, decimallatitude, month, depth, year, measurementvalue, measurementunit)) %>%
    distinct() %>%
    dplyr::filter(nb_occ >= !!SAMPLE_SELECT$MIN_SAMPLE) %>%
    collect()

  # --- 3. Disconnect from database
  dbDisconnect(db)

  # --- 4. Wrap up and save
  return(data_list)

} # END FUNCTION
