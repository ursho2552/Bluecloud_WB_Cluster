#' =============================================================================
#' @name query_abundance_biomass
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

query_abundance_biomass <- function(CALL,
                                    FOLDER_NAME = NULL,
                                    QUERY = NULL){

  # --- 1. Parameter loading
  # load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  # --- 1.1. Single query if no WORMS_CHECK
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT
  # --- 1.2. Multiple query if WORMS_CHECK = TRUE; i.e., also requesting children and synonyms
  if(CALL$WORMS_CHECK == TRUE){
    SP_SELECT <- CALL$SP_SELECT_INFO[[SP_SELECT]] %>% unlist() %>% as.numeric() %>% .[!is.na(.)] # take the vector of species to do a common query across all children and synonyms
  }

  # --- 2. Connect to database
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

  # --- 3. Extract data from our Aphia_ID of interest - abundance or biomass depending on the source
  target <- tbl(db, paste0(CALL$DATA_SOURCE, "_data")) %>%
    dplyr::filter(worms_id %in% !!SP_SELECT) %>%
    collect() %>%
    mutate(month = str_pad(month, 2, pad = "0"))

  # --- 4. Filter target according to SAMPLE_SELECT requirements
  # --- 4.1. Default univariate target
  target <- target %>%
    dplyr::filter(depth >= !!CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH &
                    depth <= !!CALL$SAMPLE_SELECT$TARGET_MAX_DEPTH &
                    year >= !!CALL$SAMPLE_SELECT$START_YEAR &
                    year <= !!CALL$SAMPLE_SELECT$STOP_YEAR &
                    measurementvalue != 0 & !is.na(measurementvalue)) %>%
    group_by(worms_id) %>%
    ungroup() %>%
    distinct() %>%
    mutate(nb_occ = n())

  # --- 4.2. Expand by target if DATA_TYPE = "proportions"
  target_proportions <- target %>%
    dplyr::select("decimallatitude","decimallongitude", "depth","year","month","measurementvalue","worms_id") %>%
    pivot_wider(names_from = "worms_id", values_from = "measurementvalue", values_fn = mean)

  # --- 5. Create Y target table
  # --- 5.1. For one target; i.e. DATA_TYPE = "continuous"
  Y <- target %>%
    dplyr::select(measurementvalue)
  # --- 5.2. For several targets; i.e. DATA_TYPE = "proportions"
  if(CALL$DATA_TYPE == "proportions"){
    Y <- target_proportions %>%
      dplyr::select(-decimallatitude, -decimallongitude, -depth, -year, -month)
    Y[is.na(Y)] <- 0
    Y <- apply(Y, 1, function(x)(x = x/sum(x))) %>% aperm(c(2,1)) %>% as.data.frame()
  }

  # --- 6. Create S sample table
  # --- 6.1. For one target; i.e. DATA_TYPE = "continuous"
  # Add the row ID in S to keep track of rows in sub-sample, train, test etc...
  S <- target %>%
    dplyr::select(-measurementvalue, -worms_id, -taxonrank, -scientificname, -nb_occ) %>%
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           month = as.numeric(month),
           ID = row_number())

  # --- 6.2. For several targets; i.e. DATA_TYPE = "proportions"
  if(CALL$DATA_TYPE == "proportions"){
    S <- target %>%
      dplyr::select(-measurementvalue, -worms_id, -taxonrank, -scientificname, -nb_occ) %>%
      unique() %>%
      mutate(decimallatitude = as.numeric(decimallatitude),
             decimallongitude = as.numeric(decimallongitude),
             month = as.numeric(month),
             ID = row_number())
  }

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
