#' =============================================================================
#' @name list_MAG
#' @description extracts available metadata, taxonomy and read counts per MAG
#' corresponding to user defined criteria, among the available data within the
#' MATOU data access service - on Bluecloud Plankton Genomics VLAB.
#' @param SAMPLE_SELECT parameter passed from the wrapper function
#' @return a complete list of metadata, samples, taxonomic
#' annotations available.

list_MAG <- function(SAMPLE_SELECT){

  # --- 1. Initialize
  # --- 1.1. Database connection
  # Connect to the public PostgreSQL database
  db <- dbConnect(
    drv=PostgreSQL(),
    host="postgresql-srv.d4science.org",
    dbname="bluecloud_demo2",
    user="bluecloud_demo2_u",
    password="6a26c54a05ec5dede958a370ca744a",
    port=5432
  )

  # --- 2. Raw query
  # --- 2.1. Extract the MAG data
  message(paste(Sys.time(), "LIST_MAGS: Start the raw query"))
  data_w_taxo <- tbl(db, "data") %>%
    mutate(MAG = str_sub(Genes, 1, 22)) %>%
    dplyr::filter(readCount > 0) %>%
    dplyr::select("Station","Phylum","Class","Order","Family","Genus","MAG") %>%
    distinct() %>%
    collect()
  message(paste(Sys.time(), "DONE"))

  # Name repair on the MAGs (to remove the last "_" eventually)
  # Due to a shift in the character number, because of 2 or 3 digit station where the MAG was referenced first
  tmp <- paste0(data_w_taxo$MAG, "tmp")
  data_w_taxo$MAG <- str_replace(tmp, "_tmp|tmp", "")

  # --- 2.2. Filter samples according to user input
  # Opening sample information and filter stations accordingly
  message(paste(Sys.time(), "LIST_MAGS: Filter stations according to defined criteria"))

  locs_w_time <- tbl(db, "locs_w_time") %>% collect() %>%
    mutate(depth = 1) %>% # MATOU data are retrieved from Tara Ocean Surface samples
    dplyr::filter(year >= SAMPLE_SELECT$START_YEAR & year <= SAMPLE_SELECT$STOP_YEAR) %>%
    dplyr::filter(depth >= SAMPLE_SELECT$TARGET_MIN_DEPTH & depth <= SAMPLE_SELECT$TARGET_MAX_DEPTH) %>%
    dplyr::select(Station)

  data_w_taxo <- data_w_taxo %>%
    inner_join(locs_w_time)

  # --- 3. Extract the taxonomic information and observations available
  # --- 3.1. Extract and format the different taxonomic ranks
  message(paste(Sys.time(), "LIST_MAGS: Formatting the output and taxonomy"))

  list_raw <- lapply(c("Phylum","Class","Order","Family","Genus","MAG"), function(x){
    id <- which(colnames(data_w_taxo) == x)
    df <- data_w_taxo %>%
      dplyr::select(all_of(1:id)) %>%
      mutate(taxonrank = x,
             scientificname = data_w_taxo[[x]])
  }) %>% bind_rows() %>%
    distinct()

  # --- 3.2. Count the number of observations for each combination of taxonomic rank and scientific name
  # Filter by number of observations and organize the data into a tidy format.
  message(paste(Sys.time(), "LIST_MAGS: Filter samples by nb. obs."))

  list_bio <- list_raw %>%
    group_by(taxonrank, scientificname) %>%
    mutate(nb_occ = n()) %>%
    ungroup() %>%
    dplyr::select(scientificname, taxonrank, nb_occ, Phylum, Class, Order, Family, Genus, MAG) %>%
    dplyr::filter(!is.na(scientificname)) %>%
    dplyr::filter(nb_occ >= SAMPLE_SELECT$MIN_SAMPLE) %>%
    distinct()

  # --- 4. Disconnect from database
  dbDisconnect(db)

  # --- 5. Wrap up and save
  return(list_bio)

} # END FUNCTION

