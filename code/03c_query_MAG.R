#' =============================================================================
#' @name query_MAG
#' @description extracts biological from the MATOU data downloaded in list_bio, 
#' according to a user provided list of species, time and depth range. 
#' The extracted data is formatted to be directly usable by the models available
#'  in this workbench.
#' @param FOLDER_NAME name of the corresponding folder
#' @param QUERY the QUERY.RData object containing the list of species
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth
#' @return Annotations: a data frame of taxonomic (maybe functional) annotation
#' for each target species.
#' @return output in a QUERY object

query_MAG <- function(FOLDER_NAME = NULL,
                         QUERY = NULL){
  
  # --- 1. Initialize
  # --- 1.1. Database access
  db <- dbConnect(
    drv=PostgreSQL(),
    host="postgresql-srv.d4science.org",
    dbname="bluecloud_demo2",
    user="bluecloud_demo2_u",
    password="6a26c54a05ec5dede958a370ca744a",
    port=5432
  )
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT
  
  # --- 1.3. Get species correspondences
  # --- 1.3.1. Raw extraction
  target_input <- CALL$LIST_BIO %>% 
    dplyr::filter(scientificname %in% SP_SELECT)
  
  # --- 1.3.2. Finding the right taxonomic rank
  # Some tags are the same across different taxonomic ranks
  # We chose the rank that has the same unique values as the SP_SELECT
  taxo_rank_nb <- target_input$taxonrank %>% 
    as.factor() %>% summary()
  trank <- which(taxo_rank_nb == length(SP_SELECT)) %>% names()
  if(length(trank) == 0){
    message("QUERY_MAG: Could not find the right taxonomic rank. Please select names that correspond to a unique taxonomic rank.")
    return(NULL)
    } # end early return

  # --- 2. Raw query
  # --- 2.1. Retrieve data for the taxa of interest
  # Do this in parallel to speed up the process
  
  # Name repair on the MAGs (to remove the last "_" eventually)
  # Due to a shift in the character number, because of 2 or 3 digit station where the MAG was referenced first
  message(paste(Sys.time(), "QUERY_MAG: querying the raw data, this might take a while"))
  target <- tbl(db, "data") %>% 
    mutate(tmp = paste0(str_sub(Genes, 1, 22), "tmp")) %>% 
    mutate(MAG = str_replace(tmp, "_tmp|tmp", "")) %>% 
    dplyr::select(readCount, Station, Latitude, Longitude, Filter, Phylum, Class, Order, Family, Genus, MAG) %>% 
    dplyr::filter_at(vars(trank), any_vars(. %in% SP_SELECT)) %>%
    mutate(taxonrank = trank, 
           scientificname = !!sym(trank)) %>% # Using mutate with bang-bang operator to create a new column from variable
    collect() %>% 
    group_by(scientificname, taxonrank, Station, Filter) %>% 
    summarize(measurementvalue = sum(readCount)) %>% 
    ungroup()
  message(paste(Sys.time(), "DONE"))
  
  # --- 2.2. Normalize by the total reads per station
  # To have the same sequencing depth
  sum_station <- tbl(db, "sum_station") %>% collect()
  
  target_norm <- target %>% 
    left_join(sum_station) %>% 
    mutate(measurementvalue = measurementvalue / sum_reads) %>% 
    group_by(scientificname, taxonrank, Station) %>% 
    summarize(measurementvalue = mean(measurementvalue)) %>% 
    ungroup()
  
  # --- 2.3. Fix the column names
  locs_w_time <- tbl(db, "locs_w_time") %>% collect()
  
  target_pretty <- target_norm %>% 
    left_join(locs_w_time) %>% 
    mutate(measurementunit = "Relative proportion of metagenomic reads; size class 0.8 to 5 micrometer",
           depth = 1,
           worms_id = scientificname) %>% 
    dplyr::select(scientificname, worms_id, Latitude, Longitude, depth, year, month, measurementvalue, measurementunit, taxonrank, Station) # re-order the columns
  names(target_pretty) <- c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","month","measurementvalue","measurementunit", "taxonrank", "station")
  
  # --- 2.4. Security
  # Discard eventual misspelling in station name resulting in NA coords
  target_pretty <- target_pretty %>% 
    dplyr::filter(!is.na(decimallatitude) | !is.na(decimallongitude))
  
  # --- 3. Build S
  S <- target_pretty %>% 
    dplyr::select(-any_of(c("measurementvalue", "worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    distinct() %>% 
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           month = as.numeric(month),
           ID = row_number())
  
  # --- 4. Build Y
  # --- 4.1. Bind the data as columns
  Y <- lapply(SP_SELECT, function(x){
    y <- target_pretty %>% 
      dplyr::select(decimallatitude, decimallongitude, depth, year, month, measurementvalue, scientificname) %>% 
      dplyr::filter(scientificname == x)
  }) %>% bind_cols() %>% 
    dplyr::select(contains("measurementvalue"))
  
  colnames(Y) <- SP_SELECT # add the names again
  
  # --- 4.2. Remove rows that sum to 0
  # We duplicate the first column to be able to apply() both in multi- univariate
  to_remove <- which(apply(cbind(Y, Y[,1]), 1, sum) == 0)
  if(length(to_remove) == 1){
    Y <- Y[-to_remove,]
    S <- S[-to_remove,]
    } # end if
  
  
  # --- 4.2. Transform to proportion data
  if(CALL$DATA_TYPE == "proportions"){
    Y <- apply(Y, 1, function(x)(x = x/sum(x))) %>% 
      aperm(c(2,1)) %>% 
      as.data.frame()
  }
  
  # --- 4.3. Transform to presence data
  if(CALL$DATA_TYPE == "presence_only"){
    Y <- Y/Y %>% 
      as.data.frame()
    names(Y) <- "measurementvalue"
  }
  
  # --- 4.4. Transform to diversity
  # In form of a Shannon index for now but we would compute more
  if(CALL$DATA_TYPE == "continuous"){
    Y <- apply(Y, 1, function(x)(x = vegan::diversity(x, "shannon"))) %>% 
      as.data.frame()
    names(Y) <- "measurementvalue"
  }
  
  # --- 5. Build Annotations
  annotations <- target_pretty %>% 
    dplyr::select(any_of(c("worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    distinct()
  
  # --- 6. Save in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  # --- 7. Disconnect from database
  dbDisconnect(db)
  
  return(QUERY)
  
} # END FUNCTION