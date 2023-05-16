#' =============================================================================
#' @name list_mgnify
#' @description extracts available metadata, taxonomy and read counts per OTU
#' corresponding to user defined criteria, among the available data within the 
#' MGnify data access service
#' @param SAMPLE_SELECT parameter passed from the wrapper function
#' @return a complete list of metadata, samples, taxonomic
#' annotations available.

list_mgnify <- function(SAMPLE_SELECT){
  
  # --- 1. Initialize
  # Create MGnify cache
  mg <- mgnify_client(usecache = T, cache_dir = paste0(project_wd, "/data/.mgnify_cache"))
  
  # --- 2. List the available studies
  message("--- LIST BIO : retrieving studies available in MGnify")
  study_list <- mgnify_query(client = mg,
                             qtype = "studies",
                             biome_name = "Marine") %>% 
    mutate(`samples-count` = as.numeric(`samples-count`)) %>% 
    dplyr::filter(`samples-count` >= SAMPLE_SELECT$MIN_SAMPLE) # Studies with less than the minimum sample per OTU are out
  
  # --- 3. Download the raw data
  # --- 3.1. Analysis accession ID
  # Analysis identifier used to retrieve metadata and phyloseq object later
  message("--- LIST BIO : retrieving accession ID from MGnify")
  analysis_accession_ID <- mgnify_analyses_from_studies(mg, study_list$accession[c(1,2,4)]) # SUBSET FOR EXAMPLE
  
  # --- 3.1.2. Retrieve metadata for each accession ID
  # We only keep the metagenomic data and latest pipeline version
  message("--- LIST BIO : retrieving studies metadata from MGnify")
  sample_metadata <- mgnify_get_analyses_metadata(mg, analysis_accession_ID) %>% 
    data.frame() %>% 
    dplyr::filter(`analysis_experiment.type` == "metagenomic" |
                    `analysis_experiment.type` == "assembly" &
                    `analysis_pipeline.version` == "5.0")
  
  # --- 3.1.3. Retrieve the corresponding phyloseq objects
  message("--- LIST BIO : retrieving studies phyloseq from MGnify")
  sample_ps <- mgnify_get_analyses_phyloseq(mg, sample_metadata$analysis_accession)
  
  # --- 4. Extract read, sample and annotation tables
  # --- 4.1. Sample location description (S)
  # Directly filter by the user selected time, depth criteria
  S <- sample_ps@sam_data %>% 
    data.frame() %>% 
    rownames_to_column(var = "sample") %>% 
    mutate(source_tbl = analysis_accession,
           decimallatitude = sample_latitude,
           decimallongitude = sample_longitude,
           depth = sample_depth,
           year = format(as.Date(`sample_collection.date`, format = "%Y-%m-%d"), "%Y"),
           month = format(as.Date(`sample_collection.date`, format = "%Y-%m-%d"), "%m"),
           measurementtype = rep("proportions", n()),
           measurementunit = rep("OTU reads", n()),
           lower_size = `sample_size.fraction.lower.threshold`,
           upper_size = `sample_size.fraction.upper.threshold`,
           ID = row_number()) %>% 
    dplyr::select(source_tbl, decimallatitude, decimallongitude, depth, year, month, measurementtype, measurementunit, lower_size, upper_size, ID) %>% 
    dplyr::filter(depth >= !!SAMPLE_SELECT$MIN_DEPTH & 
                    depth <= !!SAMPLE_SELECT$MAX_DEPTH & 
                    year >= !!SAMPLE_SELECT$START_YEAR & 
                    year <= !!SAMPLE_SELECT$STOP_YEAR)
  
  # --- 4.2. Raw reads per sample x target table (Y)
  # Not yet processed as proportions. It will be done at the query stage.
  Y <- sample_ps@otu_table %>% 
    t() %>% data.frame() %>% 
    slice(S$ID)
  names(Y) <- sub("^X","", names(Y))
  
  # --- 4.3. Taxonomic and functional annotations of targets
  # Basic taxonomic annotations at the OTU level
  annotations <- sample_ps@tax_table %>% 
    data.frame() %>% 
    rownames_to_column(var = "OTU") %>% 
    mutate(nb_occ = apply(Y, 2, sum),
           nb_station = apply(Y, 2, function(x)(x = length(which(x > 0)))))
  
  # --- 4.4. Concatenate in a data list
  data_list <- list(analysis_accession_ID = analysis_accession_ID,
                    sample_metadata = sample_metadata,
                    sample_ps = sample_ps,
                    Y = Y,
                    S = S,
                    annotations = annotations)
  return(data_list)
  
} # END FUNCTION

