#' =============================================================================
#' @name list_custom
#' @description extracts available Aphia_ID corresponding to user defined criteria
#' among the available data within the user-provided custom database
#' @param DATA_SOURCE parameter passed from the wrapper function
#' @param SAMPLE_SELECT parameter passed from the wrapper function
#' @return a list of available Worms ID or Aphia ID and number of occurrences
#' within the data type and sample criteria

list_custom <- function(DATA_SOURCE,
                        SAMPLE_SELECT){
  
  # --- 1. Open the file
  # --- 1.1. Check extension type
  ext <- tools::file_ext(DATA_SOURCE)
  
  # --- 1.2. Conditional opening
  if(ext == "xlsx"){
    df <- readxl::read_xlsx(DATA_SOURCE)
    colnames(df) <- tolower(colnames(df))
  } else if(ext == "txt" | ext == "csv"){
    df <- vroom(DATA_SOURCE)
    colnames(df) <- tolower(colnames(df))
  } else {
    message("The file extension is not recognized. It should be either .xlsx, .txt or .csv")
    return(NULL)
  }
  
  # --- 2. Check for column presence
  # --- 2.1. Mandatory columns
  names_qc <- c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","month","measurementvalue","measurementunit", "taxonrank")
  # --- 2.2. Extract corresponding names in the data
  names_df <- df %>% 
    dplyr::select(any_of(names_qc)) %>% 
    names()
  # --- 2.3. Inform the user
  if(length(names_qc) != length(names_df)){
    message("The data table must contain the following columns and one row per sample : \n
        - scientificname : name of the taxa \n
        - worms_id : aphiaID corresponding to the taxa name. Fill with another identified if none is available or not relevant \n
        - decimallatitude : decimal, latitude of the sample (-90 to +90) \n
        - decimallongitude : decimal, longitude of the sample (-180 to +180) \n
        - depth : integer, depth of the sample \n
        - year : integer, year of sampling \n
        - month : integer, month of sample \n
        - measurementvalue : numeric or string (e.g. present) \n
        - measurementunit : information on the unit of the measurement \n
        - taxonrank : taxonomic ranking corresponding to the scientific name (e.g. species, gender, order...)")
    
    return(NULL)
  }
  
  # --- 3. List of available data
  data_list <- df %>% 
    dplyr::filter(depth >= SAMPLE_SELECT$TARGET_MIN_DEPTH & depth <= SAMPLE_SELECT$TARGET_MAX_DEPTH) %>% 
    dplyr::filter(year >= SAMPLE_SELECT$START_YEAR & year <= SAMPLE_SELECT$STOP_YEAR) %>% 
    dplyr::filter(!is.na(measurementvalue)) %>% 
    distinct() %>% 
    group_by(scientificname) %>% 
    mutate(nb_occ = n()) %>% 
    ungroup()
  
  # --- 4. Wrap up
  return(data_list)
  
} # end function
