#' =============================================================================
#' @name list_bio_wrapper
#' @description wrapper around functions to list the available species within the
#' different data types and access service sources.
#' @param FOLDER_NAME name of the folder to create, corresponding to the run
#' @param DATA_SOURCE type of data to query from : string among "abundance", "occurrence" or
#' "omic" corresponding to abundance (e.g. Ecotaxa), occurrence-only (e.g. OBIS) 
#' or omics (e.g. MGnify) data; this parameter can also be a string path to a custom input table
#' @param SAMPLE_SELECT list of sample selection criteria, including :
#' MIN_SAMPLE : minimum number of geographical points to consider the selection (i.e. non-zero records)
#' MIN_DEPTH and MAX_DEPTH: minimum and maximum depth levels in m
#' START_YEAR and STOP_YEAR: start and stop years for the considered period
#' @details For custom tables, the data table must contain the following columns and one row per sample :
#' - scientificname : name of the taxa
#' - worms_id : aphiaID corresponding to the taxa name. Fill with another identified if none is available or not relevant
#' - decimallatitude : decimal, latitude of the sample (-90 to +90)
#' - decimallongitude : decimal, longitude of the sample (-180 to +180) 
#' - depth : integer, depth of the sample
#' - year : integer, year of sampling
#' - measurementvalue : numeric or string (e.g. present)
#' - measurementunit : information on the unit of the measurement 
#' - taxonrank : taxonomic ranking corresponding to the scientific name (e.g. species, gender, order...)
#' @return for "abundance" and "occurrence" data, returns a data frame with the number
#' of occurrences per worms_ID available
#' @return for "prop" data, returns a complete list of metadata, samples, taxonomic
#' annotations available.
#' @return the returned object is saved in the run file to avoid re-running the query

list_bio_wrapper <- function(FOLDER_NAME = "test_run",
                             DATA_SOURCE = "abundance",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016)){
  
  # --- 1. Initialize
  # --- 1.1. Parameter checking
  if(!is.character(DATA_SOURCE)){
    stop("The specified data source should be 'abundance', 'occurrence', 'omic' or a path to file")
  }
  
  # --- 1.2. Folder creation
  # If the directory exist, stops the creation and inform the user to manually delete the
  # previous runs, or name the directory differently
  folderpath <- paste0(project_wd,"/output/",FOLDER_NAME)
  if(file.exists(folderpath)==TRUE){
    stop("--- This foldername is already used")
  } else {
    dir.create(folderpath)
  }
  
  # --- 2. Redirection to OBIS data access
  # For abundance and occurrence source data
  if(DATA_SOURCE == "occurrence"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/01d_list_occurrence.R"))
    
    # --- 2.2. Run function
    LIST_BIO <- list_occurrence(DATA_SOURCE = DATA_SOURCE,
                                SAMPLE_SELECT = SAMPLE_SELECT)
  } # End ATLANTECO redirection
  
  # --- 3. Redirection to ATLANTECO data access
  # For abundance and occurrence source data
  if(DATA_SOURCE == "abundance"){
    # --- 3.1. Load function
    source(file = paste0(project_wd, "/code/01c_list_atlanteco.R"))
    
    # --- 3.2. Run function
    LIST_BIO <- list_atlanteco(DATA_SOURCE = DATA_SOURCE,
                               SAMPLE_SELECT = SAMPLE_SELECT)
  } # End ATLANTECO redirection
  
  # --- 4. Redirection to MGNIFY data access
  # For omic source data
  if(DATA_SOURCE == "omic"){
    # --- 4.1. Load function
    source(file = paste0(project_wd, "/code/01b_list_mgnify.R"))
    
    # --- 4.2. Run function
    LIST_BIO <- list_mgnify(SAMPLE_SELECT = SAMPLE_SELECT)
  } # End MGNIFY redirection
  
  # --- 5. Redirection to CUSTOM data access
  # For any type of data
  if(DATA_SOURCE != "omic" & DATA_SOURCE != "abundance" & DATA_SOURCE != "occurrence"){
    # --- 5.1. Load function
    source(file = paste0(project_wd, "/code/01e_list_custom.R"))
    
    # --- 5.2. Run function
    LIST_BIO <- list_custom(DATA_SOURCE = DATA_SOURCE,
                            SAMPLE_SELECT = SAMPLE_SELECT)
  } # End MGNIFY redirection
  
  # --- 6. Wrap up and save
  CALL <- list(DATA_SOURCE = DATA_SOURCE, 
               SAMPLE_SELECT = SAMPLE_SELECT, 
               LIST_BIO = LIST_BIO)
  save(CALL, file = paste0(project_wd,"/output/", FOLDER_NAME, "/CALL.RData"))
  return(CALL$LIST_BIO)
  
} # END FUNCTION
