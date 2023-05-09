#' =============================================================================
#' @name folder_init
#' @description initialize an output folder corresponding to the run call (i.e., species
#' and selection criteria) in which the output will be saved
#' @param DATA_TYPE string among "cont", "pres" or "prop" corresponding to
#' continuous (e.g. Ecotaxa), presence-only (e.g. OBIS) or proportions (e.g. 
#' omics) data
#' @param SAMPLE_SELECT list of sample selection criteria, including :
#' MIN_SAMPLE : minimum number of geographical points to consider the selection (i.e. non-zero records)
#' MIN_DEPTH and MAX_DEPTH: minimum and maximum depth levels in m
#' START_YEAR and STOP_YEAR: start and stop years for the considered period
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @return creates a folder with the call written in a file for traceback

folder_init <- function(DATA_TYPE = "cont",
                        SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016),
                        SP_SELECT = NULL,
                        FOLDER_NAME = "test_run"){

  # --- 1. Create directory
  # If the directory exist, stops the creation and inform the user to manually delete the
  # previous runs, or name the directory differently
  folderpath <- paste0(project_wd,"/output/",FOLDER_NAME)
  if(file.exists(folderpath)==TRUE){
    stop("--- This foldername is already used")
  } else {
    dir.create(folderpath)
  }
  
  # --- 2. Create the call file
  # Create a RData file in the folder root informing the general call parameters
  CALL <- list(DATA_TYPE = DATA_TYPE, 
               SAMPLE_SELECT = SAMPLE_SELECT, 
               SP_SELECT = SP_SELECT)
  save(CALL, file = paste0(folderpath, "/CALL.RData"))
  
  # --- 3. Create species sub-directories
  # Named by their ID, to be used as a list of species argument later
  for(i in SP_SELECT){
    dir.create(paste0(folderpath, "/", i))
  }
  
} # END FUNCTOIN

