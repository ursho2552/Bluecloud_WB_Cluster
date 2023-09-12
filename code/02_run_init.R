#' =============================================================================
#' @name run_init
#' @description (1) initialize an output folder corresponding to the run call (i.e., species
#' and selection criteria) in which the output will be saved.
#' @description (2) initialize all parameters of the run and stores them in
#' a CALL object

#' @param FOLDER_NAME name of the run folder we want to work in
#' @param SP_SELECT vector of IDs, corresponding to the species to parallelize on
#' @param FAST TRUE or FALSE; if TRUE, does not compute projections and plot for algorithms
#' that did not pass the Quality Checks
#' @param LOAD_FROM load a previous list_bio object from another folder to be
#' duplicated in the new FOLDER_NAME. It avoids re-doing the initial list_bio step
#' that can be long for omics data

#' @param SP_SELECT species to run the analysis for, in form of Aphia ID for traditional data
#' and OTU ID for omics
#' @param DATA_TYPE the output type of the data, which can influence the sub-folder
#' architecture. See details.

#' @param ENV_VAR a two column data frame containing (i) the folder name corresponding
#' to each environmental variable; (ii) the corresponding variable name in the .nc files
#' @param ENV_PATH string or vector of path to the root where the folders containing
#' the .nc are.
#' @param METHOD_PA method of pseudo-absence, either "mindist" or "cumdist" or "density"
#' @param NB_PA number of pseudo-absences to generate
#' @param PER_RANDOM ratio of pseudo-absences that are sampled randomly in the background
#' @param DIST_PA if METHOD_PA = "mindist", distance from presences (in meters),
#'  from which to define the background data. Expert use only.
#' @param BACKGROUND_FILTER additional background filter for finer tuning, such
#' as selecting pseudo-absences within the sampled background of a given campaign
#' or instrument deployment. Passed by the user in the form of a 2 column 
#' data frame, x = longitude and y = latitude where the pseudo-absences
#' can be sampled. Or a path to a raster object where pseudo-absences are sampled in
#' non NA cells, weighted by the cell values.

#' @param OUTLIER if TRUE, remove outliers
#' @param UNIVARIATE if true, performs a univariate predictor pre-selection
#' @param ENV_COR numeric, removes the correlated environmental values from the
#' query objects and CALL according to the defined threshold. Else NULL.

#' @param NFOLD number of folds, used defined integer
#' @param FOLD_METHOD method used to create the folds, integer between "kfold"
#' and "lon"

#' @param MODEL_LIST list of algorithms from which to compute hyperparameter
#' selection
#' @param LEVELS maximum number of parameter values to test in each of the grids

#' @param ENSEMBLE TRUE or FALSE; if TRUE, computes an ensemble at the evaluation and projection steps
#' @param N_BOOTSTRAP number of bootstrap to do for the projections and partial dependency plots
#' @param CUT numeric or NULL; if numeric, level at which the projections are considered to be 0.
#' Projection patches without observation are then removed.
#' @param PROJ_PATH (optional) path to a environmental raster, potentially 
#' different than the one given in the QUERY object. This is the case for 
#' supplementary projections in other time and climate scenarios for example. 
#' To your own risk and only for expert users !

#' @details Different data transformation between DATA_SOURCE and DATA_TYPE are implemented, including:
#' - "occurrence" to "binary" (default)
#' - "abundance" to "continuous" (default, not recommended for more than 50 targets)
#' - "abundance" to "proportions" (not recommended if the sampling stations are not the same)
#' - "omic" to "proportions" (default, not recommended for more than 50 targets)
#' - "omic" to "continuous" in form of richness
#' - "omic" to "binary" in form of presence-only
#' @return creates one subfolder per SP_SELECT and a vector of sub-directories names 
#' used as an argument for parallel computing later
#' @return all global parameters in a CALL.RData object

run_init <- function(FOLDER_NAME = "test_run",
                     SP_SELECT = NULL,
                     FAST = FALSE,
                     LOAD_FROM = NULL,
                     DATA_TYPE = NULL,
                     ENV_VAR = data.frame(name = c("oxygen_monthly_WOA18", 
                                                   "salinity_monthly_WOA18", 
                                                   "mld_monthly_WOA18", 
                                                   "silicate_monthly_WOA18",
                                                   "nitrate_monthly_WOA18",
                                                   "phosphate_monthly_WOA18",
                                                   "temp_monthly_WOA18"),
                                          ncvar = c("o_an", "s_an", "M_an", "i_an", "n_an", "p_an", "t_an")),
                     ENV_PATH = "/net/meso/work/nknecht/Masterarbeit/Data/21_10_18_environmental_data",
                     METHOD_PA = "density",
                     NB_PA = NULL,
                     PER_RANDOM = 0.25,
                     DIST_PA = NULL,
                     BACKGROUND_FILTER = NULL,
                     OUTLIER = TRUE,
                     UNIVARIATE = TRUE,
                     ENV_COR = 0.8,
                     NFOLD = 3,
                     FOLD_METHOD = "lon",
                     MODEL_LIST = c("GLM","GAM","RF","MLP","SVM","BRT"),
                     LEVELS = 3,
                     ENSEMBLE = TRUE,
                     N_BOOTSTRAP = 10,
                     CUT = 0.1,
                     PROJ_PATH = NULL){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - create file
  if(is.null(LOAD_FROM)){
    sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/log.txt"), open = "wt"),
                         START = TRUE)
    message(paste(Sys.time(), "******************** START : run_init ********************"))
  }
  
  # --- 1.2. Define the different folder path
  # --- 1.2.1. Work in an existing directory
  in_path <- out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  
  # --- 1.2.2. Duplicate to a new directory
  if(!is.null(LOAD_FROM)){
    in_path <- paste0(project_wd, "/output/", LOAD_FROM)
    out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  }
  
  # --- 1.3. Load parameters
  load(paste0(in_path, "/CALL.RData"))
  
  # --- 2. Create new directory if needed
  # If we want to use an old list_bio without changing the directory
  if(!is.null(LOAD_FROM)){
    if(file.exists(out_path)==TRUE){
      stop("--- This new foldername is already used")
    } else {
      dir.create(out_path)
      sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/log.txt"), open = "wt"),
                           START = TRUE)
      message(paste(Sys.time(), "******************** START : run_init ********************"))
    } # if exists
  } # if LOAD_FROM
  
  # --- 3. Create species sub-directories
  # Named by their Worms ID, OTU ID, or named proportion or richness depending on
  # the data source and type. Contains a QUERY object with the corresponding
  # species selection to consider in each sub folder.
  if(DATA_TYPE == "proportions"){
    dir.create(paste0(out_path, "/proportions"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/proportions/QUERY.RData"))
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "omic"){
    dir.create(paste0(out_path, "/richness"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/richness/QUERY.RData"))
  } else {
    for(i in SP_SELECT){
      dir.create(paste0(out_path, "/", i))
      QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = i))
      save(QUERY, file = paste0(out_path, "/", i, "/QUERY.RData"))
    }
  } # if DATA_TYPE and SOURCE
  
  # --- 4. Extract environmental raster from NCDF
  # Provides a list per month, containing raster stack of all variables
  ENV_DATA <- list()
  
  # --- Start loop over month
  for(m in 1:12){
    # --- Compute over variables
    stack_month <- mclapply(X = 1:nrow(ENV_VAR), FUN = function(x){
      
      # --- 4.1. Open NC
      nc <- list.dirs(ENV_PATH, full.names = TRUE) %>% 
        .[grep(pattern = ENV_VAR$name[x], x = .)] %>% 
        list.files(full.names = TRUE) %>% 
        .[grep(pattern = paste0(str_pad(string = m, pad = "0", width = 2, side = "left"),"_"), x = .)] %>% 
        nc_open()
      
      # --- 4.2. Open the variable
      ncvar <- ncvar_get(nc, ENV_VAR$ncvar[x])
      
      # --- 4.3. Get the correct depth bounds
      if(length(dim(ncvar)) == 3){
        depth_id <- ncvar_get(nc, "depth_bnds") %>% t() %>% 
          as.data.frame() %>% 
          mutate(ID = row_number()) %>% 
          dplyr::filter(V1 >= CALL$SAMPLE_SELECT$MIN_DEPTH & V2 <= CALL$SAMPLE_SELECT$MAX_DEPTH) %>% 
          .$ID
      }
      
      # --- 4.4. Filter the right layer if there is a depth dimension
      if(length(dim(ncvar)) == 3){
        ncvar <- ncvar[,,depth_id] %>% apply(c(1,2), function(x)(x = mean(x, na.rm = TRUE)))
      }
      
      # --- 4.5. Get latitudes & longitudes
      lon <- nc$dim$lon$vals %>% as.numeric()
      lat <- nc$dim$lat$vals %>% as.numeric()
      
      # --- 4.6. Build the raster
      r <- raster(t(ncvar), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat)) %>% 
        flip(direction = 'y') 
    },
    mc.cores = min(nrow(ENV_VAR), MAX_CLUSTERS)) %>% stack() %>% synchroniseNA()
    
    # --- 4.7. Pretty names
    names(stack_month) <- ENV_VAR$name
    
    # --- 4.8. Assign to list
    ENV_DATA[[m]] <- stack_month
    message(paste(Sys.time(), "--- ENV. STACK : month", m, "done \t"))
  } # end month loop
  
  # --- 5. Check for biological data homogeneity ?
  # To be implemented ?
  
  # --- 6. Update CALL object
  # --- 6.1. Append CALL with ENV_DATA
  CALL[["ENV_DATA"]] <- ENV_DATA
  
  # --- 6.2. Append CALL with DATA_TYPE
  # Define it the same as DATA_SOURCE if left blank
  if(is.null(DATA_TYPE)){DATA_TYPE <- CALL$DATA_SOURCE}
  CALL[["DATA_TYPE"]] <- DATA_TYPE
  
  # --- 6.3. Append CALL with all other objects
  CALL[["SP_SELECT"]] <- SP_SELECT
  CALL[["FAST"]] <- FAST
  CALL[["ENV_VAR"]] <- ENV_VAR$name
  CALL[["ENV_PATH"]] <- ENV_PATH
  CALL[["METHOD_PA"]] <- METHOD_PA
  CALL[["NB_PA"]] <- NB_PA
  CALL[["PER_RANDOM"]] <- PER_RANDOM
  CALL[["DIST_PA"]] <- DIST_PA
  CALL[["BACKGROUND_FILTER"]] <- BACKGROUND_FILTER
  CALL[["OUTLIER"]] <- OUTLIER
  CALL[["UNIVARIATE"]] <- UNIVARIATE
  CALL[["ENV_COR"]] <- ENV_COR
  CALL[["NFOLD"]] <- NFOLD
  CALL[["FOLD_METHOD"]] <- FOLD_METHOD
  CALL[["MODEL_LIST"]] <- MODEL_LIST
  CALL[["LEVELS"]] <- LEVELS
  CALL[["ENSEMBLE"]] <- ENSEMBLE
  CALL[["N_BOOTSTRAP"]] <- N_BOOTSTRAP
  CALL[["CUT"]] <- CUT
  CALL[["PROJ_PATH"]] <- PROJ_PATH
  
  # --- 7. Wrap up and save
  # --- 7.1. Save file(s)
  save(CALL, file = paste0(out_path, "/CALL.RData"))
  
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
  # --- 7.3. List sub-directories to return
  # To be returned as a vector, further used to define parallel runs
  if(DATA_TYPE == "proportions"){
    parallel <- "proportions"
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "omic"){
    parallel <- "richness"
  } else {
    parallel <- CALL$SP_SELECT
  }
  return(parallel)
  
} # END FUNCTOIN

