#' =============================================================================
#' @name eval_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return the MODELS object updated with evaluation metric values
#' @return outputs are saved in the MODEL.RData object

eval_wrapper <- function(FOLDER_NAME = NULL,
                         SUBFOLDER_NAME = NULL){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : eval_wrapper ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 2. Redirection to PRESENCE model evaluation
  if(CALL$DATA_TYPE == "pres"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/08b_eval_pres.R"))
    
    # --- 2.2. Run function
    MODEL<- eval_pres(QUERY = QUERY,
                      MODEL = MODEL)
  } # END if pres
  
  # --- 3. Redirection to CONTINUOUS model evaluation
  if(CALL$DATA_TYPE == "cont"){
    # --- 3.1. Load function
    source(file = paste0(project_wd, "/code/08c_eval_cont.R"))
    
    # --- 3.2. Run function
    MODEL <- eval_cont(QUERY = QUERY,
                       MODEL = MODEL)
  } # END if pres
  
  # --- 4. Redirection to PROPORTION model evaluation
  # TO BE IMPLEMENTED
  
  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
} # END FUNCTION
