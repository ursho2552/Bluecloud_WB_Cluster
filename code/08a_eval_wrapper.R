#' =============================================================================
#' @name eval_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data.
#' @param FOLDER_NAME name of the corresponding folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @param ENSEMBLE if TRUE, computes variable importance metrics for the ensemble 
#' model as well
#' @return the MODEL object updated with evaluation metric values (model performance
#' metric and variable importance metric)
#' @return outputs are saved in the MODEL.RData object

eval_wrapper <- function(FOLDER_NAME = NULL,
                         SUBFOLDER_NAME = NULL,
                         ENSEMBLE = TRUE){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - append file
  sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/log.txt"), open = "a"),
                       START = TRUE)
  message(paste(Sys.time(), "******************** START : eval_wrapper ********************"))
  
  # --- 1.2. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 1.3. Start PDF saving - variable importance
  pdf(paste0(project_wd,"/output/",FOLDER_NAME,"/",SUBFOLDER_NAME,"/04_variable_importance.pdf"))
  
  # --- 2. Redirection to PRESENCE model evaluation
  if(CALL$DATA_TYPE == "pres"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/08b_eval_pres.R"))
    
    # --- 2.2. Run function
    MODEL<- eval_pres(QUERY = QUERY,
                      MODEL = MODEL,
                      ENSEMBLE = ENSEMBLE)
  } # END if pres
  
  # --- 3. Redirection to CONTINUOUS model evaluation
  if(CALL$DATA_TYPE == "cont"){
    # --- 3.1. Load function
    source(file = paste0(project_wd, "/code/08c_eval_cont.R"))
    
    # --- 3.2. Run function
    MODEL <- eval_cont(QUERY = QUERY,
                       MODEL = MODEL,
                       ENSEMBLE = ENSEMBLE)
  } # END if pres
  
  # --- 4. Redirection to PROPORTION model evaluation
  if(CALL$DATA_TYPE == "proportions"){
    # --- 3.1. Load function
    source(file = paste0(project_wd, "/code/08d_eval_proportions.R"))
    
    # --- 3.2. Run function
    MODEL <- eval_proportions(CALL = CALL,
                              QUERY = QUERY,
                              MODEL = MODEL)
  } # END if pres
  
  # --- 5. Wrap up and save
  # --- 5.1. Save file(s)
  save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SUBFOLDER_NAME, "/MODEL.RData"))
  # --- 5.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  # --- 5.3. Stop PDF
  dev.off()
  
} # END FUNCTION
