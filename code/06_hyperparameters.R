#' =============================================================================
#' @name hyperparameter
#' @description this function creates an hyperparameter grid for each available
#' algorithm in the pipeline and returns those selected by the user
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @param MODEL_LIST list of algorithms from which to compute hyperparameter
#' selection
#' @param LEVELS maximum number of parameter values to test in each of the grids
#' @return a hyperparameter list of the chosen models, including spec and grid
#' per model
#' @return the general HP list is saved in the CALL.RData object

hyperparameter <- function(FOLDER_NAME = NULL,
                           MODEL_LIST = c("GLM","GAM","RF","MLP"),
                           LEVELS = 3){ # TO DO : double check model names

  # --- 1. Initialize the  object
  # --- 1.1. Output object
  HP <- list()
  # --- 1.2. Load CALL object
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 2. Fill up the object with a set of possible grids
  # --- 2.1. RANDOM FOREST 
  # --- 2.1.1. Define the model specifications
  HP$RF$model_spec <- rand_forest(mode = "regression",
                                  engine = "randomForest",
                                  mtry = tune(),
                                  trees = tune(),
                                  min_n = tune())
  
  # --- 2.1.2. Define the grid according to built in functions
  HP$RF$model_grid <- grid_regular(mtry(range = c(1, length(CALL$ENV_VAR))),
                                   trees(),
                                   min_n(),
                                   levels = LEVELS)
  
  # --- 2.2. GENERALIZED ADDITIVE MODELS 
  # --- 2.2.1. Define the model specifications
  if(CALL$DATA_TYPE == "pres"){
    HP$GAM$model_spec <- gen_additive_mod(mode = "regression",
                                          adjust_deg_free = tune(),
                                          select_features = tune()) %>% 
      set_engine("mgcv",
                 family = stats::binomial(link = "logit")) %>% 
      translate()
  } else {
    HP$GAM$model_spec <- gen_additive_mod(mode = "regression",
                                          engine = "mgcv",
                                          adjust_deg_free = tune(),
                                          select_features = tune())
  }
  
  # --- 2.2.2. Define the grid according to built in functions
  HP$GAM$model_grid <- grid_regular(adjust_deg_free(),
                                    select_features(),
                                    levels = LEVELS)
  
  # --- 2.3. GENERALIZED LINEAR MODELS
  # --- 2.3.1. Define the model specifications
  if(CALL$DATA_TYPE == "pres"){
    HP$GLM$model_spec <- linear_reg(mode = "regression") %>% 
      set_engine("glm", 
                 family = stats::binomial(link = "logit")) %>% 
      translate()
  } else {
    HP$GLM$model_spec <- linear_reg(mode = "regression",
                                    engine = "glm")
  }

  # --- 2.3.2. Define the grid according to built in functions
  # No parameters to tune for engine "glm"
  
  # --- 2.4. SINGLE LAYER NEURAL NETWORK
  # --- 2.4.1. Define the model specifications
  HP$MLP$model_spec <- mlp(mode = "regression",
                           engine = "nnet",
                           hidden_units = tune(),
                           penalty = tune(),
                           dropout = NULL,
                           epochs = 5,
                           activation = NULL,
                           learn_rate = 1e-2)
  
  # --- 2.4.2. Define the grid according to built in functions
  HP$MLP$model_grid <- grid_regular(hidden_units(),
                                    penalty(),
                                    levels = LEVELS)
  
  # --- 3. Hyper parameter selection
  # According to the specified model list
  HP <- HP[MODEL_LIST]
  HP$CALL$MODEL_LIST <- MODEL_LIST
  
  # --- 4. Append CALL and save
  CALL[["HP"]] <- HP
  save(CALL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
} # END FUNCTION
