#' =============================================================================
#' @name hyperparameter
#' @description this function creates an hyperparameter grid for each available
#' algorithm in the pipeline and returns those selected by the user
#' @param QUERY the query object from previous function
#' @param MODEL_LIST list of algorithms from which to compute hyperparameter
#' selection
#' @return a hyperparameter list of the chosen models, including spec and grid
#' per model

hyperparameter <- function(QUERY = query,
                           MODEL_LIST = c("GLM","GAM","RF","MLP")){ # TO DO : double check model names

  # --- 1. Generate the output object
  HP <- list()
  
  # --- 2. Fill up the object with a set of possible grids
  # --- 2.1. RANDOM FOREST 
  # --- 2.1.1. Define the model specifications
  HP$RF$model_spec <- rand_forest(mode = "regression",
                                  engine = "randomForest",
                                  mtry = tune(),
                                  trees = tune(),
                                  min_n = tune())
  
  # --- 2.1.2. Define the grid according to built in functions
  HP$RF$model_grid <- grid_regular(mtry(range = c(1, length(QUERY$CALL$ENV_VAR))),
                                   trees(),
                                   min_n(),
                                   levels = 5)
  
  # --- 2.2. GENERALIZED ADDITIVE MODELS 
  # --- 2.2.1. Define the model specifications
  if(QUERY$CALL$DATA_TYPE == "pres"){
    HP$GAM$model_spec <- gen_additive_mod(mode = "regression") %>% 
      set_engine("mgcv",
                 family = stats::binomial(link = "logit"),
                 adjust_deg_free = tune(),
                 select_features = tune()) %>% 
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
                                    levels = 5)
  
  # --- 2.3. GENERALIZED LINEAR MODELS
  # --- 2.3.1. Define the model specifications
  if(QUERY$CALL$DATA_TYPE == "pres"){
    HP$GLM$model_spec <- linear_reg(mode = "regression") %>% 
      set_engine("glm", 
                 family = stats::binomial(link = "logit"),
                 penalty = tune(),
                 mixture = tune()) %>% 
      translate()
  } else {
    HP$GLM$model_spec <- linear_reg(mode = "regression",
                                    engine = "glm",
                                    penalty = tune(),
                                    mixture = tune())
  }

  # --- 2.3.2. Define the grid according to built in functions
  HP$GLM$model_grid <- grid_regular(penalty(),
                                    mixture(),
                                    levels = 5)
  
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
                                    levels = 5)
  
  # --- 3. Hyper parameter selection
  # According to the specified model list
  HP <- HP[MODEL_LIST]
  HP$CALL$MODEL_LIST <- MODEL_LIST
  
  return(HP)
  
} # END FUNCTION