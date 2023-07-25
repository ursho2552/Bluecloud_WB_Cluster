#' =============================================================================
#' @name query_omic
#' @description extracts biological from the MGNIFY data downloaded in list_bio, 
#' according to a user provided list of species, time and depth range. 
#' The extracted data is formatted to be directly usable by the models available
#'  in this workbench.
#' @param FOLDER_NAME name of the corresponding folder
#' @param QUERY the QUERY.RData object containing the list of species
#'(depth range later - TO IMPLEMENT)
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth
#' @return Annotations: a data frame of taxonomic (maybe functional) annotation
#' for each target species.
#' @return output in a QUERY object

query_omic <- function(FOLDER_NAME = NULL,
                         QUERY = NULL){
  
  # --- 1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  SP_SELECT <- QUERY$SUBFOLDER_INFO$SP_SELECT
  
  # --- 2. Select OTU in the annotation table
  annotations <- CALL$LIST_BIO$annotations %>% 
    dplyr::filter(OTU %in% SP_SELECT)
  
  # --- 3. Propagate selection to the Y table
  # --- 3.1. Propagate selection on columns
  Y <- CALL$LIST_BIO$Y %>% 
    dplyr::select(all_of(SP_SELECT)) %>% 
    as.data.frame()
  # --- 3.2. Remove stations with no observations
  # We duplicate the first column to be able to apply() both in multi- univariate
  to_remove <- which(apply(cbind(Y, Y[,1]), 1, sum) == 0)
  Y <- Y[-to_remove,]
  
  # --- 4. Compute richness, presence abs, or prop
  # --- 4.1. Transform to proportion data
  if(CALL$DATA_TYPE == "proportions"){
    Y <- apply(Y, 1, function(x)(x = x/sum(x))) %>% 
      aperm(c(2,1)) %>% 
      as.data.frame()
  }

  # --- 4.2. Transform to presence data
  if(CALL$DATA_TYPE == "binary"){
    Y <- Y/Y %>% 
      as.data.frame()
    names(Y) <- "measurementvalue"
  }

  # --- 4.3. Transform to richness
  # In form of a Shannon index for now but we would compute more
  if(CALL$DATA_TYPE == "continuous"){
    Y <- apply(Y, 1, function(x)(x = vegan::diversity(x, "shannon"))) %>% 
      as.data.frame()
    names(Y) <- "measurementvalue"
  }

  # --- 5. Propagate selection to corresponding stations
  S <- CALL$LIST_BIO$S %>% 
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude))
  S <- S[-to_remove,]
  
  # --- 6. Save in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)
  
} # END FUNCTION