#' =============================================================================
#' @name query_custom
#' @description extracts biological from the provided custom input, according to a
#' user provided list of species, time and depth range. The extracted data is
#' formatted to be directly usable by the models available in this workbench.
#' @param FOLDER_NAME name of the corresponding folder
#' @param QUERY the QUERY.RData object containing the list of species
#' (depth range later - TO IMPLEMENT)
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth
#' @return Annotations: a data frame of taxonomic (maybe functional) annotation
#' for each target species.
#' @return the output in a QUERY object

query_custom <- function(FOLDER_NAME = NULL,
                         QUERY = NULL){
  
  # --- 1. Parameter loading
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  
  # --- 2. Create Y target table
  # --- 2.1. For presence data
  if(CALL$DATA_TYPE == "binary"){
    Y <- CALL$LIST_BIO %>% 
      dplyr::select(measurementvalue) %>% 
      apply(2, function(x)(x/x)) %>% 
      data.frame()
  }
  
  # --- 2.2. For abundance data
  if(CALL$DATA_TYPE == "continuous"){
    Y <- CALL$LIST_BIO %>% 
      dplyr::select(measurementvalue) 
  }
  
  # --- 2.3. For proportion data
  if(CALL$DATA_TYPE == "proportions"){
    target_proportions <- CALL$LIST_BIO %>% 
      dplyr::select("decimallatitude","decimallongitude", "depth","year","measurementvalue","worms_id") %>% 
      pivot_wider(names_from = "worms_id", values_from = "measurementvalue", values_fn = mean)
    
    Y <- target_proportions %>% 
      dplyr::select(-decimallatitude, -decimallongitude, -depth, -year)
    Y[is.na(Y)] <- 0
    Y <- apply(Y, 1, function(x)(x = x/sum(x))) %>% aperm(c(2,1)) %>% as.data.frame()
  }
  
  # --- 3. Create S sample table
  S <- CALL$LIST_BIO %>% 
    dplyr::select(-any_of(c("measurementvalue", "worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    mutate(decimallatitude = as.numeric(decimallatitude),
           decimallongitude = as.numeric(decimallongitude),
           ID = row_number())
  
  # --- 4. Create annotation table
  annotations <- CALL$LIST_BIO %>% 
    dplyr::select(any_of(c("worms_id", "taxonrank", "scientificname", "nb_occ"))) %>% 
    distinct()
  
  # --- 5. Save in the QUERY object
  QUERY[["Y"]] <- Y
  QUERY[["S"]] <- S
  QUERY[["annotations"]] <- annotations
  
  return(QUERY)
  
} # END FUNCTION