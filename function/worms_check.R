#' @description little utility function to check for worms ID validity, synonyms & childrens
#' @param ID an AphiA ID
#' @param MARINE_ONLY TRUE or FALSE, checks if the ID corresponds to a Marine taxa
#' @return a list for each ID, containing VALID, SYNONYM and CHILDREN vectors. 

worms_check <- function(ID, MARINE_ONLY = TRUE){
  # 1. Acceptance check - early stop
  init_check <- wm_record(id = ID)
  if(length(init_check) == 0){
    message(paste("WORMS: Species n°", ID, "is not registered in WoRMS - please make you to give relevant WoRMS ID or turn the WORMS_CHECK to FALSE if not adapted to your data \n"))
    return(NULL)
    } # end if
  
  # 2. Check for marine (optional)
  if(MARINE_ONLY == TRUE & init_check$isMarine != 1){
    message(paste("WORMS: Species n°", ID, "is not Marine - Discarded \n"))
    return(NULL)
  } # end if
  
  # 3. Get synonyms
  synonym_check <- wm_synonyms_(id = ID)
  if(length(synonym_check) == 0){
    VALID <- ID # security
    SYNONYM <- NULL
  } else {
    VALID <- synonym_check$valid_AphiaID %>% unique() %>% .[1]# security
    SYNONYM <- synonym_check$AphiaID %>% unique()
  }  # end if
  
  # 5. Return
  return(list(VALID = VALID, SYNONYM = SYNONYM))
  
} # END FUNCTION
