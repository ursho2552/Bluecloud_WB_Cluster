#' @name qc_recommantations
#' @description function to extract the recommendations according to a set of
#' predefined quality checks of the algorithms.
#' @param MODEL the model object returned by the projection step
#' @param DATA_TYPE the data type from the CALL object, 
#' as recommendations differ from one to another
#' @param ENSEMBLE compute the QC for the ensemble
#' @return a recommendation table with the quality checks and associated text

qc_recommandations <- function(QUERY, MODEL,
                               DATA_TYPE,
                               ENSEMBLE = FALSE,
                               RECOMMANDATIONS_DF = data.frame(PRE_VIP = c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1),
                                                               FIT = c(0,0,1,0,1,0,1,1,0,0,1,0,1,0,1,1),
                                                               CUM_VIP = c(0,1,0,0,0,1,1,1,0,1,0,0,0,1,1,1),
                                                               DEV = c(0,0,0,1,1,1,0,1,0,0,0,1,1,1,0,1),
                                                               Recommandation = c("Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Start by working on predictors",
                                                                                  "Do not use. Algorithms are not well fitted",
                                                                                  "Do not use. Predictors are not meaningful",
                                                                                  "Do not use. Predictors are not meaningful",
                                                                                  "Promising but predictors are not meaningful",
                                                                                  "Promising but algorithms are not fitting",
                                                                                  "Promising but projection uncertainty is high",
                                                                                  "Satisfying for proposal writing"),
                                                               COL = c("#B64A60","#B64A60","#B64A60","#B64A60","#B64A60","#B64A60","#B64A60","#B64A60",
                                                                       "#B64A60","#B64A60","#B64A60","#B64A60","#ffc800","#ffc800","#ffc800","#1F867B"))){
  
  # --- 1. Extract labels
  # --- 1.1. Algorithm names
  if(ENSEMBLE == TRUE){
    m_names <- names(MODEL) %>% 
      .[. ==  "ENSEMBLE"]
  } else {
    m_names <- names(MODEL) %>% 
      .[. != "ENSEMBLE" & . !=  "MODEL_LIST"]
  }

  # --- 1.2. Quality checks names
  qc_names <- c("PRE_VIP","FIT","CUM_VIP","DEV")
  
  # --- 2. Build quality check matrix
  # --- 2.1. Create
  qc_matrix <- matrix(NA, nrow = length(m_names), ncol = length(qc_names), 
                      dimnames = list(m_names, qc_names))
  
  # --- 2.2. Fill up
  for(m in m_names){
    tmp <- MODEL[[m]][["eval"]] %>% unlist()
    id <- names(tmp) %>% grep(pattern = "CBI|R2|CUM_VIP|NSD") # get the right model QC
    tmp <- tmp[id]
    qc_matrix[m,1] <- QUERY$eval$PRE_VIP # fill the query QC
    while(length(tmp) < 3){tmp <- c(tmp, NA)} # to add missing QC due to early discard
    qc_matrix[m,2:4] <- tmp # fill the model QC
  }
  
  # --- 2.3. Transform to 0 and 1 according to predefined criterias
  qc_matrix_01 <- qc_matrix
  for(m in m_names){
    # --- 2.3.1. A priori predictor importance
    if(qc_matrix[m,1] >= 0.05 & !is.na(qc_matrix[m,1])){qc_matrix_01[m,1] <- 1} else {qc_matrix_01[m,1] <- 0}
    # --- 2.3.2. Perdictive performance
    if(DATA_TYPE == "binary"){
      if(qc_matrix[m,2] >= 0.5 & !is.na(qc_matrix[m,2])){qc_matrix_01[m,2] <- 1} else {qc_matrix_01[m,2] <- 0}
    } else {
      if(qc_matrix[m,2] >= 0.25 & !is.na(qc_matrix[m,2])){qc_matrix_01[m,2] <- 1} else {qc_matrix_01[m,2] <- 0}
    } # if R2 or CBI - we are less stringent from R2
    # --- 2.3.3. Cumulative variable importance
    if(qc_matrix[m,3] >= 50 & !is.na(qc_matrix[m,3])){qc_matrix_01[m,3] <- 1} else {qc_matrix_01[m,3] <- 0}
    # --- 2.3.4. Projection uncertainty
    if(qc_matrix[m,4] <= 0.5 & !is.na(qc_matrix[m,4])){qc_matrix_01[m,4] <- 1} else {qc_matrix_01[m,4] <- 0}
  }
  qc_matrix_01 <- as.data.frame(qc_matrix_01) %>% 
    mutate(ID = row_number())
  
  # --- 3. Build output matrix
  # rec <- as.data.frame(qc_matrix_01) %>% 
  #   left_join(RECOMMANDATIONS_DF)
  rec <- merge(x = qc_matrix_01, y = RECOMMANDATIONS_DF, sort = FALSE) %>% 
    arrange(ID) %>% 
    dplyr::select(-ID)
  rownames(rec) <- m_names
  
  # --- 4. Wrap up and save
  return(rec)
} # END FUNCTION


