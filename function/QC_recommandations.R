#' @name qc_recommantations
#' @description function to extract the recommendations according to a set of
#' predefined quality checks of the algorithms.
#' @param MODEL the model object returned by the projection step
#' @param ENSEMBLE compute the QC for the ensemble
#' @return a recommendation table with the quality checks and associated text

qc_recommandations <- function(MODEL,
                               ENSEMBLE = FALSE,
                               RECOMMANDATIONS_DF = data.frame(FIT = c(0,0,1,0,1,0,1,1),
                                                               VIP = c(0,1,0,0,0,1,1,1),
                                                               DEV = c(0,0,0,1,1,1,0,1),
                                                               Recommandation = c("Do not use. Work on predictors",
                                                                                  "Do not use. Work on models",
                                                                                  "Do not use. Work on predictors",
                                                                                  "Do not use. Work on predictors",
                                                                                  "Promising. Work on predictors",
                                                                                  "Promising. Work on models and data",
                                                                                  "Promising. Work on models and data",
                                                                                  "Write a proposal !"),
                                                               COL = c("#B64A60","#B64A60","#B64A60","#B64A60","#ffc800","#ffc800","#ffc800","#1F867B"))){
  
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
  qc_names <- c("FIT","VIP","DEV")
  
  # --- 2. Build quality check matrix
  # --- 2.1. Create
  qc_matrix <- matrix(NA, nrow = length(m_names), ncol = length(qc_names), 
                      dimnames = list(m_names, qc_names))
  
  # --- 2.2. Fill up
  for(m in m_names){
    tmp <- MODEL[[m]][["eval"]] %>% unlist()
    id <- names(tmp) %>% grep(pattern = "CBI|R2|CUM_VIP|NSD")
    tmp <- tmp[id]
    while(length(tmp) < length(qc_names)){tmp <- c(tmp, NA)}
    qc_matrix[m,] <- tmp
  }
  
  # --- 2.3. Transform to 0 and 1 according to predefined criterias
  qc_matrix_01 <- qc_matrix
  for(m in m_names){
    if(qc_matrix[m,1] >= 0.5 & !is.na(qc_matrix[m,1])){qc_matrix_01[m,1] <- 1} else {qc_matrix_01[m,1] <- 0}
    if(qc_matrix[m,2] >= 50 & !is.na(qc_matrix[m,2])){qc_matrix_01[m,2] <- 1} else {qc_matrix_01[m,2] <- 0}
    if(qc_matrix[m,3] <= 0.5 & !is.na(qc_matrix[m,3])){qc_matrix_01[m,3] <- 1} else {qc_matrix_01[m,3] <- 0}
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


