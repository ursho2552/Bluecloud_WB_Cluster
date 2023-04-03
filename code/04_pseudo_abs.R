#' =============================================================================
#' @name pseudo_abs
#' @description computes pseudo absences added to the X and Y matrices from
#' query_env and bio according to a convexhull or geographical method.
#' @param QUERY query object resulting from query_bio and query_env
#' @param NB_PA number of pseudo-absences to generate
#' @param METHOD_PA method of pseudo-absence, either "sre" or "disk"
#' @return X updated with the pseudo-absence values (= 0)
#' @return Y updated with the environmental values corresponding
#' 

pseudo_abs <- function(QUERY = query,
                       NB_PA = nrow(QUERY$S),
                       METHOD_PA = "sre"){
  
  # =============================== BY CONVEX HULL =============================
  # ----- FAR TOO HEAVY AND LONG
  # library(geometry)
  # build_chull <- QUERY$X %>% distinct() 
  # # Attention grosse approximation car en vrai, il faudrait prendre les données env. mensuelles
  # chull <- convhulln(build_chull, options = 'Qbk:0Bk:O', output.options = TRUE)
  # 
  # features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>% 
  #   readAll() %>% 
  #   rasterToPoints()
  # 
  # check_chull <- inhulln(chull, features[,-c(1,2)])
  # summary(check_chull)
  
  # ================= BY BIOCLIMATIQUE ENVELOPPE KIND OF =======================
  # --- Calculates the environmental presence boundaries for each variable
  # --- Selects presence outside the n-th quantile for each boundary
  # --- WORKS BUT SRE IS TAKING ALL AVAILABLE SPACE WITH THAT MUCH VARIABLES...
  
  # # Define environmental boundaries
  # sre <- QUERY$X %>% 
  #   apply(2, function(x) (x = quantile(x, probs = c(0.25, 0.75))))
  # # Filter environmental background
  # features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>%
  #   readAll() %>%
  #   rasterToPoints() %>% 
  #   as.data.frame()
  # for(c in colnames(sre)){
  #   tmp <- which(features[,c] < sre[1,c] | features[,c] > sre[2,c])
  #   features <- features[tmp,]
  # }  

  # ============================== GEOGRAPHICAL DISK ===========================
  # --- Open environmental data
  features <- stack(paste0(project_wd, "/data/features_mean_from_monthly")) %>%
    readAll()
  
  # --- Create base raster
  r <- features[[1]]
  r[!is.na(r)] <- 0
  
  # --- Extract presence points
  xy <- QUERY$S %>% 
    dplyr::select(decimallongitude, decimallatitude)
  val <- data.frame(val = rep(1, nrow(xy)))
  
  # --- Calculate distance to presences in raster
  background <- rasterize(xy, r, update=TRUE)
  background[background < 1] <- NA
  background <- distance(background)
  
  background <- synchroniseNA(stack(background, r))[[1]] %>% 
    rasterToPoints() %>% 
    as.data.frame() %>% 
    dplyr::filter(layer > 100e3 & !is.na(layer))
  
  # --- Sample background data
  tmp <- sample(x = 1:nrow(background), size = nrow(QUERY$S))
  xy <- background[tmp,1:2]
  
  # --- Append the query
  X <- raster::extract(features, xy) %>% 
    as.data.frame() %>% 
    rbind(QUERY$X)
  
  Y <- data.frame(measurementvalue = c(rep(0, nrow(xy)), 
                                       rep(1, nrow(QUERY$Y))))
  
  S <- data.frame(decimallongitude = xy$x,
                  decimallatitude = xy$y,
                  measurementtype = "Pseudo-absence") %>% 
    bind_rows(QUERY$S) %>% 
    dplyr::select(colnames(QUERY$S))

  return(list(Y = Y,
              X = X,
              S = QUERY$S,
              CALL = QUERY$CALL))
  
} # END FUNCTION
