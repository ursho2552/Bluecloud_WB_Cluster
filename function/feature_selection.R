#' @description little utility function to help decision making with features selection
#' Computes both the Spearman-rank correlation and Mutual Information criteria relative
#' to a NULL feature
#' @param V1 the target, a numeric vector
#' @param V2 the feature to test, a numeric vector
#' @param BINS number of bins to descritize the variables on, herited from infotheo::discretize
#' @param METHOD the method to compute entropy, herited from infotheo::entropy
#' @return a dataframe of two column, containing the Mutual Information and Spearman values

feature_selection <- function(V1, V2, BINS = ceiling(length(V1)^(1/3)), METHOD = "emp") {
  # --- 1. Discretize target and feature
  normalized_V1 <- discretize(V1, nbins = BINS, disc = "equalwidth")
  normalized_V2 <- discretize(V2, nbins = BINS, disc = "equalwidth")
  
  # --- 2. NULL target
  # null <- rep(mean(V1), length(V1)) # create a NULL at the average of the observations
  null <- sample(x = 100, size = length(V1), replace = TRUE)
  null_disc <- discretize(null, nbins = BINS, disc = "equalwidth")
  
  # --- 3. Compute Mutual information
  # --- 3.1. Raw MI // normalize by target entropy
  mi <- mutinformation(normalized_V1, normalized_V2, method = METHOD)
  norm <- entropy(normalized_V1,  method = METHOD)
  
  # --- 3.2. Same operation with the NULL target
  mi_null <- mutinformation(normalized_V1, null_disc, method = METHOD)
  
  # --- 3.3. Difference to the NULL & normalize
  mi <- (mi-mi_null)/norm
  
  # --- 4. Compute a spearman correlation
  # --- 4.1. Raw spearman
  s <- cor(V1, V2, method = "spearman") %>% abs()
  
  # --- 4.2. Same treatment for the null
  s_null <- cor(V1, null_disc, method = "spearman") %>% abs()
  
  # --- 4.3. Difference to the null
  s <- s-s_null
  
  # --- 5. Wrap up
  out <- data.frame(mutual_information = mi, spearman = s)
  
  return(out)
} # end function
