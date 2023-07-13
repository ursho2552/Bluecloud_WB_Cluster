

regrid_env <- function(FOLDER_NAME, 
                       ENV_PATH){
  
  # --- 1. Initialize objects
  # --- 1.1. Information list
  r_info <- list(ENV_PATH = NULL, r = NULL, extent = NULL, res = NULL, na_nb = NULL)
  # --- 1.2. Rewrite trigger
  do_write <- FALSE
  
  # --- 2. Extract raster information
  # Here we check for extent, resolution and number of NA cells
  if(length(ENV_PATH > 1)){
    for(i in 1:length(ENV_PATH)){
      r <- stack(ENV_PATH[i])
      r_info$r[[i]] <- r
      r_info$ENV_PATH <- rbind(r_info$ENV_PATH, ENV_PATH[i])
      r_info$extent <- rbind(r_info$extent, as.vector(extent(r)))
      r_info$res <- rbind(r_info$res, as.vector(res(r)))
      r_info$na_nb <- rbind(r_info$na_nb, r[[1]][is.na(r[[1]])] %>% length())
    } # for
  } # if
  
  # --- 3. Homogenize
  # --- 3.1. Crop to the same extent
  if(r_info$extent %>% unique() %>% nrow() > 1){
    message("The input environmental variables have different extent \n
          --- Automatic cropping")
    
    min_ext <- apply(r_info$extent, 2, min)
    r_info$r <- lapply(r_info$r, function(x)(x = raster::crop(x, extent(min_ext))))
    do_write <- TRUE
  } # if
  
  # --- 3.2. Aggregate to the same resolution
  if(r_info$res %>% unique() %>% nrow() > 1){
    message("The input environmental variables have different resolutions \n
          --- Automatic re-gridding")
    
    max_res <- apply(r_info$res, 2, max)
    r_info$r <- lapply(r_info$r, function(x)(x = aggregate(x, fact = c(max_res/res(x), 1))))
    do_write <- TRUE
  } # if
  
  # --- 3.3. Synchronize NA between layers
  if(r_info$na_nb %>% unique() %>% length() > 1){
    message("The input environmental variables have different number of NA \n
          --- Automatic synchronizing")
    
    # --- 3.3.1. Build the NA mask
    sync_r <- stack(ENV_PATH)
    sync_r <- synchroniseNA(sync_r)[[1]]
    
    # --- 3.3.2. Perform the global sync
    r_info$r <- lapply(r_info$r, function(x)(x = synchroniseNA(stack(sync_r, x))[[-1]]))
    do_write <- TRUE
  }
  
  # --- 4. Wrap up
  if(do_write == TRUE){
    # --- 4.1. Write in the output object
    features_regrid <- stack(unlist(r_info$r))
    writeRaster(features_regrid, paste0(project_wd, "/output/", FOLDER_NAME,"/features_regrid"))
    
    # --- 4.2. Update the ENV_VAR object
    ENV_PATH <- paste0(project_wd, "/output/", FOLDER_NAME,"/features_regrid")
    
  }
  
  return(ENV_PATH)
} # END FUNCTION

