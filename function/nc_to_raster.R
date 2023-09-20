#' =============================================================================
#' @name nc_to_raster
#' @description extracts the environmental data from a .nc files and store
#' them in one raster per month considered... easier to
#' move the data around in the pipeline.
#' @param MONTH month to extract the data from
#' @param NC environmental .nc files to extract the data from
#' @param MIN_DEPTH minimum depth range
#' @param MAX_DEPTH maximum depth range

nc_to_raster <- function(MONTH, NC, MIN_DEPTH, MAX_DEPTH){
  # --- 1. Read .nc informations
  # --- 1.1. Open file
  nc <- nc_open(NC)
  
  # --- 1.2. Extract the variable raw data
  # The variable of interest has to be the first variable in the .nc file
  varname <- nc$var %>% names()
  ncvar <- ncvar_get(nc, varname[1])
  
  # --- 1.3. Extract dimensions and match their ordering with data
  nc_dimsize0 <- lapply(nc$dim, FUN = function(x)(x = length(x$vals)))
  id <- match(dim(ncvar), unlist(nc_dimsize0))
  nc_dimsize0 <- nc_dimsize0[id]
  
  # --- 1.4. Get latitudes & longitudes
  lon <- nc$dim$lon$vals %>% as.numeric()
  lat <- nc$dim$lat$vals %>% as.numeric()
  res <- abs(lon[1])-abs(lon[2])
  
  # --- 2. Transpose file and dimension to the desired order
  id <- match(c("lon","lat","time","depth"), names(nc_dimsize0)) %>% .[!is.na(.)]
  ncvar <- aperm(ncvar, id)
  nc_dimsize <- nc_dimsize0[id]
  
  # --- 3. Extract time first; if exist
  # Because we order the dim before, time is always in 3rd position
  if(!is.na(names(nc_dimsize[3])) & names(nc_dimsize)[3] == "time"){
    # --- 3.1. By default we take the first layer
    t_id <- 1
    if(nc_dimsize[[3]] > 1){
      # --- 3.2. If there is several, we take the corresponding month
      t_id <- MONTH
      # --- 3.3. Now we extract the dimension
      if(length(nc_dimsize) == 3){ncvar <- ncvar[,,t_id]
      } else {ncvar <- ncvar[,,t_id,]}
    } # if length time
    # --- 3.4. And we delete the from dimension list
    nc_dimsize <- nc_dimsize[-3]
  } # if time exist
  
  # --- 4. Extract depth; if exist
  # If time existed, it is now extracted; so depth (if exist) is 3rd dimension
  if(!is.na(names(nc_dimsize[3])) & names(nc_dimsize)[3] == "depth"){
    # --- 4.1. By default we take the first layer
    d_id <- 1
    if(nc_dimsize[[3]] > 1){
      # --- 4.2. If there is several, we take the corresponding range
      depth_bnds <- abs(depth[-length(depth)]-depth[-1])
      depth <- depth[-length(depth)]
      depth_id <- which(depth >= CALL$SAMPLE_SELECT$MIN_DEPTH & depth <= CALL$SAMPLE_SELECT$MAX_DEPTH)
      # --- 4.3. Now we extract with weighted mean across depth bounds
      ncvar <- ncvar[,,depth_id] %>% 
        apply(c(1,2), function(x)(x = sum(x*depth_bnds[depth_id]/sum(depth_bnds[depth_id]))))
    } # if length depth
  } # if depth exist
  
  # --- 5. Build the raster
  r <- raster(t(ncvar), xmn = min(lon-0.5*res), xmx = max(lon+0.5*res), 
              ymn = min(lat-0.5*res), ymx = max(lat+0.5*res))
  # --- 6. Reverse if latitudes and/or longitude are reversed
  if(lat[1]-lat[2] < 0){r <- flip(r, direction = "y")}
  if(lon[1]-lon[2] > 0){r <- flip(r, direction = "x")}
  return(list(r, nc_dimsize0))
} # End function
