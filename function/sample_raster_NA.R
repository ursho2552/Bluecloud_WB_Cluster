
# --- HOMEMADE FUNCTIONS
# To transfer to a function file
# Doing a homemade function to sample from nearest cell if there is a NA due to coarse resolution
sample_raster_NA <- function(r, xy){
  apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
}