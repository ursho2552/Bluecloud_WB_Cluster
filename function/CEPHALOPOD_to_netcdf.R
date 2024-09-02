#' =============================================================================
#' @name CEPHALOPOD_to_netcdf
#' @title Write mean and standard deviation values to NetCDF
#' @description
#'   This function saves geographical and time-series data in NetCDF format,
#'   including mean and standard deviation values for multiple features.
#' @param r0   Baseline raster
#' @param a_m  A three-dimensional array (geo_coordinates x time x feature) containing mean values.
#' @param a_sd A three-dimensional array (geo_coordinates x time x feature) containing standard deviation values.
#' @param output_file Name of the NetCDF file to be created.
#' @return
#'   A NetCDF file containing mean and standard deviation values of features at geographical and time scales,
#'   saved in the specified output folder.
#' Saved in FOLDERNAME.

CEPHALOPOD_to_netcdf <- function(r0, a_m, a_sd, names, reco_list, output_file = "output_abundance_test.nc") {
  library(RNetCDF)
  library(raster)
  
  MONTH <- as.list(1:13)
  nlon <- ncol(r0)  # Assuming r0 is defined elsewhere
  nlat <- nrow(r0)  # Assuming r0 is defined elsewhere
  ntime <- length(MONTH)
  nfeatures <- dim(a_m)[3]  # Number of variables (formerly species)
  names <- as.character(names)
  print(names)
  lon_values <- seq(-180, 180, length.out = nlon)
  lat_values <- seq(-90, 90, length.out = nlat)
  
  # Create NetCDF file
  nc_file <- create.nc(output_file)
  
  # Define dimensions
  dim.def.nc(nc_file, dimname = "longitude", nlon)
  dim.def.nc(nc_file, dimname = "latitude", nlat)
  dim.def.nc(nc_file, dimname = "time", ntime)  # Not unlimited
  dim.def.nc(nc_file, dimname = "feature", nfeatures)
  
  # Define variables with float precision
  var.def.nc(nc_file, varname = "longitude", vartype = "NC_FLOAT", dimensions = "longitude")
  var.def.nc(nc_file, varname = "latitude", vartype = "NC_FLOAT", dimensions = "latitude")
  var.def.nc(nc_file, varname = "time", vartype = "NC_FLOAT", dimensions = "time")
  var.def.nc(nc_file, varname = "mean_values", vartype = "NC_FLOAT", dimensions = c("longitude", "latitude", "time", "feature"))
  var.def.nc(nc_file, varname = "sd_values", vartype = "NC_FLOAT", dimensions = c("longitude", "latitude", "time", "feature"))
  
  # Add attributes
  att.put.nc(nc_file, variable = "longitude", name = "units", type = "NC_CHAR", value = "degrees_east")
  att.put.nc(nc_file, variable = "latitude", name = "units", type = "NC_CHAR", value = "degrees_north")
  att.put.nc(nc_file, variable = "mean_values", name = "_FillValue", type = "NC_FLOAT", value = -9999)
  att.put.nc(nc_file, variable = "mean_values", name = "long_name", type = "NC_CHAR", value = "Mean Values")
  att.put.nc(nc_file, variable = "mean_values", name = "reco", type = "NC_CHAR", value = reco_list)
  att.put.nc(nc_file, variable = "mean_values", name = "feature_names", type = "NC_CHAR", value = paste(names, collapse = ","))
  att.put.nc(nc_file, variable = "sd_values", name = "_FillValue", type = "NC_FLOAT", value = -9999)
  att.put.nc(nc_file, variable = "sd_values", name = "long_name", type = "NC_CHAR", value = "Standard Deviation Values")
  att.put.nc(nc_file, variable = "time", name = "units", type = "NC_CHAR", value = "months since 1970-01-01")
  
  var.put.nc(nc_file, variable = "longitude", data = lon_values)
  var.put.nc(nc_file, variable = "latitude", data = lat_values)
  var.put.nc(nc_file, variable = "time", data = 1:ntime)
  
  # Write data to NetCDF
  for (i in 1:nfeatures) {
    for (m in 1:ntime) {
      tmp_m <- a_m[, m, i]
      tmp_sd <- a_sd[, m, i]
      plot_scale <- quantile(a_m, 0.95, na.rm = TRUE)
      
      # Prepare rasters
      r_m <- setValues(r0, tmp_m)
      r_m[r_m > plot_scale] <- plot_scale
      r_sd <- setValues(r0, tmp_sd)
      
      # Get matrix values from rasters
      rm_values <- getValues(r_m)
      rsd_values <- getValues(r_sd)
      
      # Reshape vector into a matrix with dimensions matching r_m
      mean_values <- matrix(rm_values, nrow = nlat, ncol = nlon, byrow = FALSE)
      sd_values <- matrix(rsd_values, nrow = nlat, ncol = nlon, byrow = FALSE)
      
      # Write values to NetCDF variables
      var.put.nc(nc_file, variable = "mean_values", data = mean_values, start = c(1, 1, m, i), count = c(nlon, nlat, 1, 1))
      var.put.nc(nc_file, variable = "sd_values", data = sd_values, start = c(1, 1, m, i), count = c(nlon, nlat, 1, 1))
    }
  }
  
  # Close the NetCDF file
  close.nc(nc_file)
}
