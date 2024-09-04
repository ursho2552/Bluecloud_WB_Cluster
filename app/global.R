# --- Global launch elements

message("load global libraries")
message(paste("Current directory:", getwd()))
message(paste(list.files(), collapse = " ; "))

library(shiny)
library(raster)
library(RNetCDF)
library(leaflet)
library(fields)
library(RColorBrewer)
library(shinythemes)

message("DONE - libraries")

# --- Useful functions

mako_pal <- function(n){
  tmp <- c("black", "#0A0A32", "#142C50", "#1F4D70", "#2A6D8F", "#368DAE", "#46ABBB", "#64C5BC", "#8DD6B8", "#B7E5B4","white")
  pal <- colorRampPalette(tmp)(n)
} # END FUNCTION