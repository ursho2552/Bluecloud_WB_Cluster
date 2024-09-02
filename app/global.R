message("load global libraries")
message(paste("Current directory:", getwd()))
message(paste(list.files(), collapse = " ; "))

library(shiny)
library(raster)
library(RNetCDF)
library(leaflet)
library(shinydashboard)
library(fields)
library(viridis)
library(RColorBrewer)
library(leaflet)


message("DONE - libraries")
