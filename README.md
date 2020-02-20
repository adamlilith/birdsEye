# birdsEye
This package is a complement to the popular `dismo` package for R by Robert Hijmans. Its contains a suite of efficiency functions for preparing data, training and evaluating species distribution models and ecological niche models, and comparing ecological niches.

You can install this package in R using these commands:

`install.packages('devtools') # if you haven't done this already`  
`library(devtools)`  
`devtools::install_github('adamlilith/omnibus')`  
`devtools::install_github('adamlilith/enmSdm')`  
`devtools::install_github('adamlilith/birdsEye')`  

NB: If for some reason these commands don't work, you can install the package(s) by downloading the latest zip/tar file from the `zipTarFiles` directory and installing the package(s) manually. If you do this, you will also have to install the `omnibus` and `enmSdm` packages, which are on GitHub also under my account (`adamlilith`).

## Spatial interpolation
* `interpPolysByBuffer` and `interpPolysByDist`: Interpolate between two spatial polygons

## Utilities
* `coordsToLine`: Convert coordinates to a spatial line
* `coordsToLines`: Convert a list of coordinates to spatial lines
* `coordsToPoly`: Convert coordinates to a spatial polygon
* `countSubGeom`: Number of lines or polygons in a SpatialLines* or SpatialPolygons* object
* `geomToCoords`: Extract coordinates from a spatial object
* `subGeomFromGeom`: Extract a sub-geometry (line, polygon) from a SpatialLines* or SpatialPolygons* object

Adam