# birdsEye
This package is a complement to the popular `sp` package for R by Edzer Pebesma, Roger Bivand, and Virgilio Gomez-Rubio. Its contains a suite of functions for converting spatial objects into other spatial objects and for interpolating between spatial polygons.

Please note that `birdsEye` is in active development. Many of the tools work, but I have yet to do extensive testing. Please be cognizant of this as you use it.  I really appreciate issues brought to my attention as you find them.

You can install this package in R using these commands:

`install.packages('devtools') # if you haven't done this already`  
`library(devtools)`  
`devtools::install_github('adamlilith/omnibus')`  
`devtools::install_github('adamlilith/enmSdm')`  
`devtools::install_github('adamlilith/birdsEye')`  

NB: If for some reason these commands don't work, you can install the package(s) by downloading the latest zip/tar file from the `zipTarFiles` directory and installing the package(s) manually. If you do this, you will also have to install the `omnibus` and `enmSdm` packages, which are on GitHub also under my account (`adamlilith`).

## Spatial interpolation
* `interpPolysByBuffer` and `interpPolysByTween`: Interpolate between two spatial polygons

## Coordinate reference systems
* `makeCRS`: Return a custom proj4string (coordinate reference system string)

## Utilities
* `coordsToLine`: Convert coordinates to a spatial line
* `coordsToLines`: Convert a list of coordinates to spatial lines
* `coordsToPoly`: Convert coordinates to a spatial polygon
* `countSubGeoms`: Number of lines or polygons in a SpatialLines* or SpatialPolygons* object
* `datum`: Extract datum
* `ellipsoid`: Extract ellipsoid
* `geomToCoords`: Extract coordinates from a spatial object
* `subGeom`: Extract a sub-geometry (line, polygon) from a SpatialLines* or SpatialPolygons* object

Adam