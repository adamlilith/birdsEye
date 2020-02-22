#' Ellipsoid from a Spatial* object, CRS, or proj4 string
#'
#' These functions extract the ellipsoid or datum from a spatial object, a CRS object, or from a proj4 string.
#' @param x SpatialPoints*, SpatialPolygons*, or SpatialLines* object, or an object of class \code{CRS} or \code{character} representing a proj4 string.
#' @return Character.
#' @examples
#' crs <- enmSdm::getCRS('albersNA')
#' ellipsoid(crs)
#' datum(crs)
#' @export
ellipsoid <- function(x) {

	cl <- class(x)
	if (cl %in% c('SpatialPoints', 'SpatialPointsDataFrame', 'SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpaitaLines', 'SpatialLinesDataFrame')) {
		x <- sp::proj4string(x)
	} else if (cl %in% c('CRS')) {
		x <- as.character(x)
	}
	
	split <- strsplit(x, split=' ')[[1]]
	whichOne <- grepl(split, pattern='+ellps=')
	thisOne <- split[whichOne]
	
	out <- strsplit(thisOne, split='=')[[1]][2]
	out

}

#' @describeIn ellipsoid Datum from a Spatial* object, CRS, or proj4 string
#' @export
datum <- function(x) {

	cl <- class(x)
	if (cl %in% c('SpatialPoints', 'SpatialPointsDataFrame', 'SpatialPolygons', 'SpatialPolygonsDataFrame', 'SpaitaLines', 'SpatialLinesDataFrame')) {
		x <- sp::proj4string(x)
	} else if (cl %in% c('CRS')) {
		x <- as.character(x)
	}
	
	split <- strsplit(x, split=' ')[[1]]
	whichOne <- grepl(split, pattern='+datum=')
	thisOne <- split[whichOne]
	
	out <- strsplit(thisOne, split='=')[[1]][2]
	out

}
