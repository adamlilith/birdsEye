#' Convert coordinates to spatial polygon(s) or spatial line(s)
#'
#' This function converts coordinates to a SpatialPolygon object. To close the polygon the last set of coordinates must be equal to the first. If they are not, then the first set is added to the end. Note that it is possible to create an improper polygon (i.e., that crosses itself).
#' @param x Either a two-column matrix or data frame with coordinates (first column is longitude, second is latitude), or a list of matrices or data frames, or an object from which coordinates can be extracted using \code{link[sp]{coordinates}}.
#' @param crs A character string or an object of class \code{CRS} which species the coordinate reference system (proj4string) for the polygon. If the object in \code{x} is a spatial object the CRS will be obtained from it unless \code{crs} is specified.
#' @param closePoly Logical, if \code{TRUE} (default), then close the polygon.
#' @param id Integer or character, ID value for the object(s).
#' @seealso \code{\link[sp]{SpatialPolygons}}
#' @return A spatial polygon object.
#' @examples
#' x <- c(-116.57, -99.13, -62.95, -65.21, -90.08)
#' y <- c(62.09, 75.50, 73.29, 68.34, 66.43)
#' xy <- cbind(x, y)
#' poly <- coordsToPoly(xy, getCRS('wgs84'))
#' line <- coordsToLine(xy, getCRS('wgs84'))
#' par(mfrow=c(1, 2))
#' plot(poly)
#' plot(line)
#' @export
coordsToPoly <- function(x, crs, closePoly = TRUE, id = 1) {

	if (class(x) %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
	
		coords <- sp::coordinates(x)
		crs <- sp::proj4string(x)
	
	} else {
	
		coords <- as.matrix(x)
	
	}

	if (class(crs) != 'CRS') crs <- sp::CRS(crs)
	
	# close the polygon
	if (closePoly) {
		n <- nrow(coords)
		if (coords[1, 1] != coords[n, 1] | coords[1, 2] != coords[n, 2]) {
			coords <- rbind(coords, coords[1, ])
		}
		
	}
	
	rownames(coords) <- NULL
	
	poly <- sp::Polygon(coords)
	poly <- sp::Polygons(list(poly), id)
	poly <- sp::SpatialPolygons(list(poly), proj4string=crs)
	poly

}

#' @describeIn coordsToPoly Convert coordinates to a spatial line
#' @export
coordsToLine <- function(x, crs, id = 1) {

	if (class(x) %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
	
		coords <- sp::coordinates(x)
		crs <- sp::proj4string(x)
	
	} else {
		coords <- as.matrix(x)
	}

	if (class(crs) != 'CRS') crs <- sp::CRS(crs)

	rownames(coords) <- NULL
	
	sen <- sp::Line(coords)
	sen <- sp::Lines(sen, id)
	sen <- sp::SpatialLines(list(sen), proj4string=crs)
	sen

}

#' @describeIn coordsToPoly Convert a list of coordinates to spatial lines
#' @export
coordsToLines <- function(x, crs, id=seq_along(x)) {

	if (class(crs) != 'CRS') crs <- sp::CRS(crs)

	n <- length(x)
	
	sens <- coordsToLine(x[[1]], crs=crs, id=id[1])
	
	if (n > 1) {
	
		for (i in 2:n) {
			this <- coordsToLine(x[[i]], crs, id=id[i])
			sens <- rbind(sens, this)
		}
		
	}
	
	sens
	
}
