#' Extract a sub-polygon or coordinates thereof from a SpatialPolygons* object
#'
#' This function extracts a sub-polygon or the coordinates of vertices a sub-polygon from a SpatialPolygon or SpatialPolygonDataFrame object.
#' @param x SpatialPolygon or SpatialPolygonDataFrame.
#' @param n Positive integer indicating which sub-polygon to obtain. To see the number of sub-polygons in an object use \code{\link[enmSdm]{countSubGeoms}}.
#' @param out Character, either \code{'poly'} (default) or \code{'coords'}, indicating the returned object should either be a spatial polygon or a matrix of coordinates.
#' @seealso \code{\link[enmSdm]{countSubGeoms}}, \code{\link[enmSdm]{geomToCoords}}
#' @return A spatial polygon object or a two-column matrix.
#' @examples
#' data(mad0)
#' sub <- subGeomFromGeom(mad0, 1)
#' pts <- geomToCoords(mad0, 1)
#' plot(sub)
#' points(pts)
#' 
#' N <- countSubGeoms(mad0)
#' sub <- subGeomFromGeom(mad0, N)
#' plot(sub)
#'
#' par(mfrow=c(4, 4))
#' for (n in 1:16) {
#' sub <- subGeomFromGeom(mad0, n)
#' plot(sub)
#' }
#' @export
subGeomFromGeom <- function(x, n) {

	crs <- sp::proj4string(x)
	crs <- sp::CRS(crs)
	
	cl <- class(x)
	
	if (cl %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
	
		out <- sp::coordinates(x@polygons[[1]]@Polygons[[n]])
		out <- coordsToPoly(out, crs)
		
	} else if (cl %in% c('SpatialLines', 'SpatialLinesDataFrame')) {
	
		out <- sp::coordinates(x@lines[[n]])[[1]]
		out <- coordsToLine(out, crs)
	}
	
	out
	
}
