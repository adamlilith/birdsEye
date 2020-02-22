#' Extract a sub-polygon or coordinates thereof from a SpatialPolygons* object
#'
#' This function extracts a sub-polygon or the coordinates of vertices a sub-polygon from a SpatialPolygon or SpatialPolygonDataFrame object.
#' @param x SpatialPolygon or SpatialPolygonDataFrame.
#' @param n Positive integer(s) indicating which sub-polygon(s) to obtain. To see the number of sub-polygons in an object use \code{\link[birdsEye]{countSubGeoms}}. If \code{NULL}, then the function basically returns \code{x}, but with all sub-geometries split (i.e., a multi-part polygon/lines).
#' @param out Character, either \code{'poly'} (default) or \code{'coords'}, indicating the returned object should either be a spatial polygon or a matrix of coordinates.
#' @seealso \code{\link[birdsEye]{countSubGeoms}}, \code{\link[birdsEye]{geomToCoords}}
#' @return A spatial polygon object or a two-column matrix.
#' @examples
#' data(mad0)
#' sub <- subGeom(mad0, 1)
#' pts <- geomToCoords(mad0, 1)
#' plot(sub)
#' points(pts)
#' 
#' N <- countSubGeoms(mad0)
#' sub <- subGeom(mad0, N)
#' plot(sub)
#'
#' par(mfrow=c(4, 4))
#' for (n in 1:16) {
#' sub <- subGeom(mad0, n)
#' plot(sub)
#' }
#' @export
subGeom <- function(x, n) {

	crs <- sp::proj4string(x)
	crs <- sp::CRS(crs)
	
	cl <- class(x)

	if (is.null(n)) n <- 1:countSubGeoms(x)
	
	if (cl %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
	
		# first polygon
		out <- sp::coordinates(x@polygons[[1]]@Polygons[[n[1]]])
		out <- coordsToPoly(out, crs, id=1)
		
		# remaining polygons
		if (length(n) > 1) {
		
			for (i in n[2]:n[length(n)]) {
			
				thisOut <- sp::coordinates(x@polygons[[1]]@Polygons[[i]])
				thisOut <- coordsToPoly(thisOut, crs, id=i)
				out <- rbind(out, thisOut)
				
			}
			
		}
		
	} else if (cl %in% c('SpatialLines', 'SpatialLinesDataFrame')) {
	
		out <- sp::coordinates(x@lines[[n[1]]])[[1]]
		out <- coordsToLine(out, crs, id=1)
		
		# remaining polygons
		if (length(n) > 1) {
		
			for (i in n[2]:n[length(n)]) {
			
				thisOut <- sp::coordinates(x@lines[[i]])[[1]]
				thisOut <- coordsToLine(thisOut, crs, id=i)
				out <- rbind(out, thisOut)
				
			}
			
		}
		
	}
	
	out
	
}
