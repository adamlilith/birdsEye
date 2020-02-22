#' Extract vertex points from spatial polygon or lines object
#'
#' This function extracts coordinates of vertices of a polygon or line from a SpatialPolygon, SpatialPolygonDataFrame, SpatiaLines, or SpatialLinesDataFrame.
#' @param x SpatialPolygon, SpatialPolygonDataFrame, SpatiaLines, SpatialLinesDataFrame, SpatialPoints, or SpatialPointsDataFrame.
#' @param n \code{NULL} or positive integer(s) indicating which sub-polygon(s) or line(s) to extract. The default \code{NULL} is to extract points for all vertices of all sub-geometries. To see the number of sub-polygons or lines in an object use \code{\link[birdsEye]{countSubGeoms}}. This is ignored if the object in \code{x} is a SpatialPoints* object.
#' @param out Character, either \code{'matrix'} (default) or \code{'pts'}, indicating the returned object should either be a matrix with coordinates (default) or an object of class SpatialPoints
#' @param removeClose Logical, if \code{TRUE} (default), then remove the last pair of coordinates (for polygons only--ignored if lines). Most spatial polygons have the same set of first and last coordinates to complete the polygon.
#' @seealso \code{\link[birdsEye]{countSubGeoms}}, \code{\link[birdsEye]{subGeom}}
#' @return A two-column matrix or a SpatialPoints object.
#' @examples
#' data(mad0)
#' island <- subGeom(mad0, 1)
#' pts <- geomToCoords(island, 1)
#' plot(island)
#' points(pts)
#' @export
geomToCoords <- function(x, n = NULL, out='matrix', removeClose = TRUE) {

	cl <- class(x)
	objType <- if (cl %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
		'poly'
	} else if (cl %in% c('SpatialLines', 'SpatialLinesDataFrame')) {
		'line'
	} else if (cl %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
		'pts'
	} else {
		stop('Argument "x" must be a "Spatial*" (Lines, Points, Polygons) object.')
	}

	# if a points object
	if (objType == 'pts') {
	
		coords <- sp::coordinates(x)
		
	# if not a points object
	} else {
		
		# no sub-polygon specified
		toDo <- if (is.null(n)) {
			1:countSubGeoms(x)
		} else {
			n
		}

		# get coordinates from first sub-geometry
		sub <- subGeom(x, 1)
		if (objType == 'poly') {
			coords <- sub@polygons[[1]]@Polygons[[1]]@coords
			if (removeClose) coords <- coords[-nrow(coords), , drop=FALSE]
		} else {
			coords <- sp::coordinates(sub@lines[[1]])[[1]]
		}
		
		# if more than one sub-polygon, add points
		if (length(toDo) > 1) {
			
			for (i in 2:length(toDo)) {
			
				sub <- subGeom(x, i)
				if (objType == 'poly') {
					theseCoords <- sp::coordinates(sub@polygons[[1]]@Polygons[[1]])
					if (removeClose) theseCoords <- theseCoords[-nrow(theseCoords), , drop=FALSE]
				} else {
					theseCoords <- sp::coordinates(sub@lines[[i]])[[1]]
				}
				coords <- rbind(coords, theseCoords)
				
			}
		
		}
		
	} # if not a points object

	if (out == 'pts') {
	
		crs <- sp::proj4string(x)
		crs <- sp::CRS(crs)
		coords <- SpatialPoints(coords, proj4string=crs)
		
	}
	
	coords
	
}
