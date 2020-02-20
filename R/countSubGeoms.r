#' Number of sub-geometries in a SpatialPolygons* or SpatiaLines* object
#'
#' This function returns the number of sub-polygons or lines in a SpatialPolygons* or SpatialLines* object.
#' @param x SpatialPolygons, SpatialPolygonsDataFrame, SpatialLines, or SpatialLinesDataFrame.
#' @seealso
#' @return An integer.
#' @examples
#' data(mad0)
#' countSubGeoms(mad0)
#' @export
countSubGeoms <- function(x) {

	cl <- class(x)
	if (cl %in% c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
		length(x@polygons[[1]]@Polygons)
	} else if (cl %in% c('SpatialLines', 'SpatialLinesDataFrame')) {
		length(x@lines)
	}

}
