#' Interpolate between spatial shapes
#'
#' This function generates a SpatialPolygons object that is the spatial "average" of two other SpatialPolygons or SpatialPolygonDataFrame objects. For example, assume you have two spatial polygons for ice over 12 Kybp and 10 Kybp and you wish to interpolate ice cover for 11 Kybp. Barring other information, the most parsimonious assumption is that the border of the ice 11 Kybp was exactly between the two other layers. Note that this function is different from taking values in two layers and calculating and average of them.
#' @param x1 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param x2 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param eaCRS Character string or object of class CRS. Coordinate reference system (proj4string) for an equal-area projection.
#' @param between Numeric between 0 an 1. Indicates the relative distance from the object in \code{x1} to \code{x2} represented by the interpolated layer. For example, if this value is 0.5, then the interpolated layer will have a border that is exactly between \code{x1} and \code{x2}. If values are between 0 and 0.5, they are closer to \code{x1}, and if between 0.5 and 1 they are closer to \code{x2}.
#' @param delta Positive numeric, precision of output. The polygon that results will be "true" to the "real" polygon (in a perfect world) to within this distance typically in meters).
#' @param verbose Logical, if \code{TRUE} then display progress.
#' @return A SpatialPolygons object.
#' @seealso \code{\link[birdsEye]{interPolysByBuffer}}
#' @examples
#' \dontrun{
#' xx1 <- c(44.02, 43.5, 42.61, 42.18, 42, 42.41, 42.75, 41.75, 41.49, 43.61,
#' 46.02, 46.5, 47.5, 47.39, 48.64, 49.05, 48.46, 48.18, 47.54, 46.73)
#' yy1 <- c(-18.83, -18.67, -18.87, -19.67, -20.65, -21.64, -23.08, -24.9,
#' -26.51, -27.09, -26.74, -25.6, -25.14, -26.44, -26.46, -24.96, -23.63, -22.72,
#' -23.36, -22.29)
#' xy1 <- cbind(xx1, yy1)
#' 
#' xx2 <- c(44.53, 44.18, 44.00, 42.93, 42.29, 42.71, 43.43, 45.15, 46.08, 45.94,
#' 45.36, 45.76, 46.97, 46.87, 45.94, 45.97, 45.08, 44.50, 44.58)
#' yy2 <- c(-24.27, -23.68, -22.86, -21.88, -20.56, -19.31, -20.36, -20.53, -20.93,
#' -21.81, -21.64, -22.90, -23.44, -24.08, -24.76, -25.95, -25.88, -25.61, -24.46)
#' xy2 <- cbind(xx2, yy2)
#' 
#' x1 <- coordsToPoly(xy1, enmSdm::getCRS('wgs84'))
#' x2 <- coordsToPoly(xy2, enmSdm::getCRS('wgs84'))
#' 
#' plot(x1)
#' plot(x2, border='blue', add=TRUE)
#' 
#' eaCrs <- enmSdm::getCRS('madAlbers')
#' 
#' inter <- interpPolysByDist(
#' 	x1, x2, eaCrs=eaCrs, between = 0.3, delta=10000
#' )
#' 
#' plot(inter, lty='dotted', add=TRUE)
#' }
#' @export
interpPolysByDist <- function(
	x1,
	x2,
	eaCrs,
	between = 0.5,
	delta = 1000,
	verbose = TRUE
) {

	### BAUHAUS
	### find intersection of each polygon in x1 and x2... if any:
		### for each sub-polygon in x1 with at least one x2 sub-polygon intersecting it
			### for each point on this sub-polygon of x1
				### find shortest line between this point and all points on the x2 sub-polygon
				### *** if this line crosses a border of x1 (then the x1 sub-polygon is locally convex to the x2 sub-polygon)
					### find all lines between the x1 point and all points in x2 sub-polygon
					### if all lines cross a border of the x1 sub-polygon, then abandon and go to next point in x1 sub-polygon
					### if not, select line with shortest distance
					### find point along this line with distance from the x1 sub-polygon proportionate to between... this will be a vertex of the interpolation polygon
				### if this line crosses a border of x1 (then the x1 sub-polygon is locally concave or parallel to the x2 sub-polygon)
					### find point along this line with distance from the x1 sub-polygon proportionate to between... this will be a vertex of the interpolation polygon
			### if any points found, create a spatial polygon
				### union with x2 sub-polygon and crop with x1
				### remember the polygon
	### find sub-polygons of x2 that intersect x1 AND has a portion falling outside of x1 (x2 partially covers x1)... if any:
		### for each sub-polygon of x2
			### for each vertex point in x2 sub-polygon
				### if point falls on border of intersection of x1 and x2, keep the point as-is (ie, no change between x1 and x2 at this point)
				### if point does not fall on border, find line with shortest distance to x1 sub-polygon
				### repeat from ***
	
	if (!sp::identicalCRS(x1, x2)) stop('Objects supplied for arguments "x1" and "x2" do not have the same coordinate reference system.')
	
	crs <- sp::CRS(sp::proj4string(x1))
	if (class(eaCrs) != 'CRS') eaCrs <- sp::CRS(eaCrs)
	
	numX1Polys <- countSubGeoms(x1)
	numX2Polys <- countSubGeoms(x2)
	
	### intersection of x1 and x2
	#############################
	
	# for each sub-polygon in x1
	for (countX1SubPolys in 1:numX1Polys) {

		# convert sub-polygon to polygon
		x1SubPolySp <- subGeomFromGeom(x1, countX1SubPolys)
		x1SubPolyCoordsSp <- geomToCoords(x1, countX1SubPolys, 'pts', removeClose=TRUE)

		# for each sub-polygon in x2
		for (countX2SubPolys in 1:numX2Polys) {
		
			if (verbose) say('Intersection of x1 sub-polygon #', countX1SubPolys, ' with x2 sub-polygon #', countX2SubPolys, '...')
		
			# convert sub-polygon to polygon
			x2SubPolySp <- subGeomFromGeom(x2, countX2SubPolys)
			x2SubPolyCoordsSp <- geomToCoords(x2, countX2SubPolys, 'pts', removeClose=TRUE)
			
			# if at least part of x2 sub-polygon is in x1 sub-polygon
			commonSp <- rgeos::gIntersection(x1SubPolySp, x2SubPolySp)

			if (!is.null(commonSp)) {
				
				commonPointsSp <- geomToCoords(commonSp, 1, 'pts', removeClose=TRUE)

				for (countX1Point in seq_along(x1SubPolyCoordsSp)) {

					thisStartPointSp <- x1SubPolyCoordsSp[countX1Point]
					thisStartPoint <- coordinates(thisStartPointSp)

					segEndsSp <- rgeos::gNearestPoints(thisStartPointSp, commonPointsSp)
					segEnds <- coordinates(segEndsSp)
					
					# if point is inside/on both x1 and x2 sub-polygons
					if (segEnds[1, 1] == segEnds[2, 1] & segEnds[1, 2] == segEnds[2, 2]) {
					
						# get point on interpolated polygon
						interpPoint <- coordinates(pointsAlong[closestPointToTarget])
						
						# remember point
						interps <- if (!exists('interps', inherits=FALSE)) {
							interpPoint
						} else {
							interps <- rbind(interps, interpPoint)
						}
						
					# if point is NOT inside/on both x1 and x2 sub-polygons
					} else {
						
						# does segment cross x1? (is x1 convex?)
						segLineSp <- coordsToLine(segEndsSp, crs=crs)
						convex <- gCrosses(segLineSp, x1)
						
						# line segment does not cross outsection... use this segment
						if (!convex) {
							
							# get desired distance along segment
							segDist <- lengthLine(segLineSp)
						
							subSegDist <- between * segDist
							
							# place regular points along the line, select one closest to desired distance from start
							pointsAlong <- sp::spsample(segLineSp, n=floor(segDist / delta), type='regular')
							distsToPoints <- geosphere::distm(segEndsSp[1], pointsAlong)
							closestPointToTarget <- which.min(abs(distsToPoints - subSegDist))
							
							# get point on interpolated polygon
							interpPoint <- coordinates(pointsAlong[closestPointToTarget])
							
							# remember point
							interps <- if (!exists('interps', inherits=FALSE)) {
								interpPoint
							} else {
								interps <- rbind(interps, interpPoint)
							}
						
						# if convex, try *all* lines from this point on x1	
						} else {
						
							# which line segments do not cross x1?
							allEnds <- list()
							for (countInCommonPoint in seq_along(commonPointsSp)) {
								allEnds[[countInCommonPoint]] <- rbind(
									thisStartPoint,
									sp::coordinates(commonPointsSp[countInCommonPoint])
								)
							}
							
							allEnds <- coordsToLines(allEnds, crs=crs)
							convex <- rgeos::gCrosses(allEnds, x1, byid=TRUE)
							allEnds <- allEnds[!convex]
							
							if (length(allEnds) > 0) {
							
								# find shortest line segment
								segDists <- geosphere::lengthLine(allEnds)
								shortest <- which.min(segDists)
								
								segEndsSp <- geomToCoords(allEnds, shortest, out='pts')
								segLineSp <- coordsToLine(segEndsSp, crs=crs)
							
								# get desired distance along segment
								segDist <- segDists[shortest]
								subSegDist <- between * segDist
								
								# place regular points along the line, select one closest to desired distance from start
								pointsAlong <- sp::spsample(segLineSp, n=floor(segDist / delta), type='regular')
								distsToPoints <- geosphere::distm(segEndsSp[1], pointsAlong)
								closestPointToTarget <- which.min(abs(distsToPoints - subSegDist))
								
								# get point on interpolated polygon
								interpPoint <- coordinates(pointsAlong[closestPointToTarget])
								
								# remember point
								interps <- if (!exists('interps', inherits=FALSE)) {
									interpPoint
								} else {
									interps <- rbind(interps, interpPoint)
								}
								
							}
						
						} # original segment between x1 and intersection of x1 and x2 crossed x1
						
					} # if point is NOT inside/on both x1 and x2 sub-polygons

				} # next point in x1 sub-polygon
				
				thisInterpSp <- coordsToPoly(interps, crs)
				thisInterpSp <- rgeos::gSimplify(thisInterpSp, tol=0)

				thisInterpSp <- rgeos::gUnion(thisInterpSp, commonSp)
				thisInterpSp <- rgeos::gIntersection(thisInterpSp, x1)

				rm(interps)
				
			} # if at least part of x2 sub-polygon is in x1 sub-polygon
			
			interpSp <- if (exists('interpSp')) {
				rgeos::gUnion(thisInterpSp, interpSp)
			} else {
				thisInterpSp
			}
			
		} # next x2 polygon
		
	} # next x1 polygon

	### part of intersection of x1 and x2 that is outside x1
	########################################################
	
	# for each sub-polygon in x1
	for (countX1SubPolys in 1:numX1Polys) {
	# for (countX1SubPolys in 1) {

		# convert sub-polygon to polygon
		x1SubPolySp <- subGeomFromGeom(x1, countX1SubPolys)
		x1SubPolyCoordsSp <- geomToCoords(x1, countX1SubPolys, 'pts', removeClose=TRUE)

		# add buffer to x1 sub-polygon for checking if points lie in it (obviates problem with some points that lie on an intersection but not on a vertex of x1 and x2 being considered to be outside x1)
		x1SpEa <- sp::spTransform(x1, eaCrs)
		x1SpEaBuff <- rgeos::gBuffer(x1SpEa, width=delta)
		x1SpBuff <- sp::spTransform(x1SpEaBuff, crs)
	
		# for each sub-polygon in x2
		for (countX2SubPolys in 1:numX2Polys) {
		
			if (verbose) say('Outsection of x1 sub-polygon ', countX1SubPolys, ' with x2 sub-polygon ', countX2SubPolys, '...')
		
			# convert sub-polygon to polygon
			x2SubPolySp <- subGeomFromGeom(x2, countX2SubPolys)
			x2SubPolyCoordsSp <- geomToCoords(x2, countX2SubPolys, 'pts', removeClose=TRUE)
			
			# if at least part of x2 sub-polygon is in x1 sub-polygon
			commonSp <- rgeos::gIntersection(x1SubPolySp, x2SubPolySp)
			notCommonSp <- rgeos::gDifference(x2SubPolySp, x1SubPolySp)

			if (!is.null(commonSp) & !is.null(notCommonSp)) {

				# add buffer to x2 sub-polygon for checking if points lie in it (obviates problem with some points that lie on an intersection but not on a vertex of x2 and x2 being considered to be outside x2)
				x2SpEa <- sp::spTransform(x2, eaCrs)
				x2SpEaBuff <- rgeos::gBuffer(x2SpEa, width=delta)
				x2SpBuff <- sp::spTransform(x2SpEaBuff, crs)

				# vertices of unshared area
				notCommonPointsSp <- geomToCoords(notCommonSp, n=NULL, 'pts', removeClose=TRUE)
				
				# for points in x1 and x2 sub-polygons, remember them as they are since they represent no change
				notCommonPointsSp_inX1 <- sp::over(notCommonPointsSp, x1SpBuff)
				notCommonPointsSp_inX2 <- sp::over(notCommonPointsSp, x2SpBuff)

				noChangeIndex <- which(notCommonPointsSp_inX1 == 1 & notCommonPointsSp_inX2 == 1)
								
				# for each point in outsection
				for (countNotCommonPoint in seq_along(notCommonPointsSp)) {
				
					thisStartPointSp <- notCommonPointsSp[countNotCommonPoint]
					thisStartPoint <- coordinates(thisStartPointSp)
				
					# if this point is inside/on x1 and x2 sub-polygons, remember it as-is
					if (any(countNotCommonPoint == noChangeIndex)) {
					
						# remember point
						interps <- if (!exists('interps', inherits=FALSE)) {
							thisStartPoint
						} else {
							interps <- rbind(interps, thisStartPoint)
						}
						
					# point is not on border of outsection
					} else {
				
						segEndsSp <- rgeos::gNearestPoints(x1SubPolyCoordsSp, thisStartPointSp)
						segEnds <- sp::coordinates(segEndsSp)
										
						# does segment cross x2? (is x2 convex?)
						segLineSp <- coordsToLine(segEndsSp, crs=crs)
						convex <- gCrosses(segLineSp, x2)
					
						if (!convex) {
						
							# get desired distance along segment
							segEndsSp <- geomToCoords(segLineSp, n=1, out='pts')
							segDist <- lengthLine(segLineSp)
						
							subSegDist <- (1 - between) * segDist
							
							# place regular points along the line, select one closest to desired distance from start
							pointsAlong <- sp::spsample(segLineSp, n=floor(segDist / delta), type='regular')
							distsToPoints <- geosphere::distm(segEndsSp[1], pointsAlong)
							closestPointToTarget <- which.min(abs(distsToPoints - subSegDist))
							
							# get point on interpolated polygon
							interpPoint <- coordinates(pointsAlong[closestPointToTarget])
							
							# remember point
							interps <- if (!exists('interps', inherits=FALSE)) {
								interpPoint
							} else {
								interps <- rbind(interps, interpPoint)
							}
					
						# if convex, try *all* lines from this point on x1	
						} else {
					
							# which line segments do not cross x1?
							allEnds <- list()
							for (countX1Point in seq_along(x1SubPolyCoordsSp)) {
								allEnds[[countX1Point]] <- rbind(
									thisStartPoint,
									sp::coordinates(x1SubPolyCoordsSp[countX1Point])
								)
							}
							
							allEnds <- coordsToLines(allEnds, crs=crs)
							convex <- rgeos::gCrosses(allEnds, x1, byid=TRUE)
							allEnds <- allEnds[!convex]
							
							if (length(allEnds) > 0) {
							
								segDists <- geosphere::lengthLine(allEnds)
								shortest <- which.min(segDists)
								
								segEndsSp <- geomToCoords(allEnds, shortest, out='pts')
								segLineSp <- coordsToLine(segEndsSp, crs=crs)
							
								segDist <- segDists[shortest]
								subSegDist <- (1 - between) * segDist
								
								pointsAlong <- sp::spsample(segLineSp, n=floor(segDist / delta), type='regular')
								distsToPoints <- geosphere::distm(segEndsSp[1], pointsAlong)
								closestPointToTarget <- which.min(abs(distsToPoints - subSegDist))
								
								interpPoint <- coordinates(pointsAlong[closestPointToTarget])
								
								# remember point
								interps <- if (!exists('interps', inherits=FALSE)) {
									interpPoint
								} else {
									interps <- rbind(interps, interpPoint)
								}
					
							}
							
						} # original segment between x1 and outsection of x1 and x2 crossed x2
							
					} # point is not on border of outsection

				} # next point in outsection
				
				thisInterpSp <- coordsToPoly(interps, crs)
				if (!gIsValid(thisInterpSp)) thisInterpSp <- rgeos::gSimplify(thisInterpSp, tol=0)

				thisInterpSp <- rgeos::gUnion(thisInterpSp, commonSp)

				rm(interps)
				
			} # if at least part of x2 sub-polygon is in x1 sub-polygon
			
			interpSp <- if (exists('interpSp')) {
				rgeos::gUnion(thisInterpSp, interpSp)
			} else {
				thisInterpSp
			}
			
		} # next x2 polygon
		
	} # next x1 polygon
	
	interpSp
	
}
