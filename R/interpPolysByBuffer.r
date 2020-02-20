#' Spatial interpolation between polygon borders
#'
#' This function interpolates between the borders of overlapping spatial polygons.
#' @param x1 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param x2 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param eaCrs Character string or object of class CRS. Coordinate reference system (proj4string) for an equal-area projection.
#' @param between Numeric between 0 an 1. For areas where \code{x2} is inside \code{x1}, then this is the relative distance from \code{x1} to \code{x2} to place the interpolated border (to a precision determined by \code{delta}). For areas where \code{x2} is outside \code{x1} but part of \code{x2} overlaps \code{x1}, this is relative distance from \code{x1} to the part of \code{x2} outside \code{x1} to place the interpolated polygon border. A value of 0.5 will place the interpolated border between \code{x1} and \code{x2}.
#' @param method Either \code{'grow'} or \code{'shrink'}. A growing interpolation uses a series of increasing buffers around the outside of \code{x2}. A shrinking buffer uses a series of smaller buffers inside \code{x1}.
#' @param delta Positive numeric, represents distance (typically in meters) by which to grow the buffer at each step. Smaller values yield more accurate interpolation but increase processing time.
#' @param quadsegs Positive integer, number of line segments to use to approximate a quarter circle. Larger numbers will yield more accurate results but also require more time.
#' @param verbose Logical or integer, if 0 or \code{FALSE} then display no progress indicators. If \code{TRUE} or 1 (default) then display some, if 2 or more display even more.
#' @return A SpatialPolygons object.
#' @seealso \code{\link[birdsEye]{interPolysByDist}}
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

interpPolysByBuffer <- function(
	x1,
	x2,
	eaCrs,
	between,
	method = 'grow',
	delta = 1000,
	quadsegs = 5,
	verbose = TRUE
) {

	crs <- sp::CRS(sp::proj4string(x1))
	if (class(eaCrs) != 'CRS') eaCrs <- sp::CRS(eaCrs)

	x1Ea <- sp::spTransform(x1, eaCrs)
	x2Ea <- sp::spTransform(x2, eaCrs)
	
	x1EaMinusDelta <- rgeos::gBuffer(x1Ea, width=delta / 10, byid=TRUE)
	x1MinusDelta <- sp::spTransform(x1EaMinusDelta, crs)

	### generate buffers
	
	if (verbose >= 1) omnibus::say('Generating trajectory contours...')
	
	buffs <- list()
	continue <- TRUE
	i <- 1
	
	while (continue) {

		# growing buffer
		if (method == 'grow') {
			
			buffSpEa <- rgeos::gBuffer(x2Ea, width=i * delta, byid=TRUE, quadsegs=quadsegs)
			continue <- !(rgeos::gContainsProperly(buffSpEa, x1Ea))
			buffSpEa <- rgeos::gIntersection(buffSpEa, x1Ea)
		
		# shrinking buffer
		} else {
		
			buffSpEa <- rgeos::gBuffer(x1Ea, width=-i * delta, byid=TRUE, quadsegs=quadsegs)
			contine <- !is.null(buffSpEa)
			if (continue) buffSpEa <- rgeos::gUnion(buffSpEa, x2Ea)
			
		}
		
		if (continue) {
			buffSp <- sp::spTransform(buffSpEa, crs)
			buffs[[i]] <- buffSp
			plot(buffSp, add=TRUE, lty='dotted')
			i <- i + 1
			
		}
	
	}
	
	numBuffs <- length(buffs)

	### for each vertex in x1, create a line that intercepts every point in every buffer such that each successive point is closest to the preceding point

	if (method == 'grow') {
		
		numSubsX1 <- countSubGeoms(x1)

		for (countSubsX1 in 1:numSubsX1) {
		
			if (verbose >= 1) omnibus::say('Processing sub-polygon #', countSubsX1, ' of ', numSubsX1, '...')
			if (verbose >= 2) progress <- txtProgressBar(min=0, max=numCoordsX1, style=3, width=min(64, getOption('width')))
			
			x1SubSp <- subGeomFromGeom(x1, countSubsX1)
			x1SubCoordsSp <- geomToCoords(x1SubSp, out='pts')
		
			numCoordsX1 <- length(x1SubCoordsSp)
		
			# for each point on x1, find trajectory to x2
			for (countSubPointsX1 in 1:numCoordsX1) {
	# say(countSubPointsX1)
				startCoordsSp <- x1SubCoordsSp[countSubPointsX1]
				lineCoords <- sp::coordinates(startCoordsSp)
				
				for (countBuff in numBuffs:1) {
				
					buffPointsSp <- geomToCoords(buffs[[countBuff]], out='pts')
					dists <- geosphere::distm(startCoordsSp, buffPointsSp)
					closest <- which.min(dists)
					
					lineCoords <- rbind(lineCoords, sp::coordinates(buffPointsSp[closest]))
				
				} # next buffer

				trajectSp <- coordsToLine(lineCoords, crs)
				trajectSp <- rgeos::gSimplify(trajectSp, tol=0)
				lineCoords <- geomToCoords(trajectSp)

				# trajectory is a single point
				if (nrow(lineCoords) == 2 && all(lineCoords[1, ] == lineCoords[2, ])) {

					interpPoint <- lineCoords[1, ]
					
				# trajectory is line segment(s)
				} else {
					
					# multi-lines object
					linesCoords <- list()
					for (i in 2:nrow(lineCoords)) {
						linesCoords[[i - 1]] <- rbind(lineCoords[i - 1, ], lineCoords[i, ])
					}
					
					# locate point closest to desired target length away from start point
					trajectSp <- coordsToLines(linesCoords, crs)
					trajectSpEa <- sp::spTransform(trajectSp, eaCrs)
					trajectDist <- rgeos::gLength(trajectSpEa)
					trajectDists <- rgeos::gLength(trajectSpEa, byid=TRUE)
					cumulDists <- cumsum(trajectDists)
					
					targetLength <- trajectDist * between
					
					interpPointIndex <- which.min(abs(targetLength - cumulDists))
					interpPoint <- lineCoords[interpPointIndex, , drop=FALSE]
					
				}
				
				thisInterps <- if (exists('thisInterps', inherits = FALSE)) {
					rbind(thisInterps, interpPoint)
				} else {
					interpPoint
				}
					
			} # next point on x1 sub-geometry
				
			thisInterpsSp <- coordsToPoly(thisInterps, crs=crs)
			if (!rgeos::gIsValid(thisInterpsSp)) thisInterpsSp <- rgeos::gSimplify(thisInterpsSp, tol=0)
			thisInterpsSp <- rgeos::gIntersection(thisInterpsSp, x1)
			thisInterpsSp <- rgeos::gUnion(thisInterpsSp, x2)
			
			rm(thisInterps)
			
			if (verbose >= 2) setTxtProgressBar(progress, countSubPointsX1)

		} # next x1 sub-geometry
		
		interpsSp <- if (exists('interpsSp', inherits=FALSE)) {
			rbind(interpsSp, thisInterpsSp)
		} else {
			thisInterpsSp
		}
		
	} # if "growing"
			
	interpsSp
	
}
