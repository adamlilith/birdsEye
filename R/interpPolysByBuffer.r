#' Spatial interpolation between polygon borders using buffers
#'
#' This function recreates a stage in the morphing of one SpatialPolygon* to another. The output is a SpatialPolygons object with borders that are "between" the borders of two other SpatialPolygons* objects. This particular function uses inner/outer buffers to estimate the location of polygon vertices.
#' @param x1 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param x2 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param eaCrs Character string or object of class CRS. Coordinate reference system (proj4string) for an equal-area projection.
#' @param between Numeric between 0 and 1. For areas where \code{x2} is inside \code{x1}, then this is the relative distance from \code{x1} to \code{x2} to place the interpolated border (to a precision determined by \code{delta}). For areas where \code{x2} is outside \code{x1} but part of \code{x2} overlaps \code{x1}, this is relative distance from \code{x1} to the part of \code{x2} outside \code{x1} to place the interpolated polygon border. A value of 0.5 will place the interpolated border between \code{x1} and \code{x2}.
#' @param method Either \code{'grow'} or \code{'shrink'}. A growing interpolation uses a series of increasing buffers around the outside of \code{x2}. A shrinking buffer uses a series of smaller buffers inside \code{x1}.
#' @param delta Positive numeric, represents distance (typically in meters) by which to grow the buffer at each step. Smaller values yield more accurate interpolation but increase processing time.
#' @param quadsegs Positive integer, number of line segments to use to approximate a quarter circle. Larger numbers will yield more accurate results but also require more time.
#' @param verbose Logical or integer, if 0 or \code{FALSE} then display no progress indicators. If \code{TRUE} or 1 (default) then display some, if 2 or more display even more.
#' @details This function approximates an interpolation... It will not produce a mathematically "elegant" one. A typical application for this function would be to interpolate glacial ice layers from layers that represent ice cover at coarse temporal resolution (e.g., estimating ice coverage at 8500 Kybp given layers at 8000 and 9000 Kybp). As a result, it is assumed that all subgeometries in \code{x2} at least partially overlap a subgeometry in \code{x1} (e.g., \code{x2} represents ice more recently and \code{x1} ice longer ago, so \code{x2} is mostly a spatial subset of \code{x1}). The function operates using the following procedure:
#' \enumerate{
#' 	\item For each subpolygon of \code{x1}:
#'	\item Find all subpolygons of \code{x2} that overlap it or occur in it (if none, skip to step 8 below).
#'	\item For each these subpolygons in \code{x2}, obtain the union between it and the \code{x1} subpolygon.
#'	\item Create a series of nested buffers that can expand from \code{x2} (\code{method='grow'}; i.e., an "outside" buffer on \code{x2}) or contract from \code{x1} (\code{method='shrink'}; i.e., an "inside" buffer on \code{x1}). Buffers will be multiples of \code{delta} in distance from the starting polygon. The buffers are cropped to the intersection of the \code{x2} and \code{x2} sub-polygons. Union the buffers with \code{x2} so they always contain \code{x2} (i.e., they do not get smaller than \code{x2}).
#'	\item Create a "trajectory" for each vertex of \code{x1}: Find the closest point on the buffer nearest it. Then from this point, find the closest point on the next buffer, and so on, until all buffers have been included. Add the point on \code{x2} that is closest to the last point to the trajectory. The trajectory should trace a path from the vertex of \code{x1} to a vertex of \code{x2}.
#' 	\item Calculate the length of each segment of the trajectory, starting at the point on \code{x1}. Calculate the length of the length of the entire trajectory, then the length along the trajectory one wishes to go.  This length is given by \code{between} * trajectory length. Find the point on the trajectory (which occurs on one of the buffers or on \code{x2}) that is closest to the desired length. This point represents a vertex of the interpolated polygon. Repeat for all vertices of \code{x1}. From these points create a polygon to be unioned with any preceding polygons made in the same way. This will be unioned with polygons created in subsequent steps to create the interpolated surface.
#'	\item If part of a \code{x2} subpolygon that overlaps an \code{x1} subpolygon also partially falls outside \code{x1}, calculate a series of buffer either growing from \code{x1} or shrining from \code{x2}. Crop the buffers to the portion of \code{x2} that falls outside \code{x1}. For each point on \code{x1} calculate trajectories as described, but use \code{(1 - between)} * trajectory length as the target length. The idea here is that \code{x1} had to expand to become \code{x2}. Then calculate the interpolated polygon.
#'	\item If the subpolygon of \code{x1} does not overlap or contain any part of \code{x2}, then calculate its geographic centroid. Calculate a series of buffers that "grow" the centroid or "shrink" from the borders of the \code{x2} subpolygon. Calculate trajectories and then the interpolated polygon as before. The idea here is that the subpolygon of \code{x2} completely disappears (i.e., no portion of it becomes \code{x2}), so it has to shrink to do that.
#'	\item When finished, union the resulting polygons with \code{x2}.
#' }
#' Note: If a subpolygon is too small to contain a single buffer, then it is assumed to completely disappear.\cr
#' The function can give unexpected results if the geometries of the input polygons are complicated enough if, For example, there are holes or portions that are nearly holes (i.e., the polygon nearly closes on itself). Also note that if using \code{method = 'shrink'}, then the trajectories are not guaranteed to neatly converge on the relevant parts of \code{x2} if the latter are not near the inside-most buffer. If weird results ensue, try reducing \code{delta} or increasing \code{quadsegs} (doing either increases processing time).
#' @return SpatialPolygons object.
#' @examples
#' \donttest{
#' # create "x1": has two sub-polygons
#' x <- c(44.02, 43.5, 42.61, 42.18, 42, 42.41, 42.75, 41.75, 41.49,
#' 43.61,46.02, 46.5, 47.5, 47.39, 48.64, 49.05, 48.46, 48.18, 47.54, 46.73, 45.80, 45.59)
#' y <- c(-18.83, -18.67, -18.87, -19.67, -20.65, -21.64, -23.08, -24.9,
#' -26.51, -27.09, -26.74, -25.6, -25.14, -26.44, -26.46, -24.96, -23.63,
#'  -22.72, -23.36, -22.29, -21.45, -20.69)
#' xy1a <- cbind(x, y)
#' 
#' x <- c(40.61, 40.07, 40.23, 41.38, 41.38)
#' y <- c(-20.51, -20.49, -21.11, -21.55, -21.01)
#' xy1b <- cbind(x, y)
#' 
#' x1a <- coordsToPoly(xy1a, enmSdm::getCRS('wgs84'))
#' x1b <- coordsToPoly(xy1b, enmSdm::getCRS('wgs84'))
#' x1 <- rgeos::gUnion(x1a, x1b)
#' 
#' # create "x2"
#' x <- c(44.53, 44.18, 44.00, 42.93, 42.29, 42.71, 43.43, 47.15, 48.08,
#'  45.94,45.36, 45.76, 46.97, 46.87, 45.94, 45.97, 45.08, 44.50, 44.58)
#' y <- c(-24.27, -23.68, -22.86, -21.88, -20.56, -19.31, -20.36, -20.53,
#'  -20.93,-21.81, -21.64, -22.90, -23.44, -24.08, -24.76, -25.95, -25.88, -25.61, -24.46)
#' xy2 <- cbind(x, y)
#' 
#' x2 <- coordsToPoly(xy2, enmSdm::getCRS('wgs84'))
#' 
#' eaCrs <- enmSdm::getCRS('albersNA')
#' 
#' interBuff <- interpPolysByBuffer(
#' 	x1, x2, eaCrs=eaCrs, between = 0.4, delta=10000
#' )
#' 
#' interTween <- interpPolysByTween(
#' 	x1, x2, eaCrs='laea', between = 0.4, delta=100
#' )
#'
#' interTween <- interTween[[1]]$poly
#' 
#' plot(x1, col='gray90')
#' plot(x2, add=TRUE)
#' plot(interBuff, border='red', add=TRUE)
#' plot(interTween, border='green', add=TRUE)
#' legend('bottomleft',
#' 	legend=c('x1', 'x2', 'by buffer', 'by tween'),
#' 	fill=c('gray90', NA, NA, NA),
#' 	border=c('black', 'black', 'red', 'green'),
#' 	bty='n'
#' )
#'
#' ## multiple steps
#' between <- seq(0, 1, by=0.1)
#' interTween <- interpPolysByTween(
#' 	x1, x2, eaCrs='laea', between = between, delta=100
#' )
#' 
#' plot(x1, col='gray90')
#' plot(x2, add=TRUE)
#' for (i in seq_along(between)) {
#' 	plot(interTween[[i]]$poly, border='green', lty='dotted', add=TRUE)
#' }
#' 
#' legend('bottomleft',
#' 	legend=c('x1', 'x2', 'tweens'),
#' 	fill=c('gray90', NA, NA),
#' 	border=c('black', 'black', 'green'),
#' 	bty='n'
#' )
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

	### for each subgeometry of x1 find all x2 subgeoms that touch it
		### if any x1 subgeom and any part of x2 are proximate (intersect, contain, share border)
			### subset the x1 subgeometry and the parts of x2 that are not disjoint with the x1 subgeometry
			### for parts inside x1 calculate buffer to x2
				### calculate trajectory for each point on x1 to x2
			### for parts of x2 outside x1
				### calculate buffer inside x2 shrinking to x1
				### calculate trajectory for each point on x1 to x2
		
		### if none
			### calculate inward buffers
			### place fake x2 at centroid and calculate trajectories

		# countSubsX1 = 6 for i = 3 is singlet	
		# countSubsX1 <- 11 for multi-x2

	if (!sp::identicalCRS(x1, x2)) stop('"x1" and "x2" must have the same coordinate reference system.')
		
	crs <- sp::CRS(sp::proj4string(x1))
	if (class(eaCrs) != 'CRS') eaCrs <- sp::CRS(eaCrs)

	x1Ea <- sp::spTransform(x1, eaCrs)
	x2Ea <- sp::spTransform(x2, eaCrs)
	# x1x2Ea <- rgeos::gUnion(x1Ea, x2Ea)

	numSubsX1 <- countSubGeoms(x1)
	numSubsX2 <- countSubGeoms(x2)
	
	x2AsSubGeoms <- subGeom(x2, n=NULL)

	### for each subgeometry of x1 find all x2 subgeoms that touch it
	#################################################################
	
	for (countSubsX1 in 1:numSubsX1) {

		if (verbose) omnibus::say('Processing subpolygon ', countSubsX1, ' of ', numSubsX1, '...')

		# x1 subgeometry
		x1Sub <- subGeom(x1, countSubsX1)
		x1SubCoords <- geomToCoords(x1Sub)
		x1SubEa <- sp::spTransform(x1Sub, eaCrs)
		
		# x1SubAndX2Proximate <- rgeos::gIntersect(x1Sub, x2) | rgeos::gContains(x1Sub, x2) | rgeos::gTouches(x1Sub, x2)
		x1x2SubsProximate <- any(!is.na(sp::over(x1Sub, x2)))

		### if x1 subgeometry overlaps any part of x2
		#############################################
		
		if (x1x2SubsProximate) {

			### subset the x1 subgeometry and the parts of x2 that are not disjoint with the x1 subgeometry
			x2Proximates <- which(!is.na(sp::over(x2AsSubGeoms, x1Sub)))

			## for each subgeometry in x2 that is proximate to this x2 subgeometry
			for (thisX2Proximate in x2Proximates) {

				x2Sub <- subGeom(x2, thisX2Proximate)
				x2SubEa <- sp::spTransform(x2Sub, eaCrs)
					
				### for parts inside x1 calculate buffer from x1 to x2
				commonEa <- rgeos::gUnion(x2SubEa, x1SubEa)
				
				# generate series of buffers
				theseBuffs <- .generateBuffsPolyVsPoly(from=x1SubEa, to=x2SubEa, demesne=commonEa, method=method, delta=delta, crs=crs, quadsegs=quadsegs)
				
				### calculate trajectory for each point on x1 to x2
				
				# for each vertex of x1 subgeometry
				if (exists('thisInterps', inherits=FALSE)) rm(thisInterps)
				for (countX1Point in 1:nrow(x1SubCoords)) {

					interpPoint <- .locateInterpPoint(start=x1SubCoords[countX1Point, , drop=FALSE], end=x2Sub, buffs=theseBuffs, between=between, crs=crs)

					thisInterps <- if (exists('thisInterps', inherits=FALSE)) {
						rbind(thisInterps, interpPoint)
					} else {
						interpPoint
					}

				} # next point in subgeometry of x1
			
				# process interpolated points, convert to polygon
				thisInterps <- thisInterps[!duplicated(thisInterps), , drop=FALSE]
				
				if (nrow(thisInterps) > 1) {
					
					thisInterpsSp <- coordsToPoly(thisInterps, crs)
					thisInterpsSp <- rgeos::gSimplify(thisInterpsSp, tol=0)
					thisInterpsSp <- rgeos::gIntersection(thisInterpsSp, x1Sub)
					
					# remember
					interpsSp <- if (exists('interpsSp', inherits = FALSE)) {
						rgeos::gUnion(interpsSp, thisInterpsSp)
					} else {
						thisInterpsSp
					}
					
				}
										
				### if x2 has any part that lies outside x1 subgeometry (but also inside), "grow" x1 into x2 or "shrink" x2 to x1
				#################################################################################################################
				
				x2Outside <- rgeos::gDifference(x2Sub, x1Sub)
				
				if (!is.null(x2Outside)) {

					# for each component of x2 subgeometry outside x1 subgeometry (ie, sub-sub-geometries of x2)
					for (x2OutsideSubs in 1:countSubGeoms(x2Outside)) {
						
						x2OutsideSub <- subGeom(x2Outside, x2OutsideSubs)
						x2OutsideSubEa <- sp::spTransform(x2OutsideSub, eaCrs)
						x2OutsideSubCoords <- geomToCoords(x2OutsideSub)

						# # generate series of buffers
						# theseBuffs <- if (method %in% c('grow', 'auto')) {
							# .generateBuffsPolyVsPoly(from=x1SubEa, to=x2OutsideSubEa, demesne=x2OutsideSubEa, method='grow', delta=delta, quadsegs=quadsegs)
						# } else if (method == 'shrink') {
							# .generateBuffsPolyVsPoly(from=x2OutsideSubEa, to=x1SubEa, demesne=x2OutsideSubEa, method='shrink', delta=delta, quadsegs=quadsegs)
						# }
						
						theseBuffs <- .generateBuffsPolyVsPoly(from=x1SubEa, to=x2OutsideSubEa, demesne=x2OutsideSubEa, method='grow', delta=delta, crs=crs, quadsegs=quadsegs)

						### calculate trajectory for each point on x2 to x1
						
						# for each vertex of x2 subgeometry
						if (exists('thisInterps', inherits=FALSE)) rm(thisInterps)
						for (countX2Point in 1:nrow(x2OutsideSubCoords)) {

							interpPoint <- .locateInterpPoint(start=x2OutsideSubCoords[countX2Point, , drop=FALSE], end=x1Sub, buffs=theseBuffs, between=1 - between, crs=crs)

							thisInterps <- if (exists('thisInterps', inherits=FALSE)) {
								rbind(thisInterps, interpPoint)
							} else {
								interpPoint
							}

						} # next point in subgeometry of x2
					
						# process interpolated points, convert to polygon
						thisInterps <- thisInterps[!duplicated(thisInterps), , drop=FALSE]
						
						if (nrow(thisInterps) > 1) {
							
							thisInterpsSp <- coordsToPoly(thisInterps, crs)
							thisInterpsSp <- rgeos::gSimplify(thisInterpsSp, tol=0)
							# thisInterpsSp <- rgeos::gUnion(thisInterpsSp, x2Sub)
							
							# remember
							interpsSp <- if (exists('interpsSp', inherits = FALSE)) {
								rgeos::gUnion(interpsSp, thisInterpsSp)
							} else {
								thisInterpsSp
							}
							
						}
						
					} # for each component of x2 subgeometry outside x1 subgeometry (ie, sub-sub-geometries of x2)
				
				
				} # if any portion of this x2 subgeometry falls outside this x1 subgeometry
				
			} ## for each subgeometry in x2 that is proximate to this x2 subgeometry
			
		### if x1 subgeometry has no overlap with any part of x2
		########################################################
		
		} else if (!x1x2SubsProximate) {

			### generate series of buffers
			centX1SpEa <- rgeos::gCentroid(x1SubEa)
			centX1Sp <- sp::spTransform(centX1SpEa, crs)
			startAt <- if (method == 'grow') { 'pt' } else if (method %in% c('auto', 'shrink')) { 'poly' }
			theseBuffs <- .generateBuffsPointVsPoly(pt=centX1SpEa, poly=x1SubEa, startAt=startAt, method=method, delta=delta, crs=crs, quadsegs=quadsegs)
		
			# if the subgeometry was not too small to place buffers, calculate interpolation... otherwise, assume it disappears
			if (length(theseBuffs) > 0) {
			
				# for each vertex of x1 subgeometry
				if (exists('thisInterps', inherits=FALSE)) rm(thisInterps)
				for (countX1Point in 1:nrow(x1SubCoords)) {

					# note: using 1 minus between!
					interpPoint <- .locateInterpPoint(start=x1SubCoords[countX1Point, , drop=FALSE], end=centX1Sp, buffs=theseBuffs, between=between, crs=crs)

					thisInterps <- if (exists('thisInterps', inherits=FALSE)) {
						rbind(thisInterps, interpPoint)
					} else {
						interpPoint
					}

				} # next point in subgeometry of x1
			
				# process interpolated points, convert to polygon
				thisInterps <- thisInterps[!duplicated(thisInterps), , drop=FALSE]
				
				if (nrow(thisInterps) > 1) {
					
					thisInterpsSp <- coordsToPoly(thisInterps, crs)
					thisInterpsSp <- rgeos::gSimplify(thisInterpsSp, tol=0)
					thisInterpsSp <- rgeos::gIntersection(thisInterpsSp, x1Sub)
					
					# remember
					interpsSp <- if (exists('interpsSp', inherits = FALSE)) {
						rgeos::gUnion(interpsSp, thisInterpsSp)
					} else {
						thisInterpsSp
					}
					
				}
			
			} # x1 subgeometry is not associated with any x2
			
		} # if the subgeometry was not too small to place buffers, calculate interpolation... otherwise, assume it disappears
			
	} # next x1 subgeometry

	# interpsSp <- rgeos::gUnion(interpsSp, x2)
	interpsSp <- rgeos::gSimplify(interpsSp, tol=0)
	x1x2 <- rgeos::gUnion(x1, x2)
	interpsSp <- rgeos::gIntersection(interpsSp, x1x2)
	interpsSp
	
}



	### generate a set of buffers inside one polygon based on another polygon
	.generateBuffsPolyVsPoly <- compiler::cmpfun(function(from, to, demesne, method, delta, crs, quadsegs=5) {

		# from		spatial object in equal-area CRS, growing/shrinking *from* this object
		# to		spatial object in equal-area CRS, growing/shrinking *to* this object
		# demesne	area buffers need to encompass before halting (if growing)
		# method	"grow" or "auto" (from "from") OR "shrink" (from "from")
		# delta		distance multiplier for each buffer
		# crs		object of class CRS
		# quadsegs	number of segments used to approximate a quarter-circle

		buffs <- list()
		continue <- TRUE
		i <- 1

		# demesne <- rgeos::gBuffer(demesne, width=-delta / 10)
		
		while (continue) {

			# growing buffer
			if (method %in% c('grow', 'auto')) {
				
				buffSpEa <- rgeos::gBuffer(to, width=i * delta, quadsegs=quadsegs)
				continue <- !(rgeos::gContainsProperly(buffSpEa, demesne))
				buffSpEa <- rgeos::gIntersection(buffSpEa, demesne)
				continue <- continue & !is.null(buffSpEa)
			
			# shrinking buffer
			} else if (method == 'shrink') {
			
				buffSpEa <- rgeos::gBuffer(from, width=-i * delta, quadsegs=quadsegs)
				continue <- !is.null(buffSpEa)
				if (continue) buffSpEa <- rgeos::gUnion(buffSpEa, to)
				
			} else {
				stop('Argument "method" must be "grow", "shrink", or "auto".')
			}
			
			if (continue && !is.null(buffSpEa)) {
				buffSp <- sp::spTransform(buffSpEa, crs)
				buffs[[i]] <- buffSp
		# plot(buffSp, add=TRUE, lty='dotted')
				i <- i + 1
				
			}
		
		}
		
		# if growing the buffer then buffers were created in reverse order (from inset to container), so reverse them
		if (method %in% c('grow', 'auto')) buffs <- rev(buffs)
		
		buffs
		
	})

	### generate a set of buffers inside a polygon relative to a point
	.generateBuffsPointVsPoly <- compiler::cmpfun(function(pt, poly, startAt, method, delta, crs, quadsegs=5) {

		# pt		spatial point in equal-area CRS, must occur *inside* or on the polygon
		# poly		spatial polygon(s) in equal-area CRS
		# startAt	"poly" (shrink from poly), OR "pt" or "auto" (grow from point)
		# method	'grow', 'shrink', 'auto'
		# delta		distance multiplier for each buffer
		# crs		object of class CRS
		# quadsegs	number of segments used to approximate a quarter-circle

		buffs <- list()
		continue <- TRUE
		i <- 1

		while (continue) {

			# growing buffer
			if (startAt %in% c('pt', 'auto')) {
				
				buffSpEa <- rgeos::gBuffer(pt, width=i * delta, byid=TRUE, quadsegs=quadsegs)
				continue <- !(rgeos::gContainsProperly(buffSpEa, poly))
				buffSpEa <- rgeos::gIntersection(buffSpEa, poly)
			
			# shrinking buffer
			} else if (startAt == 'poly') {
			
				buffSpEa <- rgeos::gBuffer(poly, width=-i * delta, byid=TRUE, quadsegs=quadsegs)
				continue <- !is.null(buffSpEa)
				
			} else {
				stop('Argument "startAt" must be "poly" or "pt",')
			}
			
			if (continue) {
				buffSp <- sp::spTransform(buffSpEa, crs)
				buffs[[i]] <- buffSp
		# plot(buffSp, add=TRUE, lty='dotted')
				i <- i + 1
				
			}
		
		}
		
		if (method %in% c('grow', 'auto')) buffs <- rev(buffs)
		buffs
		
	})

	### for a given starting point and a set of buffers, calculate the coordinate of the interpolation
	.locateInterpPoint <- compiler::cmpfun(function(start, end, method, buffs, between, crs) {
	
		# start			2-column matrix of the *coordinates* of the starting point
		# end			polygon/points of end point candidates (unprojected)
		# buffs			list of buffers
		# between		proportionate distance along trajectory where to locate interpolated point
		# crs			object of class CRS
		
		# get trajectory coordinates
		traject <- start
		endCoords <- geomToCoords(end)

		# while trajectory has not ended on the end geometry
		countBuff <- 1
		continue <- TRUE
		stopOnEnd <- FALSE
		
		# if there are any buffers
		if (length(buffs) > 0) {
			
			while (continue) {
			
				# all distances from this point to nearest buffer
				thisBuff <- buffs[[countBuff]]
				thisBuffCoords <- geomToCoords(thisBuff)
				distsToBuff <- geosphere::distm(traject[countBuff, , drop=FALSE], thisBuffCoords)
				distsToEnd <- geosphere::distm(traject[countBuff, , drop=FALSE], endCoords)
				
				# closest point is on a buffer
				if (min(distsToBuff) < min(distsToEnd)) {
				
					closest <- which.min(distsToBuff)
					traject <- rbind(traject, thisBuffCoords[closest, , drop=FALSE])
					countBuff <- countBuff + 1
					continue <- (countBuff <= length(buffs))
					stopOnEnd <- FALSE
				
				# closest point is on the end geometry
				} else {
					
					closest <- which.min(distsToEnd)
					traject <- rbind(traject, endCoords[closest, , drop=FALSE])
					continue <- FALSE
					stopOnEnd <- TRUE
					
				}
				
			} # next buffer
			
		}

		# if trajectory did not end on the end geometry (ended on a buffer), find last point on end geometry
		if (!stopOnEnd) {
			
			endPoints <- geomToCoords(end)
			dists <- geosphere::distm(traject[nrow(traject), , drop=FALSE], endPoints)
			closest <- which.min(dists)
			traject <- rbind(traject, endPoints[closest, , drop=FALSE])
			
		}
			
		traject <- traject[!duplicated(traject), , drop=FALSE]
		
		# trajectory is a single point
		if (nrow(traject) == 1) {
		
			out <- traject
			
		# trajectory is a path
		} else {
			
			# convert trajectory to multi-lines object
			trajectLines <- list()
			for (i in 2:nrow(traject)) {
				trajectLines[[i - 1]] <- rbind(traject[i - 1, ], traject[i, ])
			}
			
			trajectSp <- coordsToLines(trajectLines, crs)
			
			trajectSpEa <- sp::spTransform(trajectSp, eaCrs)
			trajectDists <- rgeos::gLength(trajectSpEa, byid=TRUE)
			trajectDist <- sum(trajectDists)
			cumulDists <- cumsum(trajectDists)
			cumulDists <- c(0, cumulDists)
			
			targetLength <- trajectDist * between
			
			interpPointIndex <- which.min(abs(targetLength - cumulDists))
			out <- traject[interpPointIndex, , drop=FALSE]
			
		} # if trajectory is a path
		
		out
			
	})
