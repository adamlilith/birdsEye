#' Spatial interpolation between polygon borders using tweening
#'
#' This function recreates stages in the morphing of one SpatialPolygon* to another. The output is a SpatialPolygons object with borders that are "between" the borders of two other SpatialPolygons* objects. This particular function uses spatial tweening to estimate vertices of polygons.
#' @param x1 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param x2 SpatialPolygon or SpatialPolygonDataFrame object in an unprojected (WGS84) coordinate reference system.
#' @param eaCrs This is either a proj4 string, an object of class \code{CRS}, or an abbreviated name of equal-area projection to use. The polygons will be projected to this coordinate reference system before interpolation. Ideally, the center point of the projection should be in the center of the polygons to minimize distortion. Options include:
#' \itemize{
#'	\item \code{'laea'}: Lambert azimuthal equal-area
#'  \item \code{'mollweide'}: Mollweide (equal-area)
#' }
#' @param between Numeric between 0 and 1. This is the relative distance from \code{x1} to \code{x2} to place the interpolated border (higher values of \code{delta} increase precision but also increase computational time). A value of 0 should return a polygon the same as \code{x1} and a value of 1 should return a polygon the same as \code{x2}.
#' @param method method used to optimize tweened polygons. Either \code{'cubic-in-out'} or \code{'linear'}.
#' @param delta Positive numeric, represents distance (typically in meters) by which to grow the buffer at each step. Smaller values yield more accurate interpolation but increase processing time.
#' @param verbose Logical. If \code{TRUE} then display progress indicators.
#' @examples
#' \dontrun{
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
#' }
#' @export

interpPolysByTween <- function(
	x1,
	x2,
	eaCrs,
	between,
	delta = 20,
	method = 'cubic-in-out',
	verbose = TRUE
) {

	if (!sp::identicalCRS(x1, x2)) stop('"x1" and "x2" must have the same coordinate reference system.')

	### create equal-area CRS with centroid of polygons as its center
	#################################################################

		# x1 <- subGeom(x1, 11)
		# x2 <- subGeom(x2, 1:countSubGeoms(x2))
	
		x1x2 <- rgeos::gUnion(x1, x2)
		centSp <- geosphere::centroid(x1x2)
		
		crs <- sp::proj4string(x1)
		
		long0 <- sp::coordinates(centSp)[1, 1]
		lat0 <- sp::coordinates(centSp)[1, 2]
		if (class(eaCrs) != 'CRS') {
			if (eaCrs %in% c('laea', 'mollweide')) {
				ell <- ellipsoid(crs)
				dat <- datum(crs)
				eaCrs <- enmSdm::makeCRS(eaCrs, long0=long0, lat0=lat0, dat=dat, ell=ell, asCRS=TRUE)
			} else {
				eaCrs <- sp::CRS(eaCrs)
			}
		}
		
		# project polygons to equal-area
		x1Ea <- sp::spTransform(x1, eaCrs)
		x2Ea <- sp::spTransform(x2, eaCrs)
		x1x2Ea <- sp::spTransform(x1x2, eaCrs)
		x2Ea <- subGeom(x2Ea, n=NULL)

		ext <- raster::extent(x1x2Ea)
		maxExt <- max(abs(ext@xmax - ext@xmin), abs(ext@ymax - ext@ymin))

		
	### tween each subpolygon of x1 to proximate subpolygons of x2
	##############################################################
	
		# for each x1 subgeometry
		numX1Subs <- countSubGeoms(x1Ea)
		for (countX1 in 1:numX1Subs) {
		
			if (verbose) omnibus::say('Processing subpolygon ', countX1, ' of ', numX1Subs, '...')

			# this subgeometry of x1
			x1SubEa <- subGeom(x1Ea, countX1)
	
			# x2 subgeometries that are proximate to this x1 subgeometry
			x2Prox <- which(!is.na(sp::over(x2Ea, x1SubEa)))
	
			# if no x2 subgeometries are proximate to this x1 subgeometry, tween to the x1 subgeometry centroid
			# NB add a small buffer to the point to make it a polygon
			if (length(x2Prox) == 0) {
			
				toSpEa <- rgeos::gCentroid(x1Ea)
				toSpEa <- rgeos::gBuffer(toSpEa, width=maxExt / 1E6)
				to <- geomToCoords(toSpEa)
				to <- as.data.frame(to)
				names(to) <- c('x', 'y')
				to$id <- 1
			
			# tween to x2 subgeometries that are proximate with this x1 subgeometry
			} else {

				# # add points to x2 subgeometry to make for better tweening
				# x2OnX1SubEa <- x2Ea[x2Prox]
				# x2OnX1SubEa <- rgeos::gUnaryUnion(x2OnX1SubEa)
				# numX2SubCoords <- nrow(geomToCoords(x2OnX1SubEa))
				# x2OnX1Sub <- sp::spTransform(x2OnX1SubEa, crs)
				# x2OnX1SubLength <- geosphere::lengthLine(x2OnX1Sub)
				# diff <- numX1Coords - numX2SubCoords
				# if (diff > 0) {
					# x2OnX1Sub <- geosphere::makePoly(x2OnX1Sub, x2OnX1SubLength / diff, sp=TRUE)
					# x2OnX1SubEa <- sp::spTransform(x2OnX1Sub, eaCrs)
				# }
			
				# get x2 subgeometries proximate to this x1 subgeometry and convert to data frame with id
				x2SubSpEa <- x2Ea[x2Prox[1]]
				
				to <- geomToCoords(x2SubSpEa)
				to <- as.data.frame(to)
				names(to) <- c('x', 'y')
				to$id <- 1
				
				if (length(x2Prox) > 1) {
					for (i in 2:length(x2Prox)) {
					
						this <- rbind(x2SubSpEa, x2Ea[x2Prox[i]])
						this <- geomToCoords(this)
						this <- as.data.frame(this)
						names(this) <- c('x', 'y')
						this$id <- i
						
						to <- rbind(to, this)
					
					}
				
				}
				
			}

			# coordinates of starting subgeometry of x1
			from <- geomToCoords(x1SubEa)
			from <- as.data.frame(from)
			names(from) <- c('x', 'y')
			from$id <- 1

			### interpolate by tweening
			tween <- transformr::tween_polygon(from, to, ease=method, nframes=delta, id=NULL, match = FALSE)
			# tween2 <- transformr::tween_polygon(from, to, ease=method, nframes=delta, id=NULL, match = FALSE) %>% keep_state(10)
			names(tween) <- c('x', 'y', 'id1', 'id2', 'phase', 'frame')
			
			whichFrame <- max(1, round(between * delta))
			thisTween <- tween[tween$frame == whichFrame, , drop=FALSE]
			
			### convert to spatial object(s)
			ids <- sort(unique(thisTween$id1))
			thisOut <- coordsToPoly(thisTween[thisTween$id1 == ids[1], c('x', 'y')], eaCrs)
			thisOut <- rgeos::gSimplify(thisOut, tol=0)
			
			if (length(ids) > 1) {
			
				ids <- sort(unique(thisTween$id1))
				uniqueIds <- seq_along(ids)
				for (countId in 2:length(ids)) {
				
					id <- ids[countId]
					thisThisOut <- coordsToPoly(thisTween[thisTween$id1 == id, c('x', 'y')], eaCrs, id=countId)
					thisThisOut <- rgeos::gSimplify(thisThisOut, tol=0)
					thisOut <- rgeos::gUnion(thisOut, thisThisOut)
				
				}
			
			}
			
			# remember
			interpSpEa <- if (exists('interpSpEa', inherits=FALSE)) {
				rgeos::gUnion(interpSpEa, thisOut)
			} else {
				thisOut
			}
	
		} # next subgeometry of x1
		
	### post-processing
	interpSpEa <- rgeos::gSimplify(interpSpEa, tol=0)
	interpSpEa <- rgeos::gIntersection(interpSpEa, x1Ea)
	interpSpEa <- rgeos::gUnion(interpSpEa, x2Ea)
	interpSp <- sp::spTransform(interpSpEa, crs)
	
	interpSp

}
