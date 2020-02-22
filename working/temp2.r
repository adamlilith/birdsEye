# source('C:/Ecology/Drive/R/birdsEye/working/temp2.r')
# source('D:/Ecology/Drive/R/birdsEye/working/temp2.r')
	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	gc()
	options(stringsAsFactors=FALSE)

	library(sp)
	library(transformr)
	library(scales)
	
	library(omnibus) # Adam's custom library on GitHub: adamlilith/omnibus
	library(enmSdm) # Adam's custom library on GitHub: adamlilith/enmSdm

	drive <- 'C:'
	# drive <- 'D:'
	source(paste0(drive, '/Ecology/Drive/R/birdsEye/R/coordsToPoly.r'))
	source(paste0(drive, '/Ecology/Drive/R/birdsEye/R/countSubGeoms.r'))
	source(paste0(drive, '/ecology/Drive/R/birdsEye/R/geomToCoords.r'))
	source(paste0(drive, '/ecology/Drive/R/birdsEye/R/subGeom.r'))
	source(paste0(drive, '/ecology/Drive/R/birdsEye/R/ellipsoid.r'))
	source(paste0(drive, '/ecology/Drive/R/enmSdm/R/makeCRS.r'))

	datesCrosswalk <- read.csv(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke & Peltier/Dates Crosswalk.csv'))
	
	
	startIce <- shapefile(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke & Peltier/IceMaps_1ka_All/07000icell'))
	endIce <- shapefile(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke & Peltier/IceMaps_1ka_All/08000icell'))
	

	x1 <- endIce
	x2 <- startIce
	# method <- 'cubic-in-out' # 'linear'
	method <- 'cubic-in-out' # 'linear'
	steps <- 100
	between <- 0.5
	# between <- 1
	verbose = TRUE

	# plot(endIce)
	# plot(startIce, border='orange', lwd=2, add=TRUE)

	# plot(ice, col=alpha('blue', 0.2), add=TRUE)

	### create equal-area CRS with centroid of polygons as its center
	#################################################################

		# x1 <- subGeom(x1, 11)
		# x2 <- subGeom(x2, 1:countSubGeoms(x2))
	
		x1x2 <- rgeos::gUnion(x1, x2)
		centSp <- geosphere::centroid(x1x2)
		
		crs <- sp::proj4string(x1)
		ell <- ellipsoid(crs)
		dat <- datum(crs)
		
		long0 <- sp::coordinates(centSp)[1, 1]
		lat0 <- sp::coordinates(centSp)[1, 2]
		eaCrs <- makeCRS('laea', long0=long0, lat0=lat0, datum=dat, ell=ell)
		
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
		for (countX1 in 1:countSubGeoms(x1Ea)) {
		# for (countX1 in 11) {

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
			tween <- transformr::tween_polygon(from, to, ease=method, nframes=steps, id=NULL, match = FALSE)
			# tween2 <- transformr::tween_polygon(from, to, ease=method, nframes=steps, id=NULL, match = FALSE) %>% keep_state(10)
			names(tween) <- c('x', 'y', 'id1', 'id2', 'phase', 'frame')
			
			whichFrame <- max(1, round(between * steps))
			thisTween <- tween[tween$frame == whichFrame, , drop=FALSE]
			
			### convert to spatial object(s)
			ids <- sort(unique(thisTween$id1))
			thisOut <- coordsToPoly(thisTween[thisTween$id1 == ids[1], c('x', 'y')], eaCrs)
			thisOut <- rgeos::gSimplify(thisOut, tol=0)
			
			if (length(ids) > 1) {
			
				ids <- sort(unique(tween$id1))
				uniqueIds <- seq_along(ids)
				for (countId in 2:length(ids)) {
				
					id <- ids[countId]
					thisThisOut <- coordsToPoly(tween[tween$id1 == id, c('x', 'y')], eaCrs, id=countId)
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
		
		interpSpEa <- rgeos::gSimplify(interpSpEa, tol=0)
		interpSpEa <- rgeos::gIntersection(interpSpEa, x1Ea)
		interpSpEa <- rgeos::gUnion(interpSpEa, x2Ea)
		interpSpEa <- spatialEco::remove.holes(interpSpEa)
		interpSp <- sp::spTransform(interpSpEa, crs)

	plot(x1Ea)
	plot(x2Ea, border='orange', lwd=2, add=TRUE)

		plot(interpSpEa, col=alpha('blue', 0.2), add=TRUE)

	
