# source('C:/Ecology/Drive/R/birdsEye/working/temp.r')
	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	gc()
	options(stringsAsFactors=FALSE)

	library(sp)
	library(raster)
	library(dismo)
	library(rgeos)
	
	library(omnibus) # Adam's custom library on GitHub: adamlilith/omnibus
	library(enmSdm) # Adam's custom library on GitHub: adamlilith/enmSdm

	source('C:/Ecology/Drive/R/birdsEye/R/coordsToPoly.r')
	source('C:/Ecology/Drive/R/birdsEye/R/countSubGeoms.r')
	source('C:/ecology/Drive/R/birdsEye/R/geomToCoords.r')
	source('C:/ecology/Drive/R/birdsEye/R/subGeomFromGeom.r')

	drive <- 'C:'
	datesCrosswalk <- read.csv(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke & Peltier/Dates Crosswalk.csv'))
	
	
	startIce <- shapefile(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke & Peltier/IceMaps_1ka_All/07000icell'))
	endIce <- shapefile(paste0(drive, '/Ecology/Drive/Data/Paleoglaciers of North America - Dyke & Peltier/IceMaps_1ka_All/08000icell'))
	
	# between <- 0.7
	between <- 1

	x1 <- endIce
	x2 <- startIce
	eaCrs <- getCRS('albersNA')
	method <- 'grow'
	delta = 30000
	quadsegs = 5
	verbose = TRUE

	plot(endIce)
	plot(startIce, border='orange', lwd=2, add=TRUE)

source('C:/Ecology/Drive/R/birdsEye/R/interpPolysByBuffer.r')
	ice <- interpPolysByBuffer(
		x1=endIce,
		x2=startIce,
		eaCrs=getCRS('mollweide'),
		between=between,
		delta=20000
	)
	
	plot(ice, col=alpha('blue', 0.2), add=TRUE)
	
	