
library(rgeos)
library(geosphere)
library(enmSdm)
library(sp)
source('C:/Ecology/Drive/R/enmSdm/working/subGeomFromGeom.r')
source('C:/Ecology/Drive/R/enmSdm/working/coordsToPoly.r')
source('C:/Ecology/Drive/R/enmSdm/working/geomToCoords.r')
source('C:/Ecology/Drive/R/enmSdm/working/countSubGeoms.r')

load('C:/Ecology/Drive/R/enmSdm/working/x1.rda')
load('C:/Ecology/Drive/R/enmSdm/working/x2.rda')

x1tox2 <- 0.3

plot(x1)
plot(x2, border='blue', add=TRUE)
	
# plot(commonSp, add=TRUE, col='gray', border=NA)			

crs <- proj4string(x1)
crs <- CRS(crs)

eaCrs <- getCRS('albersNA', TRUE)

# plot(segLine, add=T, lty='dotted')



# plot(allEnds, add=T, lty='dashed')
# plot(segEndsSp, add=T, lty='dashed', col='green')

# points(interpPoint, col='red')

# plot(interpSp, add=T, border='purple', lwd=3)

plot(notCommonSp, add=T, border='orange', lwd=3)


