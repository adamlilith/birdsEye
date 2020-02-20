
library(rgeos)
library(geosphere)
library(omnibus)
library(enmSdm)
library(sp)
library(scales)
source('C:/Ecology/Drive/R/enmSdm/working/subGeomFromGeom.r')
source('C:/Ecology/Drive/R/enmSdm/working/coordsToPoly.r')
source('C:/Ecology/Drive/R/enmSdm/working/geomToCoords.r')
source('C:/Ecology/Drive/R/enmSdm/working/countSubGeoms.r')
source('C:/Ecology/Drive/R/enmSdm/working/interpolateSpatialShapesByBuffer.r')
source('C:/Ecology/Drive/R/enmSdm/working/interpolateSpatialShapes.r')

load('C:/Ecology/Drive/R/enmSdm/working/x1.rda')
load('C:/Ecology/Drive/R/enmSdm/working/x2.rda')
# can <- readRDS9('C:/Ecology/!Scratch/gadm36_CAN_0_sp.rds')

# eaCrs <- getCRS('albersNA')
eaCrs <- getCRS('mollweide')
between <- 0.3
delta <- 100000
method <- 'grow'

plot(x1, col=alpha('blue', 0.1))
plot(x2, border='orange', lwd=4, add=TRUE)

# x1 <- makePoly(x1, 500000)
# x2 <- makePoly(x2, 500000)
interBuff <- interpolateSpatialShapesByBuffer(x1, x2, eaCrs, delta=delta, between=between, quadsegs=5, verbose=2)
plot(interBuff, border='purple', lwd=2, add=TRUE)

inter <- interpolateSpatialShapes(x1, x2, eaCrs, between = between, delta=delta, verbose=TRUE)
plot(inter, border='green', lwd=2, add=TRUE)



