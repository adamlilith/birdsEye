
library(rgeos)
library(geosphere)
library(omnibus)
library(enmSdm)
library(sp)
source('C:/Ecology/Drive/R/enmSdm/working/subPolyFromPoly.r')
source('C:/Ecology/Drive/R/enmSdm/working/coordsToPoly.r')
source('C:/Ecology/Drive/R/enmSdm/working/polyToCoords.r')
source('C:/Ecology/Drive/R/enmSdm/working/countSubPolys.r')
source('C:/Ecology/Drive/R/enmSdm/working/interpolateSpatialShapes.r')

load('C:/Ecology/Drive/R/enmSdm/working/x1.rda')
load('C:/Ecology/Drive/R/enmSdm/working/x2.rda')


x1tox2 <- 0.3
tolDist <- 10000

inter <- interpolateSpatialShapes(
	x1=x1,
	x2=x2,
	eaCrs=getCRS('albersNA', TRUE),
	x1tox2 = x1tox2,
	tolDist = tolDist,
	verbose = TRUE
)

plot(x1)
plot(x2, border='blue', add=TRUE)
plot(inter, border='purple', add=TRUE)


