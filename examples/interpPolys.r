
source('C:/Ecology/Drive/R/birdsEye/R/interpPolysByBuffer.r')
source('C:/Ecology/Drive/R/birdsEye/R/interpPolysByTween.r')

# create "x1": has two sub-polygons
x <- c(44.02, 43.5, 42.61, 42.18, 42, 42.41, 42.75, 41.75, 41.49,
43.61,46.02, 46.5, 47.5, 47.39, 48.64, 49.05, 48.46, 48.18, 47.54, 46.73, 45.80, 45.59)
y <- c(-18.83, -18.67, -18.87, -19.67, -20.65, -21.64, -23.08, -24.9,
-26.51, -27.09, -26.74, -25.6, -25.14, -26.44, -26.46, -24.96, -23.63,
 -22.72, -23.36, -22.29, -21.45, -20.69)
xy1a <- cbind(x, y)

x <- c(40.61, 40.07, 40.23, 41.38, 41.38)
y <- c(-20.51, -20.49, -21.11, -21.55, -21.01)
xy1b <- cbind(x, y)

x1a <- coordsToPoly(xy1a, enmSdm::getCRS('wgs84'))
x1b <- coordsToPoly(xy1b, enmSdm::getCRS('wgs84'))
x1 <- rgeos::gUnion(x1a, x1b)

# create "x2"
x <- c(44.53, 44.18, 44.00, 42.93, 42.29, 42.71, 43.43, 47.15, 48.08,
 45.94,45.36, 45.76, 46.97, 46.87, 45.94, 45.97, 45.08, 44.50, 44.58)
y <- c(-24.27, -23.68, -22.86, -21.88, -20.56, -19.31, -20.36, -20.53,
 -20.93,-21.81, -21.64, -22.90, -23.44, -24.08, -24.76, -25.95, -25.88, -25.61, -24.46)
xy2 <- cbind(x, y)

x2 <- coordsToPoly(xy2, enmSdm::getCRS('wgs84'))

eaCrs <- enmSdm::getCRS('albersNA')

interBuff <- interpPolysByBuffer(
	x1, x2, eaCrs=eaCrs, between = 0.4, delta=10000
)

interTween <- interpPolysByTween(
	x1, x2, eaCrs='laea', between = 0.4, delta=100
)

plot(x1, col='gray90')
plot(x2, add=TRUE)
plot(interBuff, border='red', add=TRUE)
plot(interTween, border='green', add=TRUE)
legend('bottomleft',
	legend=c('x1', 'x2', 'by buffer', 'by tween'),
	fill=c('gray90', NA, NA, NA),
	border=c('black', 'black', 'red', 'green'),
	bty='n'
)
