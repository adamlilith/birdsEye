#' Create a custom "PROJ.4" string for a datum and projection
#'
#' This function returns a custom "PROJ.4" string for a particular datum and possibly projection. Currently only the Lambert azimuthal equal-area projection and Albers conic equal-area projections are supported.
#' @param x Character. Name of PROJ.4 string to return. Possible values include:
#'	\itemize{
#'	\item \code{albers} CRS for the Albers equal-area projection: \code{+proj=aea +lat_1=lat1 +lat_2=lat2 +lat_0=lat0 +lon_0=long0 +x_0=0 +y_0=0 +ellps=ell +datum=datum +units=m}, where \code{ell} is the ellipsoid; \code{datum} the datum; the pair \code{long0} and \code{lat0} the longitude and latitude of the center point of the projection (in degrees); and \code{lat1} and \code{lat2} the first and second standard parallels. This projection is good for mid-latitudes.
#'	\item \code{lcc} CRS for the Lambert conformal conic projection: \code{'+proj=lcc +lat_1=lat1 +lat_2=lat2 +lat_0=lat0 +lon_0=long0 +x_0=0 +y_0=0 +ellps=ell +datum=dat +units=m'}, where \code{ell} is the ellipsoid; the pair \code{long0} and \code{lat0} are the longitude and latitude of the center point of the projection (in degrees); and and \code{lat1} and \code{lat2} are the first and second standard parallels.
#'	\item \code{laea} CRS for the Lambert azimuthal equal-area projection: \code{+proj=laea +lat_0=lat0 +lon_0=long0 +x_0=0 +y_0=0 +ellps= +datum=dat +units=m}, where \code{ell} is the ellipsoid; and the pair \code{long0} and \code{lat0} are the longitude and latitude of the center point of the projection (in degrees). Defaults: \code{ell = 'GRS80'}. This projection is good for high latitudes/polar regions.
#'	\item \code{mollweide} CRS for the Mollweide equal-area projection: \code{+proj=moll +lon_0=long0 +x_0=0 +y_0=0 +ellps=ell +datum=dat +units=m}, where \code{ell} is the ellipsoid; and \code{long0} the center meridian of the projection (in degrees). This projection is good for world maps.
#'	\item \code{naturalEarth} or \code{ne} CRS for the Natural Earth "compromise" projection: \code{+proj=natearth +lon_0=long0 +ellps=ell +datum=dat +units=m}, where \code{long0} is the central meridian, \code{ell} is the ellipsoid; and \code{long0} the centeral meridian of the projection (in degrees). This projection is good for world maps.
#' 	\item \code{oae} or \code{aeqd} CRS for the oblique azimuthal equidistant projection: \code{'+proj=aeqd +lat_0=lat0 +lon_0=long0 +ellps=ell +datum=dat'} where \code{long0} and \code{lat0} are the center point of the projection (in degrees). This projection is good for any region of the world.
#' 	\item \code{winkelTriple} or \code{wt} CRS for the Winkel Triple: \code{'+proj=wintri +lon_0=long0'} where \code{long0} is the central meridian of the projection (in degrees). This projection is good for the world.
#' }
#' @param long0,lat0 Numeric, specify central point or meridian of the projection. Units are degrees. Both or just one may be needed for a given projection.
#' @param long1,long2,lat1,lat2 Numeric, specify additional standard parallels or meridians, depending on the projection. These are not needed for all projections.
#' @param dat Character, name of datum.  If \code{NULL}, then default values are used. Typically, if you specify the datum then the ellipsoid (\code{ell}) need not be specified since a particular ellipsoid is used for each datum.
#' @param ell Character, name of ellipse. If \code{NULL}, then default values are used.  
#' @param asCRS Logical. If \code{TRUE} then return object is of class \code{CRS}. If \code{FALSE} (default) then returned object is of class \code{character}.
#' @details The function does its best to fill in information. So if the datum is any of \code{WGS84}, \code{NAD83}, or \code{NAD27}, the the string \code{+no_defs} is attached to the PROJ.4 string. Definitions are also attached if they are specified in \code{rgdal::projInfo('datum')}.
#' @return Object of class \code{CRS} or \code{character}.
#' @seealso \code{\link[enmSdm]{getCRS}}
#' @examples
#' # Lambert azimuthal equal-area
#' makeCRS('laea', long0=10, lat0=52) # for Europe
#' makeCRS('laea', long0=-96, lat0=37.5) # for North America
#' # Albers equal-area conic for North America
#' makeCRS('albers', long0=-96, lat0=37.5, lat1=29.5, lat2=45.5)
#' @export
makeCRS <- function(
	x,
	long0 = NULL,
	lat0 = NULL,
	long1 = NULL,
	lat1 = NULL,
	long2 = NULL,
	lat2 = NULL,
	dat = 'WGS84',
	ell = 'WGS84',
	asCRS = FALSE
) {

	datsDefs <- rgdal::getInfo('datum')

	# default datum/ellipsoid
	if (is.null(dat)) {

		dat <- 'WGS84'
		if (is.null(ell)) ell <- 'WGS84'
	
	# get ellipsoid and definitions from datum
	} else if (!is.null(dat) & is.null(ell)) {
	
		thisDat <- which(tolower(dat) == tolower(datsDefs$name))
		if (length(thisDat) == 1) {
			ell <- datsDefs$ellipse[thisDat]
			defs <- paste0(' +', datsDefs$definition[thisDat])
		} else {
			stop('Can find no matching ellipsoid for this datum.')
		}
	
	# both datum and ellipsoid specified
	} else if (!is.null(dat) & !is.null(ell)) {
	
		thisDat <- which(tolower(dat) == tolower(datsDefs$name))
		defs <- if (length(thisDat) == 1) {
			paste0(' +', datsDefs$definition[thisDat])
		} else {
			''
		}
		
	}
	
	noDefs <- if (dat %in% c('WGS84', 'NAD83', 'NAD27')) {
		' +no_defs'
	} else {
		''
	}
	
	x <- tolower(x)

	out <- if (x == 'albers') {
		paste0('+proj=aea +lat_1=', lat1,' +lat_2=', lat2, '+lat_0=', lat0, ' +lon_0=', long0, ' +x_0=0 +y_0=0 +ellps=', ell, ' +datum=', dat, ' +units=m', noDefs, defs)
	} else if (x == 'laea') {
		paste0('+proj=laea +lat_0=', lat0, ' +lon_0=', long0, ' +x_0=4321000 +y_0=3210000 +ellps=', ell, ' +datum=', dat, noDefs, defs)
	} else if (x == 'lcc') {
		paste0('+proj=lcc +lat_1=', lat1, ' +lat_2=', lat2, ' +lat_0=', lat0, ' +lon_0=', long0, ' +x_0=0 +y_0=0 +ellps=', ell, ' +datum=', dat, ' +units=m', noDefs, defs)
	} else if (x == 'mollweide') {
		paste0('+proj=moll +lon_0=', long0, ' +x_0=0 +y_0=0 +ellps=', ell, ' +datum=', dat, ' +units=m', noDefs, defs)
	} else if (x %in% c('ne', 'naturalearth')) {
		paste0('+proj=natearth +lon_0=, ', long0, ' +ellps=ell +datum=dat +units=m', noDefs, defs)
	} else if (x %in% c('oae', 'aeqd')) {
		paste0('+proj=aeqd +lat_0=', lat0, ' +lon_0=', long0, ' +ellps=', ell, ' +datum=', dat, noDefs, defs)
	} else if (x %in% c('wt', 'winkeltriple')) {
		paste0('+proj=wintri +lon_0=', long0, ' +ellps=', ell, ' +datum=', dat, noDefs, defs)
	} else {
		NA
	}

	if (asCRS) out <- sp::CRS(out)
	out

}
