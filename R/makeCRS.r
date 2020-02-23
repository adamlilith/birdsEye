#' Create a custom "proj4" string for a datum and projection
#'
#' This function returns a custom "proj4" string for a particular datum and possibly projection. Currently only the Lambert azimuthal equal-area projection and Albers conic equal-area projections are supported.
#' @param x Character. Name of proj4 string to return. Possible values include:
#'	\itemize{
#'	\item \code{albers} CRS for the Albers equal-area projection: \code{+proj=aea +lat_1=lat1 +lat_2=lat2 +lat_0=lat0 +lon_0=long0 +x_0=0 +y_0=0 +ellps=ell +datum=datum +units=m +no_defs}, where \code{ell} is the ellipsoid; \code{datum} the datum; the pair \code{long0} and \code{lat0} the longitude and latitude of the center point of the projection (in degrees); and \code{lat1} and \code{lat2} the first and second standard parallels. This projection is good for mid-latitudes.
#'	\item \code{lcc} CRS for the Lambert conformal conic projection: \code{'+proj=lcc +lat_1=lat1 +lat_2=lat2 +lat_0=lat0 +lon_0=long0 +x_0=0 +y_0=0 +ellps=ell +datum=dat +units=m no_defs'}, where \code{ell} is the ellipsoid; the pair \code{long0} and \code{lat0} are the longitude and latitude of the center point of the projection (in degrees); and and \code{lat1} and \code{lat2} are the first and second standard parallels.
#'	\item \code{laea} CRS for the Lambert azimuthal equal-area projection: \code{+proj=laea +lat_0=lat0 +lon_0=long0 +x_0=0 +y_0=0 +ellps= +datum=dat +units=m no_defs}, where \code{ell} is the ellipsoid; and the pair \code{long0} and \code{lat0} are the longitude and latitude of the center point of the projection (in degrees). Defaults: \code{ell = 'GRS80'}. This projection is good for high latitudes/polar regions.
#'	\item \code{mollweide} CRS for the Mollweide equal-area projection: \code{+proj=moll +lon_0=long0 +x_0=0 +y_0=0 +ellps=ell +datum=dat +units=m no_defs}, where \code{ell} is the ellipsoid; and \code{long0} the center meridian of the projection (in degrees). This projection is good for world maps.
#' 	\item \code{oae} or \code{aeqd} CRS for the oblique azimuthal equidistant projection: \code{'+proj=aeqd +lat_0=', lat0, ' +lon_0=', long0} where \code{long0} and \code{lat0} are the center point of the projection (in degrees).
#' }
#' @param long0,lat0 Numeric, specify central point or meridian of the projection. Units are degrees. Both or just one may be needed for a given projection.
#' @param long1,long2,lat1,lat2 Numeric, specify additional standard parallels or meridians, depending on the projection. These are not needed for all projections.
#' @param ell Character, name of ellipse. If \code{NULL}, then default values are used.
#' @param dat Character, name of datum.  If \code{NULL}, then default values are used.
#' @param asCRS Logical. If \code{TRUE} then return object is of class \code{CRS}. If \code{FALSE} (default) then returned object is of class \code{character}.
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
	ell = 'WGS84',
	dat = 'WGS8',
	asCRS = FALSE
) {

	x <- tolower(x)
	if (!is.null(ell)) ell <- toupper(ell)
	if (!is.null(datum)) dat <- toupper(dat)

	if (x == 'albers') {
		out <- paste0('+proj=aea +lat_1=', lat1,' +lat_2=', lat2, '+lat_0=', lat0, ' +lon_0=', long0, ' +x_0=0 +y_0=0 +ellps=', ell, ' +datum=', dat, ' +units=m +no_defs')
	} else if (x == 'laea') {
		out <- paste0('+proj=laea +lat_0=', lat0, ' +lon_0=', long0, ' +x_0=4321000 +y_0=3210000 +ellps=', ell, ' +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
	} else if (x == 'lcc') {
		out <- paste0('+proj=lcc +lat_1=', lat1, ' +lat_2=', lat2, ' +lat_0=', lat0, ' +lon_0=', long0, ' +x_0=0 +y_0=0 +ellps=', ell, ' +datum=', dat, ' +units=m no_defs')
	} else if (x == 'mollweide') {
		out <- paste0('+proj=moll +lon_0=', long0, ' +x_0=0 +y_0=0 +ellps=', ell, ' +datum=', dat, ' +units=m no_defs')
	} else if (x == 'oae' | x == 'aeqd') {
		out <- paste0('+proj=aeqd +lat_0=', lat0, ' +lon_0=', long0)
	} else {
		NA
	}

	if (asCRS) out <- sp::CRS(out)
	out

}
