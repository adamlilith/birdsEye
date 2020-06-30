#' Convert coordinates in decimal degrees to degrees-minutes-seconds format
#'
#' This function converts decimal degrees into degrees-minutes-seconds format.
#' @param x Numeric or numeric vector.
#' @return Either a numeric vector with four named elements or a matrix with four columns (degrees, minutes, seconds, hemisphere). "Hemisphere" is indicated by 1 (north or east) or -1 (south or west).
#' @examples
#' @export
decimalToDms <- function(x) {

	out <- matrix(NA, nrow=length(x), ncol=4)
	colnames(out) <- c('degrees', 'minutes', 'seconds', 'hemisphere')
	
	for (i in seq_along(x)) {
	
		
		out[i, 'degrees'] <- abs(trunc(x[i]))

		afterDec <- abs(x[i] - trunc(x[i]))
		out[i, 'minutes'] <- trunc(60 * afterDec)
		out[i, 'seconds'] <- 3600 * (afterDec - trunc(afterDec * 100) / 100)

		out[i, 'hemisphere'] <- if (x[i] < 0) { -1 } else { 1 }
	
	}

	out
	
}
