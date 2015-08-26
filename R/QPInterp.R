#' Quasi-potential interpolation
#'
#' This function estimates the quasi-potential value for any x- and y-values
#' @param X the value of the x to interpolate
#' @param Y the value of the y to interpolate
#' @param x.bound two-element vector with respective minimum and maximum x values.
#' @param y.bound two-element vector with respective minimum and maximum y values.
#' @param surface the surface to interpolated, from \code{\link{QPGlobal}}
#' @keywords interpolation
#' @details this function uses simple bilinear interpolation for estimattion of any x- and y-value.
########################################################################
# COMMENTED OUT FOR PASSING devtools::check() - Need a srface from QPGlobal
########################################################################
#' @examples
#' QPInterp(X=3,Y=3,x.bound=c(-0.5,20),y.bound=c(-0.5,20),surface=e1.global)

	QPInterp <- function(X, Y, x.bound, y.bound, surface) {
		mat.row.min <- 1
		mat.row.max <- nrow(surface)
		mat.col.min <- 1
		mat.col.max <- ncol(surface)


		X.range <- max(x.bound) - min(y.bound)
		Y.range <- max(y.bound) - min(y.bound)
		x.index <- ((X - min(x.bound))/X.range)*mat.row.max
		y.index <- ((Y - min(y.bound))/Y.range)*mat.col.max	

		x.1 <- floor(x.index)
		x.2 <- ceiling(x.index)
		x <- x.index
	
		y.1 <- floor(y.index)
		y.2 <- ceiling(y.index)
		y <- y.index

		Q.11 <- surface[x.1 , y.1]
		Q.12 <- surface[x.2 , y.1]
		Q.21 <- surface[x.1 , y.2]
		Q.22 <- surface[x.2 , y.2]	

		if(x.1 == x.2 && y.1 != y.2) {
			print("X SAME!")			
		val <- (1/(y.2-y.1)) * ((Q.11*(y.2-y)) + (Q.22*(y-y.1)))
		return(val)
		}
		
 		if (y.1 == y.2 & x.1 != x.2) {
			print("Y SAME!")
		val <- (1/(x.2-x.1)) * ((Q.11*(x.2-x)) + (Q.22*(x-x.1)))
		return(val)		
		}

		if (x.1 == x.2 & y.1 == y.2) {
		return(surface[x,y])
		}

		if (x.1 != x.2 & y.1 != y.2) {
		val <- (1/((x.2-x.1)*(y.2-y.1))) * ((Q.11*(x.2-x)*(y.2-y)) + (Q.21*(x-x.1)*(y.2-y)) + (Q.12*(x.2-x)*(y-y.1)) + (Q.22*(x-x.1)*(y-y.1)))
		return(val)
		}
	}