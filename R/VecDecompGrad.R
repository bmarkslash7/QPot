#' Vector decomposition and remainder fields
#'
#' This function calculates the gradient field.
#' @param surface matrix output from QPGlobal.
#' @keywords vector field decompoosition, gradient field
#' 
########################################################################
# THIS EXAMPLE FAILS WITH
# > VecDecomp(z)
# Error in seq(min(x.lim), max(x.lim), length.out = qpc) : 
#  argument "x.lim" is missing, with no default
########################################################################
# @examples
# x.limits <- c(-10,10)
# y.limits <- c(-10,10)
# length.xy <- c(15,15)
# mesh.xy <- c(20,20)
# x <- seq(x.limits[1], x.limits[2], length.out = length.xy[1]) + 0.5
# y <- seq(y.limits[1], y.limits[2], length.out = length.xy[2]) + 0.5
# eqns <- list("(y^2)/(x)" , "(x^2)/(y)")
# f <- function(x, y) { x <- (y^2)/(x) ; y <- (x^2)/(y) }
# z <- outer(x, y, f)
# VecDecomp(z)
# VecDecompGrad(z)

VecDecompGrad <- function(surface){
	qpr <- nrow(surface)
	qpc <- ncol(surface)

	# column derivative
	r.dc <- (surface[,qpc] - surface[,(qpc-1)])
	l.dc <- (surface[,2] - surface[,1])
	int.dc <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in (qpc-1):2){
		int.dc[,i] <- (surface[,(i+1)]-surface[,(i-1)])/2
	}
	int.dc[,1] <- l.dc
	int.dc[,qpc] <- r.dc
	dc <- -int.dc

	# row derivative
	t.dr <- (surface[1,] - surface[2,])
	b.dr <- (surface[(qpr-1),] - surface[qpr,])
	int.dr <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in 2:(qpr-1)){
		int.dr[i,] <- (surface[(i-1),]-surface[(i+1),])/2
	}
	int.dr[1,] <- t.dr
	int.dr[qpr,] <- b.dr
	dr <- int.dr

	array(data=c(dr,dc),dim=c(qpr,qpc,2))
}
