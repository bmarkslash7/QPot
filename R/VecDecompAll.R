#' Vector decomposition and remainder fields
#'
#' This function calculates the vector, gradient, and remainder fields.
#' @param surface matrix output from QPGlobal.
#' @param equations a two-element list with equations as strings, one for each equation.
#' @param x.bound two-element vector with respective minimum and maximum x values.
#' @param y.bound two-element vector with respective minimum and maximum y values.
#' @keywords vector field, remainder field
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
# VecDecomp(z,eqns,x.bound,y.bound)

VecDecompAll <- function(surface,equations,x.bound,y.bound){
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

	# vector and remainder fields
	x.val <- seq(min(x.bound),max(x.bound),length.out=qpc)
	y.val <- seq(min(y.bound),max(y.bound),length.out=qpr)
	n.eq <- length(equations)
	z.list <- vector(mode="list" , length = n.eq)
	for(i in 1:n.eq){
		z.list[[i]] <- 
		outer(x.val,y.val,function(x,y){eval(parse(text=equations[i]))})
		}

		# verctor field
		vx <- z.list[[1]]
		vy <- z.list[[2]]

		#remainder field
		rx <- z.list[[1]]+dr
		ry <- z.list[[2]]+dc

	array(data=c(vx,vy,dr,dc,rx,ry),dim=c(qpr,qpc,6))
}
