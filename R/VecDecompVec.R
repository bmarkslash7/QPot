#' Vector decomposition and remainder fields
#'
#' This function calculates the vector field.
#' @param x.num.steps COPY FROM QPotential
#' @param y.num.steps COPY FROM QPotential
#' @param x.bound two-element vector with respective minimum and maximum x values.
#' @param y.bound two-element vector with respective minimum and maximum y values.
#' @keywords vector field decompoosition, vector field
#' 
########################################################################
# THIS EXAMPLE FAILS WITH
# > VecDecomp(z)
# Error in seq(min(x.lim), max(x.lim), length.out = qpc) : 
#  argument "x.lim" is missing, with no default
########################################################################
# @examples
# x.donmain.boundaries <- c(-10,10)
# y.donmain.boundaries <- c(-10,10)
# length.xy <- c(15,15)
# mesh.xy <- c(20,20)
# x <- seq(x.limits[1], x.limits[2], length.out = length.xy[1]) + 0.5
# y <- seq(y.limits[1], y.limits[2], length.out = length.xy[2]) + 0.5
# eqns <- list("(y^2)/(x)" , "(x^2)/(y)")
# f <- function(x, y) { x <- (y^2)/(x) ; y <- (x^2)/(y) }
# z <- outer(x, y, f)
# VecDecompVec(z)

VecDecompVec <- function(x.num.steps,y.num.steps,equations,x.bound,y.bound){

	x.val <- seq(min(x.bound),max(x.bound),length.out=x.num.steps)
	y.val <- seq(min(y.bound),max(y.bound),length.out=y.num.steps)
	n.eq <- length(equations)
	z.list <- vector(mode="list" , length = n.eq)
	for(i in 1:n.eq){
		z.list[[i]] <- 
		outer(x.val,y.val,function(x,y){eval(parse(text=equations[i]))})
		}

	array(data=c(z.list[[1]],z.list[[2]]),dim=c(x.num.steps,y.num.steps,2))
}