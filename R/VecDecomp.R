#' Vector decomposition and remainder fields
#'
#' This function calculates the vector and remainder fields.
#' @param surface matrix output from QPGlobal.
#' @param equations a two-element list with equations as strings, one for each equation.
#' @param mesh two-element vector with respective x and y mesh sizes.
#' @param x.lim two-element vector with respective minimum and maximum x values.
#' @param y.lim two-element vector with respective minimum and maximum y values.
#' @param remainder if \code{TRUE}, calculates the remainder field, default is \code{FALSE}.
#' @keywords vector field, remainder field
#' @export
#' @examples
#' x.limits <- c(-10,10)
#' y.limits <- c(-10,10)
#' length.xy <- c(15,15)
#' mesh.xy <- c(20,20)
#' x <- seq(x.limits[1], x.limits[2], length.out = length.xy[1]) + 0.5
#' y <- seq(y.limits[1], y.limits[2], length.out = length.xy[2]) + 0.5
#' eqns <- list("(y^2)/(x)" , "(x^2)/(y)")
#' f <- function(x, y) { x <- (y^2)/(x) ; y <- (x^2)/(y) }
#' z <- outer(x, y, f)
#' VecDecomp(z)
#' VecDecomp(z,eqns,mesh.xy,x.limits,y.limits,remainder=T)

VecDecomp <- function(surface,equations,x.lim,y.lim,remainder=FALSE){
	qpr <- nrow(surface)
	qpc <- ncol(surface)

	r.dx <- (surface[,qpc] - surface[,(qpc-1)])
	l.dx <- (surface[,2] - surface[,1])
	int.dx <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in (qpc-1):2){
		int.dx[,i] <- (surface[,(i+1)]-surface[,(i-1)])/2
	}
	int.dx[,1] <- l.dx
	int.dx[,qpc] <- r.dx
	dx <- -1*int.dx

	t.dy <- (surface[1,] - surface[2,])
	b.dy <- (surface[(qpr-1),] - surface[qpr,])
	int.dy <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in 2:(qpr-1)){
		int.dy[i,] <- (surface[(i-1),]-surface[(i+1),])/2
	}
	int.dy[1,] <- t.dy
	int.dy[qpr,] <- b.dy
	dy <- -1*int.dy

	if(remainder == FALSE){
		list(dx,dy)
		} else {
			x.val <- seq(min(x.lim),max(x.lim),length.out=qpc)
			y.val <- seq(min(y.lim),max(y.lim),length.out=qpr)
			n.eq <- length(equations)
			z.list <- vector(mode="list" , length = n.eq)
			for(i in 1:n.eq){
			z.list[[i]] <- outer(x.val,y.val,function(x,y){eval(parse(text=equations[[i]]))})
			}
			rx <- z.list[[1]] + dx
			ry <- z.list[[2]] + dy
			list(dx,dy,rx,ry)
		}
}