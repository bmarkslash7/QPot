#' Plotting function for vector decomposition and remainder fields
#'
#' This function calculates the vector and remainder fields.
#' @param field list output from \code{\link{VecDecomp}}.
#' @param density two-element vector respectively specifying the number of arrows in the x and y directions.
#' @param x.lim two-element vector for the x variable's minimum and maximum values
#' @param y.lim two-element vector for the y variable's minimum and maximum values
#' @param arrow.type sets the type of line segments plotted. If set to "proportional" the length of the line segments reflects the magnitude of the derivative. If set to "equal" the line segments take equal lengths, simply reflecting the gradient of the derivative(s). Defaults to "equal".

#' @keywords vector field plot, remainder field plot
#' @export
#' @examples
#' x.limits <- c(-10,10)
#' y.limits <- c(-10,10)
#' length.xy <- c(15,15)
#' x <- seq(x.limits[1], x.limits[2], length.out = length.xy[1]) + 0.5
#' y <- seq(y.limits[1], y.limits[2], length.out = length.xy[2]) + 0.5
#' eqns <- list("(y^2)/(x)" , "(x^2)/(y)")
#' f <- function(x, y) { x <- (y^2)/(x) ; y <- (x^2)/(y) }
#' z <- outer(x, y, f)
#' VecDecomp(z)
#' VecDecomp(z,eqns,mesh.xy,x.limits,y.limits,remainder=T)

VecDecompPlot <- function(field, density, x.lim, y.lim ){
	sub.dx.x <- seq(1, ncol(field[[1]]), length.out=density[1])
	sub.dx.y <- seq(1, ncol(field[[1]]), length.out=density[2])
	sub.dy.x <- seq(1, ncol(field[[2]]), length.out=density[1])
	sub.dy.y <- seq(1, ncol(field[[2]]), length.out=density[2])

	dx.sub <- field[[1]][sub.dx.x, sub.dx.y]
	dy.sub <- field[[2]][sub.dy.x, sub.dy.y]

	dx.rel <- (dx.sub/max(((dx.sub^2)+(dy.sub^2))^0.5, na.rm = T))
	dy.rel <- (dy.sub/max(((dx.sub^2)+(dy.sub^2))^0.5, na.rm = T))

	# dx.even <- -1*dx.sub/((dx.sub^2)+(dy.sub^2))^0.5
	# dy.even <- -1*dy.sub/(((dx.sub^2)+(dy.sub^2))^0.5)

	rowmin <- min(which(dx.sub != 0 , arr.ind = T)[,1])
	rowmax <- max(which(dx.sub != 0 , arr.ind = T)[,1])
	colmin <- min(which(dx.sub != 0 , arr.ind = T)[,2])
	colmax <- max(which(dx.sub != 0 , arr.ind = T)[,2])

	x.range <- max(x.lim)-min(x.lim)
	y.range <- max(y.lim)-min(y.lim)

	x.pretty <- pretty(c(((colmin-min(x.lim))/x.range)*density[1],((colmin-min(x.lim))/x.range)*density[1]))
	x.pretty.at <- ((x.pretty-min(x.lim))/x.range)*density[1]
	y.pretty <- pretty(c(((rowmin-min(y.lim))/y.range)*density[2],((rowmin-min(y.lim))/y.range)*density[2]))
	y.pretty.at <- ((y.pretty-min(y.lim))/y.range)*density[2]

	row.min <- c(rowmin,rowmax)
	col.min <- c(colmin,colmax)

	qpr <- nrow(dx.rel)
	qpc <- ncol(dx.rel)

	plot(0 , type = "n" , xlim = row.min , ylim = col.min , las = 1 , xlab = expression(italic(X)) , ylab = expression(italic(Y)) , xaxt = "n" , yaxt = "n")
for (i in 1:(qpr)){
		for (j in 1:qpc){
			x1 <- i + (dy.rel[i,j]/2)
			x2 <- i - (dy.rel[i,j]/2)
			y1 <- j - (dx.rel[i,j]/2)
			y2 <- j + (dx.rel[i,j]/2)
			arrows(x1,y1,x2,y2 , length = 0.03 , lwd = 1.1)
		}
	}
axis(1,at=x.pretty.at,labels=x.pretty)
axis(2,at=y.pretty.at,labels=y.pretty,las=1)
}