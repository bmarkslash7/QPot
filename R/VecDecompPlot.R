#' Plotting function for vector decomposition and remainder fields
#'
#' This function calculates the vector and remainder fields.
#' @param field list output from \code{\link{VecDecompAll}}.
#' @param dens two-element vector respectively specifying the number of respective arrows in the x and y directions.
#' @param x.bound two-element vector for the x domain boundries used for the quasi-potential simulation.
#' @param y.bound two-element vector for the y domain boundries used for the quasi-potential simulation.
#' @param x.lim I DO NOT KNOW BUT NOW devtools::check() IS OK
#' @param y.lim I DO NOT KNOW BUT NOW devtools::check() IS OK
#' @param arrow.type sets the type of line segments plotted. If set to "proportional" the length of the line segments reflects the magnitude of the derivative. If set to "equal" the line segments take equal lengths, simply reflecting the gradient of the derivative(s). Defaults to "equal".
#' @param tail.length multiplies the current length of the tail (both proportional and equal arrow.types) by the specified factor.  The argument defaults to 1, which is length of the longest vector within the domain boundaries (i.e., the entire field).
#' @param ... passes arguments to both \code{\link{plot}} and \code{\link{arrows}}.
#' @keywords vector field plot, remainder field plot
#' 
#' @examples
#' #Example 1 from article
#' equationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
#' equationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
#' # 0.6.1 vector field
#' VDV <- VecDecompVec(x.num.steps=4100, y.num.steps=4100, x.rhs=equationx, y.rhs=equationy, x.bound=c(-0.5,20), y.bound=c(-0.5,20))
#' VecDecompPlot(field=list(VDV[,,1],VDV[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))
#'
########################################################################
# This is off until we decide what to do about e1.global
########################################################################
# # 0.6.2 gradient field	
# VDG <- VecDecompGrad(e1.global)
# VecDecompPlot(field=list(VDG[,,1],VDG[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))
#
# # 0.6.3 remainder field
# VDR <- VecDecompRem(surface=e1.global, x.rhs=equationx, y.rhs=equationy, x.bound=c(-0.5,20), y.bound=c(-0.5,20))
# VecDecompPlot(field=list(VDR[,,1],VDR[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))


VecDecompPlot <- function(field, dens, x.bound, y.bound, x.lim, y.lim, arrow.type="equal", tail.length=1, ...){
		x.range <- max(x.bound)-min(x.bound)
		y.range <- max(y.bound)-min(y.bound)

		row.range <- nrow(field[[1]])-1
		col.range <- ncol(field[[1]])-1

	if(missing(x.lim) | missing(y.lim)) {
		row.min <- min(which(field[[1]] != 0 , arr.ind = T)[,1])
		row.max <- max(which(field[[1]] != 0 , arr.ind = T)[,1])
		col.min <- min(which(field[[1]] != 0 , arr.ind = T)[,2])
		col.max <- max(which(field[[1]] != 0 , arr.ind = T)[,2])

		x.min <- ((row.min-1)/row.range)*x.range + min(x.bound)
		x.max <- ((row.max-1)/row.range)*x.range + min(x.bound)
		y.min <- ((col.min-1)/col.range)*y.range + min(y.bound)
		y.max <- ((col.max-1)/col.range)*y.range + min(y.bound)

		x.win <- c(x.min,x.max)
		y.win <- c(y.min,y.max)

		print(x.win)
		print(y.win)
	} else {
		x.win <- x.lim
		y.win <- y.lim
	
		row.min <- (min(x.win)-min(x.bound))/x.range*row.range + 1
		row.max <- (max(x.win)-min(x.bound))/x.range*row.range + 1
		col.min <- (min(y.win)-min(y.bound))/y.range*col.range + 1
		col.max <- (max(y.win)-min(y.bound))/y.range*col.range + 1
	}

 	sub.x <- seq(row.min, row.max, length.out=dens[1])
	sub.y <- seq(col.min, col.max, length.out=dens[2])

	sub.x.val <- ((sub.x-1)/row.range)*x.range + min(x.bound)
	sub.y.val <- ((sub.y-1)/col.range)*y.range + min(y.bound)

	dx.sub <- field[[1]][sub.x, sub.y]
	dy.sub <- field[[2]][sub.x, sub.y]

	if(arrow.type=="proportional"){
	dx.rel <- (dx.sub/max(((dx.sub^2)+(dy.sub^2))^0.5, na.rm = T))
	dy.rel <- (dy.sub/max(((dx.sub^2)+(dy.sub^2))^0.5, na.rm = T))
	dx.plot <- dx.rel*tail.length
	dy.plot <- dy.rel*tail.length
	}

	if(arrow.type=="equal"){
	dx.even <- dx.sub/((dx.sub^2)+(dy.sub^2))^0.5
	dy.even <- dy.sub/(((dx.sub^2)+(dy.sub^2))^0.5)
	dx.plot <- dx.even*tail.length
	dy.plot <- dy.even*tail.length
	}


	if(arrow.type != "proportional" & arrow.type != "equal"){
	dx.plot <- dx.sub*tail.length
	dy.plot <- dy.sub*tail.length
	}

 	qpr <- nrow(dx.plot)
	qpc <- ncol(dy.plot)

	plot(0 , type = "n" , xlim = x.win , ylim = y.win , las = 1, ...)
	for (j in 1:qpr){
			for (i in 1:qpc){
				x0 <- sub.x.val[j] - (dx.plot[j,i]/2)
				x1 <- sub.x.val[j] + (dx.plot[j,i]/2)
				y0 <- sub.y.val[i] - (dy.plot[j,i]/2)
				y1 <- sub.y.val[i] + (dy.plot[j,i]/2)
				arrows(x0,y0,x1,y1, ...)
			}
		}
}
