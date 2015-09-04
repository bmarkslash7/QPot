#' Contour plot of quasi-potential surfaces
#'
#' This function allows you to plot quasi-potential surfaces
#' @param surface the surface to be plotted, from \code{\link{QPGlobal}}
#' @param dens vector respectively for the number of \code{x} and \code{y} points to be plotted.
#' @param x.bound two-element vector for the surface's minimum and maximum x values
#' @param y.bound two-element vector for the surface's minimum and maximum y values
#' @param n.filled.contour numeric value for the nubmber of breaks in the filled contour.
#' @param n.contour.lines numeric value for the number of breaks in the contour lines.
#' @param c.parm contour line adjustment (see details).
#' @param col.contour a vector of colors used for the filled contour region.
#' @param contour.lines if TRUE, then contour lines plotted over filled contour; vice versa if FALSE.
#' @param ... passes arguments to both \code{\link{plot}} and \code{\link{contour}}.
#' @details Because, in general, capturing the topological features can be subtle, we implemented a feature in \code{\link{QPContour}} to keep the filled.contour region while changing the contour lines.  Specifically, filled.contour takes the range of the surface values (\eqn{\phi}), divides by the number of the specified contours (i.e., \code{n.filled.contour}), and creates a contour at each break, which happenes to be equal across the range.  But because visualizing some topology may (i) require looking between contour breaks and (ii) adding contour lines would overload the plot with lines, we use an equation to modify the distribution of contour lines.  Namely, adjusting the \code{c} argument in the \code{\link{QPContour}} function adjusts the \eqn{c} paramter in the following equation: \deqn{max_\phi \times (\frac{x}{n-1})^c.}.  This allows the user to keep the same number of contour lines (i.e., specified with \code{n.contour.lines}), but focus them toward the troughs or peaks of the surfaces. At \eqn{c=1}, the contour lines correspond to the filled.contour breaks.  If \eqn{c > 1}, then the contour lines become more concentrated towards the trough.  Similarly, if \eqn{c < 1}, then the contour lines are more focused at the peaks of the surface.  As an example: \cr \figure{Example3.png}.
#' 
########################################################################
# COMMENTED OUT FOR PASSING devtools::check() - DECLARE local.1 and local.2
########################################################################
# @examples
# # First, use a surface (example from QPGlobal)
# global.qp <- QPGlobal(list(local.1,local.2),c(0,4),c(0,4),c(-1,5),c(-1,5))
#
# # Second, input that surface into QPContour
# QPContour(surface=global.qp, dens=c(100,100), y.bound=c(-0.5,20), 
# y.bound=c(-0.5,20), n.filled.contour=20, n.contour.lines=20,
# col.contour=c("red", "white", "blue"), contour.lines = TRUE)

QPContour <- function(surface, dens, x.bound, y.bound, n.filled.contour=25, n.contour.lines=25, c.parm=1, col.contour, contour.lines = TRUE, ...){
	x.range <- max(x.bound)-min(x.bound)
	y.range <- max(y.bound)-min(y.bound)

	row.range <- nrow(surface)-1
	col.range <- ncol(surface)-1

	row.min <- min(which(surface != 0 , arr.ind = T)[,1])
	row.max <- max(which(surface != 0 , arr.ind = T)[,1])
	col.min <- min(which(surface != 0 , arr.ind = T)[,2])
	col.max <- max(which(surface != 0 , arr.ind = T)[,2])

	x.min <- ((row.min-1)/row.range)*x.range + min(x.bound)
	x.max <- ((row.max-1)/row.range)*x.range + min(x.bound)
	y.min <- ((col.min-1)/col.range)*y.range + min(y.bound)
	y.max <- ((col.max-1)/col.range)*y.range + min(y.bound)

	x.win <- c(x.min,x.max)
	y.win <- c(y.min,y.max)

 	sub.x <- seq(row.min, row.max, length.out=dens[1])
	sub.y <- seq(col.min, col.max, length.out=dens[2])

	sub.x.val <- ((sub.x-1)/row.range)*x.range + min(x.bound)
	sub.y.val <- ((sub.y-1)/col.range)*y.range + min(y.bound)

	eq.sub <- surface[sub.x, sub.y]

	plot(0 , type = "n" , xlim = x.win , ylim = y.win , las = 1, ...)
	min.eq.sub <- min(eq.sub , na.rm = T)
	max.eq.sub <- max(eq.sub , na.rm = T)
	contour.breaks <- seq(min.eq.sub , max.eq.sub , length = n.filled.contour)
	eq.max <- max(surface,na.rm = T)
	line.contour.breaks <- (eq.max)*(((0:n.contour.lines)/(n.contour.lines-1)))^c.parm
	myRamp <- if(missing(col.contour)){colorRampPalette(c("#FDE725FF","#E3E418FF","#C7E020FF","#ABDC32FF","#8FD744FF","#75D054FF","#5DC963FF","#47C06FFF","#35B779FF","#28AE80FF","#20A486FF","#1F9A8AFF","#21908CFF","#24868EFF","#287C8EFF","#2C728EFF" ,"#31688EFF","#355D8DFF","#3B528BFF","#404688FF","#443A83FF","#472D7BFF","#481F71FF","#471163FF","#440154FF"))(n.filled.contour)}else{colorRampPalette(col.contour)(n.filled.contour)}
	.filled.contour(sub.x.val , sub.y.val , eq.sub , levels = contour.breaks , col = myRamp)
	if(contour.lines==TRUE){contour(sub.x.val , sub.y.val , eq.sub , levels = line.contour.breaks, drawlabels = F ,  add = TRUE , col = "black" , ...)}
	}