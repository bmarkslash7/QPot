#' Contour plot of quasi-potential surfaces
#'
#' This function allows you to plot quasi-potential surfaces
#' @param surface the surface to be plotted, from \code{\link{QPGlobal}}
#' @param density vector respectively for the number of \code{x} and \code{y} points to be plotted.
#' @param x.lim two-element vector for the x variable's minimum and maximum values
#' @param y.lim two-element vector for the y variable's minimum and maximum values
#' @param n.filled.contour numeric value for the nubmber of breaks in the filled contour.
#' @param n.contour.lines numeric value for the nubmber of breaks in the contour lines.
#' @param col vector of colors to be used in the plot.
#' @param c contour line adjustment (see details).
#' @param contour.lines if TRUE, then contour lines plotted over filled contour; vice versa if FALSE.
#' @details Because, in gerneal, capturing the topological features can be subtle, we implimented a feature in \code{\link{QPContour}} to keep the filled.contour region while changing the contour lines.  Specifically, filled.contour takes the range of the surface values (\eqn{\phi}), divides by the number of the specified contours (i.e., \code{n.filled.contour}), and creates a contour at each break, which happenes to be equal across the range.  But because visualizing some topology may (i) require looking between contour breaks and (ii) adding contour lines would overload the plot with lines, we use an equation to modify the distribution of contour lines.  Namely, adjusting the \code{c} argument in the \code{\link{QPContour}} function adjusts the \eqn{c} paramter in the following equation: \deqn{max_\phi \times (\frac{x}{n-1})^c.}.  This allows the user to keep the same number of contour lines (i.e., specified with \code{n.contour.lines}), but focus them toward the troughs or peaks of the surfaces. At \eqn{c=1}, the contour lines correspond to the filled.contour breaks.  If \eqn{c > 1}, then the contour lines become more concentrated towards the trough.  Similarly, if \eqn{c < 1}, then the contour lines are more focused at the peaks of the surface.  As an example: \cr \figure{Example3.png}.

QPContour <- function(surface, density, x.lim, y.lim, n.filled.contour, n.contour.lines, col, contour.lines = TRUE){
	x.range <- max(x.lim)-min(x.lim)
	y.range <- max(y.lim)-min(y.lim)
	eq <- global.qp
	plot.mesh.xy <- c(25,25)
	sub.x <- round(seq(1,ncol(eq),length.out=plot.mesh.xy[1]))
	sub.y <- round(seq(1,nrow(eq),length.out=plot.mesh.xy[2]))
	eq.sub <- eq[sub.x,sub.y]

################# 4.0 grrrrraph ##################
	rowmin <- min(which(eq.sub > 0 , arr.ind = T)[,1])
	rowmax <- max(which(eq.sub > 0 , arr.ind = T)[,1])
	colmin <- min(which(eq.sub > 0 , arr.ind = T)[,2])
	colmax <- max(which(eq.sub > 0 , arr.ind = T)[,2])

	# x.pretty <- pretty(c(rowmin/plot.mesh.xy[1]*x.range-0.5,rowmax/plot.mesh.xy[1]*x.range-0.5))
	# x.pretty.at <- (x.pretty-xlim[1])/x.range*plot.mesh.xy[1]
	# y.pretty <- pretty(c(colmin/plot.mesh.xy[2]*y.range-0.5,colmax/plot.mesh.xy[2]*y.range-0.5))
	# y.pretty.at <- (y.pretty-ylim[1])/y.range*plot.mesh.xy[2]

	x.pretty <- pretty(c(rowmin/plot.mesh.xy[1]*x.range+xlim[1],rowmax/plot.mesh.xy[1]*x.range+xlim[2]))
	x.pretty.at <- (x.pretty-xlim[1])/x.range*plot.mesh.xy[1]
	y.pretty <- pretty(c(colmin/plot.mesh.xy[2]*y.range+ylim[1],colmax/plot.mesh.xy[2]*y.range+ylim[2]))
	y.pretty.at <- (y.pretty-ylim[1])/y.range*plot.mesh.xy[2]


	plot(0 , type = "n" , xlim = c(rowmin,rowmax) , ylim = c(colmin,colmax) , las = 1 , xlab = expression(italic(X)) , ylab = expression(italic(Y)) , xaxt = "n" , yaxt = "n")

	axis(1 , at = x.pretty.at , labels = x.pretty)
	axis(2 , at = y.pretty.at , labels = y.pretty , las = 1)

	min.eq.sub <- min(eq.sub , na.rm = T)
	max.eq.sub <- max(eq.sub , na.rm = T)
	contour.levels <- 25
	contour.breaks <- seq(min.eq.sub , max.eq.sub , length = contour.levels)
	eq.max <- max(eq,na.rm = T)
	line.contour.breaks <- (eq.max)*(((0:contour.levels)/(contour.levels-1)))^2.5 # 3	
	myRmap <- colorRampPalette(c("#FDE725FF","#E3E418FF","#C7E020FF","#ABDC32FF","#8FD744FF","#75D054FF","#5DC963FF","#47C06FFF","#35B779FF","#28AE80FF","#20A486FF","#1F9A8AFF","#21908CFF","#24868EFF","#287C8EFF","#2C728EFF" ,"#31688EFF","#355D8DFF","#3B528BFF","#404688FF","#443A83FF","#472D7BFF","#481F71FF","#471163FF","#440154FF"))(contour.levels)

	.filled.contour(colmin:colmax, rowmin:rowmax , eq.sub[rowmin:rowmax , colmin:colmax] , levels = contour.breaks , col = myRmap)
	contour(1:ncol(eq.sub) , 1:nrow(eq.sub) , eq.sub , levels = line.contour.breaks, drawlabels = F ,  add = TRUE , col = "black" , lwd = 0.5)
	}
