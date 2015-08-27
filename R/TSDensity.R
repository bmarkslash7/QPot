#' Density plot from simulation of two-dimensional stochastic differential equations
#'
#' This function allows you to create density plots for the simulation of two-dimensional stochastic differential equations in \code{\link{TSTraj}}
#' @param mat a matrix output from \code{\link{TSTraj}}.
#' @param dim dimensions of the plot; \code{dim=1} plots simple density histogram or \code{dim=2} plots the density in state space (i.e., \code{X} and \code{Y} respectively on the abscissa and ordinate axes).
#' @param contour.levels the number of contour levels for the two-dimensional plots (i.e., when \code{dim=2}).
#' @param col2d vector of colors to be used in the plot.
#' @param contour.lwd line width of contour lines if \code{contour.lines=TRUE}.
#' @param contour.lines if TRUE, then black countour lines added to the graph.
#' @keywords Density plot of stochastic simulations
#' 
#' @examples
#' model.state <- c(x=1 , y=2)
#' model.parms <- c(alpha=1.54, beta=10.14, delta=1, kappa=1, gamma=0.476, mu=0.112509)
#' model.sigma <- 0.05
#' model.time <- 100
#' model.deltat <- 0.2
#'
#' test.eqn.x = "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))"
#' test.eqn.y = "((gamma*(x^2)*y)/(kappa + (x^2))) - mu*(y^2)"
#'
#' ts.out.ex1 <- TSTraj(y0= model.state, time=model.time, deltat=model.deltat, x.rhs=test.eqn.x, y.rhs= test.eqn.y, parms=model.parms, sigma=model.sigma)
#'
#' # Plot as one-dimensional plot
#' TSDensity(ts.out.ex1, dim=1)
#' # . . . or plot as two-dimensional plot
#' TSDensity(ts.out.ex1, dim=2)

	TSDensity <- function(mat , dim = 1 , contour.levels = 15 ,  col2d = c("blue" , "yellow" , "orange" , "red") , contour.lwd = 0.5 ,  contour.lines = TRUE){
		if (dim ==1){
			densA <- density(mat[,2] , na.rm = T)
			densB <- density(mat[,3] , na.rm = T)
			y.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
			x.lims <- c(min(c(densA$x, densB$x)) , max(c(densA$x, densB$x)))
			plot(0 , type = "n" , ylab = "Density" , xlab = "State variables" , las = 1 , xlim = x.lims , ylim = y.lims)
			polygon(densA$x , densA$y, col = rgb(255,0,0,75,maxColorValue=255) , border = rgb(255,0,0,130,maxColorValue=255))
			polygon(densB$x , densB$y, col = rgb(0,0,255,75,maxColorValue=255) , border = rgb(0,0,255,130,maxColorValue=255))
			}
		if (dim ==2) {
#			require("MASS")	#Mass called in DESCRIPTION, Depends
			kern.2d <- MASS::kde2d(mat[,2] , mat[,3])
			x.max <- length(kern.2d$x)
			y.max <- length(kern.2d$y)
			x.range <- 1:x.max
			y.range <- 1:y.max
			contour.breaks <- seq(min(kern.2d$z) , max(kern.2d$z), length = contour.levels)
			myRmap <- colorRampPalette(col2d)(contour.levels)
			plot(0 , type = "n" , xlim = c(1 , x.max)  , ylim = c(1 , y.max) , las = 1  , xlab = colnames(mat)[2] , ylab = colnames(mat)[3] , xaxt = "n" , yaxt = "n")
			.filled.contour(x.range , y.range , kern.2d$z , levels = contour.breaks , col = myRmap)
			contour(x.range , y.range , kern.2d$z , levels = contour.breaks , col = myRmap , add = T ,  drawlabels = F)
			if (contour.lines == T) {contour(x.range , y.range , kern.2d$z , levels = contour.breaks, drawlabels = F ,  add = TRUE , col = "black" , lwd = contour.lwd)}
			par(new = TRUE)
			plot(0, type = "n" , xlim = c(min(kern.2d$x), max(kern.2d$x)) , ylim = c(min(kern.2d$y), max(kern.2d$y)) , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n")
			axis(1 , las = 1)
			axis(2 , las = 1)
		}
		if (dim != 1 & dim !=2 & dim != "both") { warning("Incorrect number of dimensions") }
		if (dim == "both") {
			par(mfrow=c(1,2))
			densA <- density(mat[,2] , na.rm = T)
			densB <- density(mat[,3] , na.rm = T)
			y.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
			x.lims <- c(min(c(densA$x, densB$x)) , max(c(densA$x, densB$x)))
			plot(0 , type = "n" , ylab = "Density" , xlab = "State variables" , las = 1 , xlim = x.lims , ylim = y.lims)
			polygon(densA$x , densA$y, col = rgb(255,0,0,75,maxColorValue=255) , border = rgb(255,0,0,130,maxColorValue=255))
			polygon(densB$x , densB$y, col = rgb(0,0,255,75,maxColorValue=255) , border = rgb(0,0,255,130,maxColorValue=255))

#			require("MASS")
			kern.2d <- MASS::kde2d(mat[,2] , mat[,3])
			x.max <- length(kern.2d$x)
			y.max <- length(kern.2d$y)
			x.range <- 1:x.max
			y.range <- 1:y.max
			contour.breaks <- seq(min(kern.2d$z) , max(kern.2d$z), length = contour.levels)
			myRmap <- colorRampPalette(col2d)(contour.levels)
			plot(0 , type = "n" , xlim = c(1 , x.max)  , ylim = c(1 , y.max) , las = 1  , xlab = colnames(mat)[2] , ylab = colnames(mat)[3] , xaxt = "n" , yaxt = "n")
			.filled.contour(x.range , y.range , kern.2d$z , levels = contour.breaks , col = myRmap)
			contour(x.range , y.range , kern.2d$z , levels = contour.breaks , col = myRmap , add = T , drawlabels=F)
			if (contour.lines == T) {contour(x.range , y.range , kern.2d$z , levels = contour.breaks, drawlabels = F ,  add = TRUE , col = "black" , lwd = contour.lwd)}
			par(new = TRUE)
			plot(0, type = "n" , xlim = c(min(kern.2d$x), max(kern.2d$x)) , ylim = c(min(kern.2d$y), max(kern.2d$y)) , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n")
			axis(1 , las = 1)
			axis(2 , las = 1)
		}
	}
