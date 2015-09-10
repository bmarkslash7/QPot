#' Plot simulation of two-dimensional stochastic differential equations
#'
#' This function allows you to plot the simulation of simulate two-dimensional stochastic differential equations from \code{\link{TSTraj}}
#' @param mat a matrix output from \code{\link{TSTraj}}.
#' @param deltat numeric value indicating the frequency of stochastic perturbation, as \eqn{\Delta t}, used in the function to recaluculate axes if applicable.
#' @param dim dimensions of the plot; \code{dim=1} to plot a timeseries with \code{X} and \code{Y} on the ordinate axis or \code{dim=2} to plot the trjectories in state space (i.e., \code{X} and \code{Y} respectively on the abscissa and ordinate axes).
#' @param y.lim for \code{dim=1}, allows user to specify the range of the y-axis as a two-element vector.
#' @param x.lab for \code{dim=1}, allows user to specify the axis as "time" or "steps," with steps being \eqn{time \times \Delta t}
#' @param dens if \code{dens=TRUE}, plots a horizontal one-dimensional density plot adjacent to the timerseries.
#' @param lwd line width.
#' @param line.alpha transparency of lines from 0--255.
#' @param zero.axes if TRUE, then axes plotted at \code{X=0} and \code{Y=0}.
#' @keywords plot stochastic simulations
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
#' ts.out.ex1 <- TSTraj(y0= model.state, time=model.time, deltat=model.deltat, 
#'   x.rhs=test.eqn.x, y.rhs= test.eqn.y, parms=model.parms, sigma=model.sigma)
#' # Plot as one-dimensional plot
#' TSPlot(ts.out.ex1, deltat=model.deltat)
#' # . . . or plot as two-dimensional plot
#' TSPlot(ts.out.ex1, deltat=model.deltat, dim=2)

	TSPlot <- function(mat, deltat , dim = 1 , y.lim = NA , x.lab = "time" , dens = TRUE , lwd = 2 , line.alpha = 130 , zero.axes = TRUE, ...) {
		if (dim ==1) {
			if (table(is.infinite(mat))["FALSE"] != nrow(mat)*3 ) { # if Inf values in the timeseries
					message("Simulation -> Inf. Try: (i) set exact y-axis limits using the y.lim argument and (ii) dens = FALSE")
				if (is.numeric(y.lim) == TRUE) {
					y.coords <- y.lim
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new=FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], mat[,2], type = "l" , ylim = y.coords , las = 1 , xlab = x.label , ylab = "State variables" , col = rgb(255,0,0, line.alpha,maxColorValue=255) , lwd = lwd , xaxt = "n")
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,maxColorValue=255) , lwd = lwd)
					if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat)) , ylim = y.coords , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n")
							axis(1)
						}
						if(dens == TRUE) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						y.lims <- y.coords
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , ylim = y.lims , yaxt = "n")
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,maxColorValue=255) , border = rgb(255,0,0,130,maxColorValue=255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,maxColorValue=255) , border = rgb(0,0,255,130,maxColorValue=255))
						axis(4 , las = 1)}
				} else {
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new=FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], mat[,3], type = "l"  , las = 1 , xlab = x.label , ylab = "State variables" , col = rgb(255,0,0,line.alpha,maxColorValue=255) , lwd = lwd , xaxt = "n")
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,maxColorValue=255) , lwd = lwd)
						if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat)) , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n")
							axis(1)
						}
						if(dens == T) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , yaxt = "n")
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,maxColorValue=255) , border = rgb(255,0,0,130,maxColorValue=255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,maxColorValue=255) , border = rgb(0,0,255,130,maxColorValue=255))
						axis(4 , las = 1)}
				}
			} else {
					ifelse(is.numeric(y.lim) == TRUE , y.coords <- y.lim , y.coords <- c(min(mat[,2:3],na.rm = T) , max(mat[,2:3],na.rm = T)) )
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new = FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], mat[,2], type = "l" , ylim = y.coords , las = 1 , ylab = "State variables" , xlab = x.label , col = rgb(255,0,0,line.alpha,maxColorValue=255) , lwd = lwd , xaxt = "n")
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,maxColorValue=255) , lwd = lwd)
					if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat)) , ylim = y.coords , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n")
							axis(1)
						}
						if(dens == TRUE) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						y.lims <- y.coords
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , ylim = y.lims , yaxt = "n")
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,maxColorValue=255) , border = rgb(255,0,0,130,maxColorValue=255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,maxColorValue=255) , border = rgb(0,0,255,130,maxColorValue=255))
						axis(4 , las = 1)}
			}
		} else { # dim != 1
			x.lim <- c(min(mat[,2]) , max(mat[,2]))
			y.lim <- c(min(mat[,3]) , max(mat[,3]))
			plot(mat[,2] , mat[,3] , type = "n", las = 1 , ylim = y.lim , xlim = x.lim, ...)
			if (zero.axes == TRUE) {abline(v = 0 , h = 0 , col = "grey75" , lwd = 0.75 , lty = 1)}
			lines(mat[,2] , mat[,3] , col = rgb(50,50,50,line.alpha,maxColorValue=255) , lwd = lwd)
		}
		}
