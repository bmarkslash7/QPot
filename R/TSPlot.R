TSPlot <- function(mat, deltat, dim = 1, xlim = 'NULL', ylim = 'NULL', x.lab = "time", dens = TRUE, lwd = 2, line.alpha = 130, zero.axes = TRUE, ...) {
		if (missing(deltat) == TRUE) {stop("deltat is missing and needed to compute steps and 1-D hist.  Please specify.")}
		global.min <- min(mat[,2:3], na.rm = T)
		global.max <- max(mat[,2:3], na.rm = T)
		x.min <- min(mat[,2], na.rm = T)
		x.max <- max(mat[,2], na.rm = T)
		y.min <- min(mat[,3], na.rm = T)
		y.max <- max(mat[,3], na.rm = T)
		if (dim == 1) {
			if (table(is.infinite(mat))["FALSE"] != nrow(mat)*3 ) { # if Inf values in the timeseries
					message("Simulation -> Inf. Try: (i) set exact y-axis limits using the y.lim argument and (ii) dens = FALSE")
				if (missing(ylim)) {
					ylim = c(global.min, global.max)
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new=FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], type = "l" , ylim = ylim , las = 1 , xlab = x.label , ylab = "State variables" , col = rgb(255, 0, 0, line.alpha, NULL, 255) , lwd = lwd , xaxt = "n", ...)
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,NULL,255) , lwd = lwd)
					if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat)) , ylim = ylim , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n", ...)
							axis(1)
						}
						if(dens == TRUE) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						y.lims <- ylim
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , ylim = y.lims , yaxt = "n")
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,NULL,255) , border = rgb(255,0,0,130,NULL,255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,NULL,255) , border = rgb(0,0,255,130,NULL,255))
						axis(4 , las = 1)}
				} else {
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new=FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], type = "l", las = 1 , ylim = ylim, xlab = x.label , ylab = "State variables" , col = rgb(255,0,0,line.alpha,NULL,255) , lwd = lwd , xaxt = "n", ...)
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,NULL,255) , lwd = lwd)
						if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat)), ylim = ylim , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n", ...)
							axis(1)
						}
						if(dens == T) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , yaxt = "n", ylim = ylim)
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,NULL,255) , border = rgb(255,0,0,130,NULL,255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,NULL,255) , border = rgb(0,0,255,130,NULL,255))
						axis(4 , las = 1)}
				}
			} else { # No Inf values in timeseries
					if (missing(ylim)) {ylim <- c(global.min, global.max)}
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new = FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], type = "n", ylim = ylim , las = 1 , ylab = "State variables" , xlab = x.label , lwd = lwd , xaxt = "n", ...)
					lines(mat[,1] , mat[,2] , col = rgb(255,0,0,line.alpha,NULL,255) , lwd = lwd)
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,NULL,255) , lwd = lwd)
					if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat))  , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n", ...)
							axis(1)
						}
						if(dens == TRUE) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , ylim = ylim , yaxt = "n")
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,NULL,255) , border = rgb(255,0,0,130,NULL,255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,NULL,255) , border = rgb(0,0,255,130,NULL,255))
						axis(4 , las = 1)}
			}
		} else { # dim != 1
			if(missing(xlim)) {xlim = c(min(mat[,2]) , max(mat[,2]))}
			if(missing(ylim)) {ylim = c(min(mat[,3]) , max(mat[,3]))}
			plot(mat[,2] , mat[,3] , type = "n", las = 1 , ylim = ylim , xlim = xlim, ...)
			if (zero.axes == TRUE) {abline(v = 0 , h = 0 , col = "grey75" , lwd = 0.75 , lty = 1)}
			lines(mat[,2] , mat[,3] , col = rgb(50,50,50,line.alpha,NULL,255) , lwd = lwd)
		}
		}
