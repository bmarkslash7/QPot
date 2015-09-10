# Below is a pipe every 10 characters.  The R Journal suggests that code is wrapped around 80 characters
#*********|*********|*********|*********|*********|*********|*********|*********

##### 0.1 preperation #####
	# clean up and set seed
	rm(list=ls())
	your.favourite.number <- 3818919 #chris
	set.seed(your.favourite.number)

	#load libraries
	library(QPot)

##### 0.2 stochasticbounds.x simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		var.eqn.x <- "x*((1+alpha1)-x*x-x*y-y*y)"
		var.eqn.y <- "y*((1+alpha2)-x*x-x*y-y*y)"

	# 0.2.1 parameters
		model.state <- c(x = 0.5, y = 0.5)
		model.parms <- c(alpha1 = 1.25, alpha2 = 2)
		model.sigma <- 0.8
		model.time <- 5000
		model.deltat <- 0.01

	# 0.2.2 time series
	ts.ex3 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)

	# 0.2.3 time series plots
		TSPlot(ts.ex3, deltat = model.deltat)
		TSPlot(ts.ex3, deltat = model.deltat, dim = 2 , line.alpha = 5)
		TSDensity(ts.ex3, dim = 1)
		TSDensity(ts.ex3, dim = 2 , contour.levels = 20 , contour.lwd = 0.1)

		# # plots for the paper figures
			# print.wd <- "/Users/christophermoore/DropBox/QPRPackage/QPotPaper/Figures/"
			# # Ex3_TS_1D.png
			# png(paste(print.wd,"Ex3_TS_1D.png",sep=""), width = 720, height = 480)
			# TSPlot(ts.ex3, deltat = model.deltat)
			# dev.off()
			# # Ex3_TS_2D.png
			# png(paste(print.wd,"Ex3_TS_2D.png",sep=""), width = 500, height = 500)
			# TSPlot(ts.ex3, deltat = model.deltat, dim = 2, line.alpha = 25, lwd = 1, xlab = expression(italic(x)), ylab = expression(italic(y)), zero.axes = F)
			# dev.off()
			# # Ex3_Dens_2D.png
			# png(paste(print.wd,"Ex3_Dens_2D.png",sep=""), width = 500, height = 500)
			# TSDensity(ts.ex3, deltat = model.deltat, dim = 2, line.alpha = 50, lwd = 1, xlab = expression(italic(x)), ylab = expression(italic(y)), contour.lwd = 0.75)
			# dev.off()


##### 0.3 local quasi-potential!!! #####
		equation.x = "x*((1+1.25)-x*x-x*y-y*y)"
		equation.y = "y*((1+2)-x*x-x*y-y*y)"
		bounds.x = c(-3, 3)
		bounds.y = c(-3, 3)
		step.number.x = 6000
		step.number.y = 6000
		eq1.x = 0
		eq1.y = -1.73205
		eq2.x = 0
		eq2.y = 1.73205

eq1.local <- QPotential(x.rhs = equation.x, x.start = eq1.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq1.y, y.bound = bounds.y, y.num.steps = step.number.y)
eq2.local <- QPotential(x.rhs = equation.x, x.start = eq2.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq2.y, y.bound = bounds.y, y.num.steps = step.number.y)


##### 0.4 global quasi-potential!!! #####
ex3.global <- QPGlobal(local.surfaces = list(eq1.local, eq2.local),unstable.eq.x = c(0, -1.5, 1.5), unstable.eq.y = c(0, 0, 0), x.bound = bounds.x, y.bound = bounds.y)


##### 0.5 quasi-potential vizualization!!! #####
QPContour(ex3.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 1)

	# figure for the paper
		# print.wd <- "/Users/christophermoore/DropBox/QPRPackage/QPotPaper/Figures/"
		# # Ex3_QP_contour.png
		# png(paste(print.wd,"Ex3_QP_contour.png",sep=""), width = 500, height = 500)
		# QPContour(ex3.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 5, xlab=expression(italic(x)), ylab=expression(italic(y)))
		# dev.off()


##### 0.6 vector field decompisition!!! #####
	# 0.6.0 all fields
	VDAll <- VecDecomAll(surface = ex3.global, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
	VecDecomPlot(field = list(VDAll[,,1], VDAll[,,2]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "equal", tail.length = 0.2, head.length = 0.03)
	VecDecomPlot(field = list(VDAll[,,3], VDAll[,,4]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 1/3, head.length = 0.03)
	VecDecomPlot(field = list(VDAll[,,5], VDAll[,,6]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.5, head.length = 0.03)

	# 0.6.1 vector field
	VDV <- VecDecomVec(x.num.steps = step.number.y, y.num.steps= step.number.y, x.rhs=equation.x, y.rhs=equation.y, x.bound= bounds.x, y.bound= bounds.x)
		VecDecomPlot(field=list(VDV[,,1],VDV[,,2]), dens=c(50,50), x.bound= bounds.x, y.bound=bounds.y , tail.length = 0.2 , length = 0.03)

	# 0.6.2 gradient field	
	VDG <- VecDecomGrad(e3.global)
		VecDecomPlot(field=list(VDG[,,1],VDG[,,2]), dens=c(25,25), x.bound= bounds.x, y.bound= bounds.y , tail.length = 0.5 , length = 0.03 , arrow.type = "proportional")

	# 0.6.3 remainder field
	VDR <- VecDecomRem(surface=e3.global, x.rhs=equation.x, y.rhs=equation.y, x.bound=bounds.x, y.bound=bounds.y)
		VecDecomPlot(field=list(VDR[,,1],VDR[,,2]), dens=c(25,25), x.bound= bounds.y, y.bound=bounds.y , tail.length = 0.5 , length = 0.03 , arrow.type = "proportional")

##### 0.7 3D graphs!!! #####
	library(rgl)
	dens.sub <- c(1000,1000)
	global.sub <- ex3.global[round(seq(1,nrow(ex3.global),length.out=dens.sub[1])) , round(seq(1,ncol(ex3.global),length.out=dens.sub[2]))]
	# global.sub[global.sub > 2] <- NA # to cut off large values for a better view of the basin
	persp3d(x = global.sub, col = "orange", expand = 1.1)