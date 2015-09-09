# Below is a pipe every 10 characters.  The R Journal suggests that code is wrapped around 80 characters
#*********|*********|*********|*********|*********|*********|*********|*********

##### 0.1 preperation #####
	# clean up and set seed
	rm(list=ls())
	your.favourite.number <- 3818919 #chris
	set.seed(your.favourite.number)

	#load libraries
	library(QPot)


##### 0.2 stochastic simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		var.eqn.x <- "-(y-beta) + mu*(x-alpha)*(1-(x-alpha)^2-(y-beta)^2) "
		var.eqn.y <- "(x-alpha) + mu*(y-beta)*(1-(x-alpha)^2-(y-beta)^2)"

	# 0.2.1 parameters
		model.state <- c(x = 3, y = 3)
		model.parms <- c(alpha = 4, beta = 5, mu = 0.2)
		model.sigma <- 0.1
		model.time <- 2500
		model.deltat <- 0.005

	# 0.2.2 time series
ts.ex2 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)

	# 0.2.3 time series plots
		TSPlot(ts.ex2, deltat = model.deltat)
		quartz() # new quartz device window
		TSPlot(ts.ex2, deltat = model.deltat, dim = 2, line.alpha = 25)
		TSDensity(ts.ex2, dim = 1)
		TSDensity(ts.ex2, dim = 2)


##### 0.3 local quasi-potential!!! #####
		equation.x = "-(y-5) + 0.2*(x-4)*(1-(x-4)^2-(y-5)^2)"
		equation.y = "(x-4) + 0.2*(y-5)*(1-(x-4)^2-(y-5)^2)"
		bounds.x = c(-0.5, 7.5)
		bounds.y = c(-0.5, 7.5)
		step.number.x = 4000
		step.number.y = 4000
		xinit = 4.15611
		yinit = 5.987774

eq1.qp <- QPotential(x.rhs = equation.x, x.start = xinit, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = yinit, y.bound = bounds.y, y.num.steps = step.number.y)


##### 0.4 global quasi-potential!!! #####
	#same as the local quasi-potential, calculated above


##### 0.5 quasi-potential vizualization!!! #####
QPContour(eq1.qp, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 10)


##### 0.6 vector field decompisition!!! #####
	# 0.6.1 vector field
	VDV <- VecDecompVec(x.num.steps = step.number.x, y.num.steps = step.number.y, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
		VecDecompPlot(field = list(VDV[,,1],VDV[,,2]), dens = c(25,25), x.bound = bounds.x, y.bound = bounds.y, tail.length = 0.5, head.length = 0.03, arrow.type = "proportional", x.lim = c(2, 6), y.lim = c(3, 7), xlab = expression(italic(X)), ylab = expression(italic(Y)))

	# 0.6.2 gradient field	
	VDG <- VecDecompGrad(eq1.qp)
		VecDecompPlot(field = list(VDG[,,1], VDG[,,2]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, tail.length = 0.2, head.length = 0.05 , arrow.type = "proportional", x.lim = c(2, 6), y.lim = c(3, 7), xlab = expression(italic(X)), ylab = expression(italic(Y)))

	# 0.6.3 remainder field
	VDR <- VecDecompRem(surface = eq1.qp, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.x)
		VecDecompPlot(field = list(VDR[,,1], VDR[,,2]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.x , tail.length = 0.2 , head.length = 0.03 , arrow.type = "proportional", x.lim = c(2, 6), y.lim = c(3, 7), xlab = expression(italic(X)), ylab = expression(italic(Y)))
