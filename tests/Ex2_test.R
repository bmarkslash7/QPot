##### 0.1 preperation #####
	# clean up and set seed
	rm(list=ls())
	your.favourite.number <- 6174
	set.seed(your.favourite.number)
	def.par <- par() # for deefault par settings for resetting plottting functions

	#load libraries
	library(QPot)


##### 0.2 stochastic simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		test.eqn.x = "-(y-beta) + mu*(x-alpha)*(1 - (x-alpha)^2 - (y-beta)^2) "
		test.eqn.y = "(x-alpha) + mu*(y-beta)*(1 - (x-alpha)^2 - (y-beta)^2)"

	# 0.2.1 parameters
		model.state <- c(x=3, y=3)
		model.parms <- c(alpha=4, beta=5, mu=0.2)
		model.sigma <- 0.1
		model.time <- 250
		model.deltat <- 0.005

	# 0.2.2 time series
		ts.out.ex2 <- TSTraj(y0= model.state, time=model.time, deltat=model.deltat, x.rhs=test.eqn.x, y.rhs= test.eqn.y, parms=model.parms, sigma=model.sigma)

	# 0.2.3 time series plots
		TSPlot(ts.out.ex2, deltat=model.deltat)
		par(def.par)
		TSPlot(ts.out.ex2, deltat=model.deltat, dim=2)
		TSDensity(ts.out.ex2, dim=1)
		TSDensity(ts.out.ex2, dim=2)


##### 0.3 local quasi-potential!!! #####
		testequationx = "-(y-5) + 0.2*(x-4)*(1 - (x-4)^2 - (y-5)^2)"
		testequationy = "(x-4) + 0.2*(y-5)*(1 - (x-4)^2 - (y-5)^2)"
		xbounds = c(-0.5, 7.5)
		ybounds = c(-0.5, 7.5)
		xstepnumber = 2000
		ystepnumber = 2000
		xinit = 4.15611
		yinit = 5.987774

		exampletwo_eq1 <- QPotential(x.rhs = testequationx, x.start = xinit, x.bound = xbounds, x.num.steps = xstepnumber, y.rhs = testequationy, y.start = yinit, y.bound = ybounds, y.num.steps = ystepnumber)


##### 0.4 global quasi-potential!!! #####
	e2.global <- exampletwo_eq1
	rm(exampletwo_eq1)
	gc()


##### 0.5 quasi-potential vizualization!!! #####
	QPContour(e2.global,c(1000,1000),c(-0.5,20),c(-0.5,20),c.parm=5)


##### 0.6 vector field decompisition!!! #####
	# 0.6.1 vector field
	VDV <- VecDecompVec(x.num.steps=xstepnumber, y.num.steps=ystepnumber, x.rhs=testequationx, y.rhs=testequationy, x.bound=xbounds, y.bound=ybounds)
		VecDecompPlot(field=list(VDV[,,1],VDV[,,2]), dens=c(25,25), x.bound=xbounds, y.bound=ybounds , tail.length = 0.5 , length = 0.03 , arrow.type="proportional", x.lim=c(2,6), y.lim=c(3,7))

	# 0.6.2 gradient field	
	VDG <- VecDecompGrad(e2.global)
		VecDecompPlot(field=list(VDG[,,1],VDG[,,2]), dens=c(25,25), x.bound=xbounds, y.bound= ybounds , tail.length = 0.1 , length = 0.03 , arrow.type = "proportional", x.lim=c(2,6), y.lim=c(3,7))

	# 0.6.3 remainder field
	VDR <- VecDecompRem(surface=e2.global, x.rhs=testequationx, y.rhs=testequationy, x.bound=xbounds, y.bound=ybounds)
		VecDecompPlot(field=list(VDR[,,1],VDR[,,2]), dens=c(25,25), x.bound=xbounds, y.bound=xbounds , tail.length = 0.2 , length = 0.03 , arrow.type = "proportional", x.lim=c(2,6), y.lim=c(3,7))
