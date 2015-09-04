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
		test.eqn.x = "x*((1+alpha1) - x*x - x*y - y*y)"
		test.eqn.y = "y*((1+alpha2) - x*x - x*y - y*y)"

	# 0.2.1 parameters
		model.state <- c(x=0.5, y=0.5)
		model.parms <- c(alpha1=1.25, alpha2=2)
		model.sigma <- 0.8
		model.time <- 2500
		model.deltat <- 0.01

	# 0.2.2 time series
	ts.out.ex3 <- TSTraj(y0= model.state, time=model.time, deltat=model.deltat, x.rhs=test.eqn.x, y.rhs= test.eqn.y, parms=model.parms, sigma=model.sigma)

	# 0.2.3 time series plots
		TSPlot(ts.out.ex3, deltat=model.deltat)
		par(def.par)
		TSPlot(ts.out.ex3, deltat=model.deltat, dim=2 , line.alpha = 5)
		TSDensity(ts.out.ex3, dim=1)
		TSDensity(ts.out.ex3, dim=2 , contour.levels = 20 , contour.lwd=0.1)

##### 0.3 local quasi-potential!!! #####
		testequationx = "x*((1+1.25) - x*x - x*y - y*y)"
		testequationy = "y*((1+2) - x*x - x*y - y*y)"
		xbounds = c(-3, 3)
		ybounds = c(-3, 3)
		xstepnumber = 3000
		ystepnumber = 3000
		xinit1 = 0
		yinit1 = -1.73205
		xinit2 = 0
		yinit2 = 1.73205

		examplethree_eq1 <- QPotential(x.rhs = testequationx, x.start = xinit1, x.bound = xbounds, x.num.steps = xstepnumber, y.rhs = testequationy, y.start = yinit1, y.bound = ybounds, y.num.steps = ystepnumber)
		examplethree_eq2 <- QPotential(x.rhs = testequationx, x.start = xinit2, x.bound = xbounds, x.num.steps = xstepnumber, y.rhs = testequationy, y.start = yinit2, y.bound = ybounds, y.num.steps = ystepnumber)


##### 0.4 global quasi-potential!!! #####

	e3.global <- QPGlobal(local.surfaces = list(examplethree_eq1, examplethree_eq2),unstable.eq.x=c(0,-1.5,1.5), unstable.eq.y = c(0,0,0), x.bound=xbounds, y.bound=ybounds)

##### 0.5 quasi-potential vizualization!!! #####

	QPContour(e3.global,c(1000,1000), x.bound=xbounds,y.bound=ybounds,c.parm=1,lwd=0.5,xlab="X",ylab="Y")

##### 0.6 vector field decompisition!!! #####
	# 0.6.1 vector field
	VDV <- VecDecompVec(x.num.steps= xstepnumber, y.num.steps= ystepnumber, x.rhs=testequationx, y.rhs=testequationy, x.bound= xbounds, y.bound= xbounds)
		VecDecompPlot(field=list(VDV[,,1],VDV[,,2]), dens=c(50,50), x.bound= xbounds, y.bound=ybounds , tail.length = 0.1 , length = 0.03)

	# 0.6.2 gradient field	
	VDG <- VecDecompGrad(e3.global)
		VecDecompPlot(field=list(VDG[,,1],VDG[,,2]), dens=c(25,25), x.bound= xbounds, y.bound= ybounds , tail.length = 0.5 , length = 0.03 , arrow.type = "proportional")

	# 0.6.3 remainder field
	VDR <- VecDecompRem(surface=e3.global, x.rhs=testequationx, y.rhs=testequationy, x.bound=xbounds, y.bound=ybounds)
		VecDecompPlot(field=list(VDR[,,1],VDR[,,2]), dens=c(25,25), x.bound= xbounds, y.bound=ybounds , tail.length = 0.5 , length = 0.03 , arrow.type = "proportional")