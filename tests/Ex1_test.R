##### 0.1 preperation #####
	# 0.1.0 clean up and set seed
	rm(list=ls())
	your.favourite.number <- 6174
	set.seed(your.favourite.number)


##### 0.2 stochastic simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		test.eqn.x = "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))"
		test.eqn.y = "((gamma*(x^2)*y)/(kappa + (x^2))) - mu*(y^2)"

	# 0.2.1 parameters
		model.state <- c(x=1 , y=2)
		model.parms <- c(alpha=1.54, beta=10.14, delta=1, kappa=1, gamma=0.476, mu=0.112509)
		model.sigma <- 0.05
		model.time <- 1000
		model.deltat <- 0.2

	# 0.2.2 time series
	ts.out.ex1 <- TSTraj(y0= model.state, time=model.time, deltat=model.deltat, x.rhs=test.eqn.x, y.rhs= test.eqn.y, parms=model.parms, sigma=model.sigma)

	# 0.2.3 time series plots
		TSPlot(ts.out.ex1, deltat=model.deltat)
		quartz() # new window with a single plotting area
		TSPlot(ts.out.ex1, deltat=model.deltat, dim=2)
		TSDensity(ts.out.ex1, dim=1)
		TSDensity(ts.out.ex1, dim=2)


##### 0.3 local quasi-potential!!! #####
	equation.x = "1.54*x*(1.0-(x/10.14))-(y*x*x)/(1.0+x*x)"
	equation.y = "((0.476*x*x*y)/(1+x*x))-0.112590*y*y"
	bounds.x = c(-0.5, 20.0)
	bounds.y = c(-0.5, 20.0)
	step.number.x = 4100
	step.number.y = 4100
	eq.1.1.x = 1.40491
	eq.1.1.y = 2.80808
	eq.1.2.x = 4.9040
	eq.1.2.y = 4.06187

	eq.1.1.local <- QPotential(x.rhs = equation.x, x.start = eq.1.1.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq.1.1.y,  y.bound = bounds.y, y.num.steps = step.number.y)
	eq.1.2.local <- QPotential(x.rhs = equation.x, x.start = eq.1.2.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq.1.2.y, y.bound = bounds.y, y.num.steps = step.number.y)


##### 0.4 global quasi-potential!!! #####

	eq.1.global <- QPGlobal(local.surfaces = list(eq.1.1.local,eq.1.2.local),unstable.eq.x = c(0,4.2008),unstable.eq.y = c(0,4.0039), x.bound = xbounds, y.bound = ybounds)


##### 0.5 quasi-potential vizualization!!! #####
	QPContour(surface=eq.1.global,dens=c(1000,1000), x.bound = xbounds, y.bound = ybounds,c.parm=5)


##### 0.6 vector field decompisition!!! #####
	# 0.6.0 all fields
	VDAll <- VecDecompAll(surface=eq.1.global, x.rhs=testequationx, y.rhs=testequationy, x.bound=xbounds, y.bound=ybounds)

	# 0.6.1 vector field
	VDV <- VecDecompVec(x.num.steps=step.number.x, y.num.steps=step.number.y, x.rhs=testequationx, y.rhs=testequationy, x.bound=xbounds, y.bound=ybounds)
		VecDecompPlot(field=list(VDV[,,1],VDV[,,2]), dens=c(50,50), x.bound=xbounds, y.bound=ybounds, x.lim=c(0,11), y.lim=c(0,6), arrow.type="proportional", tail.length=0.75, length=0.03)

	# 0.6.2 gradient field	
	VDG <- VecDecompGrad(eq.1.global)
		VecDecompPlot(field=list(VDG[,,1],VDG[,,2]), dens=c(50,50), x.bound=xbounds, y.bound=ybounds, arrow.type="proportional", length=0.03, tail.length=0.5)

	# 0.6.3 remainder field
	VDR <- VecDecompRem(surface=eq.1.global, x.rhs=testequationx, y.rhs=testequationy, x.bound=xbounds, y.bound=ybounds)
		VecDecompPlot(field=list(VDR[,,1],VDR[,,2]), dens=c(50,50), x.bound=xbounds, y.bound=ybounds, arrow.type="proportional", tail.length=0.5, length=0.03)

##### 0.7 3D graphs!!! #####
	dens.sub <- c(1000,1000)
	global.sub <- eq.1.global[round(seq(1,nrow(eq.1.global),length.out=dens.sub[1])) , round(seq(1,ncol(eq.1.global),length.out=dens.sub[2]))]
	global.sub[global.sub > 0.1] <- NA
	persp(global.sub, xlab="", ylab="y", zlab="", theta=0, phi=42.5, d=5, expand=.6, zlim=c(0,0.11), xlim=c(0,.4), ylim=c(0.15,.3), lphi=30, ltheta=0, col="orange", shade=1.25, ticktype="detailed", border=NA, r=1)