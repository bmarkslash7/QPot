##### 0.0 prepreperation #####
	# install.packages("fortunes")
	# library("fortunes")
	# fortune(344) # or use fortune() to be distracted for way too long

##### 0.1 preperation #####
	# 0.1.0 clean up and set seed
	rm(list=ls())
	your.favourite.number <- 6174
	set.seed(your.favourite.number) # see fortune(306)
	def.par <- par() # for deefault par settings for resetting plottting functions

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
		par(def.par)
		TSPlot(ts.out.ex1, deltat=model.deltat, dim=2)
		TSDensity(ts.out.ex1, dim=1)
		TSDensity(ts.out.ex1, dim=2)

##### 0.3 local quasi-potential!!! #####
		testequationx = "1.54*x*(1.0-(x/10.14))-(y*x*x)/(1.0+x*x)"
		testequationy = "((0.476*x*x*y)/(1+x*x))-0.112590*y*y"
		xbounds = c(-0.5, 20.0)
		ybounds = c(-0.5, 20.0)
		xstepnumber = 4100
		ystepnumber = 4100
		xinit1 = 1.40491
		yinit1 = 2.80808
		xinit2 = 4.9040
		yinit2 = 4.06187

		QPotential(x.rhs = testequationx, x.start = xinit1, x.bound = xbounds, x.num.steps = xstepnumber, 
			y.rhs = testequationy, y.start = yinit1, y.bound = ybounds, y.num.steps = ystepnumber, 
			filename = 'exampleone_eq1.txt')
		QPotential(x.rhs = testequationx, x.start = xinit2, x.bound = xbounds, x.num.steps = xstepnumber, 
			y.rhs = testequationy, y.start = yinit2, y.bound = ybounds, y.num.steps = ystepnumber, 
			filename = 'exampleone_eq2.txt')
	# 0.3.0 read in the results using fread() in the package data.table
		require(data.table)		#		requireNamespace("data.table")
		e1.1.raw <- fread('exampleone_eq1.txt')
		e1.1.local <- t(data.matrix(e1.1.raw))
		rm(e1.1.raw)
		e1.2.raw <- fread('exampleone_eq2.txt')
		e1.2.local <- t(data.matrix(e1.2.raw))
		rm(e1.2.raw)
#	#read in the tables using read.table to test for speed and compare to fread + transpose
#		e1.1.readtest <- read.table('exampleone_eq1.txt', sep="\t", header = FALSE)
#		e1.2.readtest <- read.table('exampleone_eq2.txt', sep="\t", header = FALSE)

		e1.1.local[e1.1.local == 5e+05] = NA
		e1.2.local[e1.2.local == 5e+05] = NA

		if (isTRUE(all.equal(t(data.matrix(e1.1.readtest)), e1.1.local))) {print("SUCCESS - Code from e1.1 get same answer")} else {print("FAIL - code from e1.1 produced different matrices")} #PASSES
		if (isTRUE(all.equal(t(data.matrix(e1.2.readtest)), e1.2.local))) {print("SUCCESS - Code from e1.2 get same answer")} else {print("FAIL - code from e1.2 produced different matrices")} #PASSES

##### 0.4 global quasi-potential!!! #####

	e1.global <- QPGlobal(local.surfaces = list(e1.1.local,e1.2.local),unstable.eq.x = c(0,4.2008),unstable.eq.y = c(0,4.0039), x.bound = c(-0.5,20),y.bound = c(-0.5,20))

##### 0.5 quasi-potential vizualization!!! #####

	QPContour(e1.global,c(1000,1000),c(-0.5,20),c(-0.5,20),c.parm=5)

##### 0.6 vector field decompisition!!! #####

	# 0.6.0 all fields
	VDAll <- VecDecompAll(surface=e1.global, x.rhs=testequationx, y.rhs=testequationy ,x.bound=c(-0.5,20), y.bound=c(-0.5,20))

	# 0.6.1 vector field
	VDV <- VecDecompVec(x.num.steps=4100, y.num.steps=4100, x.rhs=testequationx, y.rhs=testequationy, x.bound=c(-0.5,20), y.bound=c(-0.5,20))
		VecDecompPlot(field=list(VDV[,,1],VDV[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))

	# 0.6.2 gradient field	
	VDG <- VecDecompGrad(e1.global)
		VecDecompPlot(field=list(VDG[,,1],VDG[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))

	# 0.6.3 remainder field
	VDR <- VecDecompRem(surface=e1.global, x.rhs=testequationx, y.rhs=testequationy, x.bound=c(-0.5,20), y.bound=c(-0.5,20))
		VecDecompPlot(field=list(VDR[,,1],VDR[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))
