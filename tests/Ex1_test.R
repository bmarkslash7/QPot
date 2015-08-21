##### 0.0 prepreperation #####
	# install.packages("fortunes")
	# library("fortunes")
	# fortune(344) # or use fortune() to be distracted for way too long

##### 0.1 preperation #####
	# clean up and set seed
	rm(list=ls())
	your.favourite.number <- 
	set.seed(your.favourite.number) # see fortune(306)
	def.par <- par() # for deefault par settings for resetting plottting functions

	#load libraries

##### 0.2 stochastic simulations!!! #####
	# model function (same as other DE packages)
		model.ex1 <- function(t, state, parms) {
		with(as.list(c(state, parms)), {
		dx <- (alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))
		dy <- ((gamma*(x^2)*y)/(kappa + (x^2))) - mu*(y^2)
		list(c(dx,dy))
		})
		}

	# parameters
		model.state <- c(x = 1 , y = 2)
		model.parms <- c(alpha = 1.54 , beta = 10.14 , delta = 1 , kappa = 1 , gamma = 0.476 , mu = 0.112509)
		model.sigma <- 0.05
		model.time <- 1000
		model.deltat <- 0.2

	# time series
		ts.out.ex1 <- TSTraj(y0= model.state,time=model.time,deltat=model.deltat,func=model.ex1,parms=model.parms,sigma=model.sigma)
			# you can futz around with the emat
				apply(ts.out.ex1[,2:3],2,mean)

	# time series plots
		TSPlot(ts.out.ex1,deltat=model.deltat)
		par(def.par)
		TSPlot(ts.out.ex1,deltat=model.deltat,dim=2)
		TSDensity(ts.out.ex1,dim=1)
		TSDensity(ts.out.ex1,dim=2)

##### 0.3 local quasi-potential!!! #####

##### 0.4 global quasi-potential!!! #####

	e1.global <- QPGlobal(list(e1.1.local,e1.2.local),c(0,4.2008),c(0,4.0039),c(-0.5,20),c(-0.5,20))

##### 0.5 quasi-potential vizualization!!! #####

	QPContour(e1.global,c(1000,1000),c(-0.5,20),c(-0.5,20),c=5)

##### 0.6 vector field decompisition!!! #####

	VecDecompAll(e1.global, equations,c(-0.5,20),c(-0.5,20))

	# all fields
	VDAll <- VecDecompAll(e1.global, equations,c(-0.5,20),c(-0.5,20))

	# vector field
	VDV <- VecDecompVec(4100,4100, equations,c(-0.5,20),c(-0.5,20))
		VecDecompPlot(list(VDV[,,1],VDV[,,2]),density=c(50,50),x.bound=c(-0.5,20),y.bound=c(-0.5,20))

	# gradient field	
	VDG <- VecDecompGrad(e1.global)
		VecDecompPlot(list(VDG[,,1],VDG[,,2]),density=c(50,50),x.bound=c(-0.5,20),y.bound=c(-0.5,20))

	# remainder field
	VDR <- VecDecompRem(e1.global, equations,c(-0.5,20),c(-0.5,20))
		VecDecompPlot(list(VDR[,,1],VDR[,,2]),density=c(50,50),x.bound=c(-0.5,20),y.bound=c(-0.5,20))
