# Below is a pipe every 10 characters.  The R Journal suggests that code is wrapped around 80 characters
#*********|*********|*********|*********|*********|*********|*********|*********
require(QPot)
require(rgl)

##### 0.1 preperation #####
	# 0.1.0 clean up and set seed
	rm(list=ls())
	your.favourite.number <- 3818919 #chris
	set.seed(your.favourite.number)


##### 0.2 stochastic simulations!!! #####
	# 0.2.0 equations without parameters for easier paramter manipulation in a list (in 0.2.1)
		var.eqn.x <- "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa+(x^2)))"
		var.eqn.y <- "((gamma*(x^2)*y)/(kappa+(x^2))) - mu*(y^2)"

	# 0.2.1 parameters
		model.state <- c(x = 1, y = 2)
		model.parms <- c(alpha = 1.54, beta = 10.14, delta = 1, gamma = 0.476, kappa = 1, mu = 0.112509)
		model.sigma <- 0.05
		model.time <- 12500
		model.deltat <- 0.025

	# 0.2.2 time series
	ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)

	# 0.2.3 time series plots
		TSPlot(ts.ex1, deltat = model.deltat)
		TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
		TSDensity(ts.ex1, dim = 1)
		TSDensity(ts.ex1, dim = 2)

		# plots for the paper figures
			# print.wd <- "/Users/christophermoore/DropBox/QPRPackage/QPotPaper/Figures/"
			# # Ex1_TS_1D.png
			# png(paste(print.wd,"Ex1_TS_1D.png",sep=""), width = 720, height = 480)
			# TSPlot(ts.ex1, deltat = model.deltat)
			# dev.off()
			# # Ex1_TS_2D.png
			# png(paste(print.wd,"Ex1_TS_2D.png",sep=""), width = 500, height = 500)
			# TSPlot(ts.ex1, deltat = model.deltat, dim = 2, line.alpha = 50, lwd = 1, xlab = expression(italic(x)), ylab = expression(italic(y)))
			# dev.off()
			# # Ex1_Dens_2D.png
			# png(paste(print.wd,"Ex1_Dens_2D.png",sep=""), width = 500, height = 500)
			# TSDensity(ts.ex1, deltat = model.deltat, dim = 2, line.alpha = 50, lwd = 1, xlab = expression(italic(x)), ylab = expression(italic(y)), contour.lwd = 0.75)
			# dev.off()


##### 0.3 local quasi-potential!!! #####
	equation.x = "1.54*x*(1.0-(x/10.14))-(y*x*x)/(1.0+x*x)"
	equation.y = "((0.476*x*x*y)/(1+x*x))-0.112590*y*y"
	bounds.x = c(-0.5, 20.0)
	bounds.y = c(-0.5, 20.0)
	step.number.x = 4100
	step.number.y = 4100
	eq1.x = 1.40491
	eq1.y = 2.80808
	eq2.x = 4.9040
	eq2.y = 4.06187

	eq1.local <- QPotential(x.rhs = equation.x, x.start = eq1.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq1.y,  y.bound = bounds.y, y.num.steps = step.number.y)
	eq2.local <- QPotential(x.rhs = equation.x, x.start = eq2.x, x.bound = bounds.x, x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq2.y, y.bound = bounds.y, y.num.steps = step.number.y)


##### 0.4 global quasi-potential!!! #####
	ex1.global <- QPGlobal(local.surfaces = list(eq1.local, eq2.local), unstable.eq.x = c(0, 4.2008), unstable.eq.y = c(0, 4.0039), x.bound = bounds.x, y.bound = bounds.y)


##### 0.5 quasi-potential vizualization!!! #####
	QPContour(surface = ex1.global, dens = c(1000, 1000), x.bound = bounds.x, y.bound = bounds.y, c.parm = 5)
	# A paper-like graph
		# par(mfrow = c(1, 2), mar = c(4, 2, 0.5, 0.5), oma = rep(3,4))
		# QPContour(surface = ex1.global, dens = c(100, 100), x.bound = bounds.x, y.bound = bounds.y, c.parm = 5, xlab = "", ylab = "", lwd = 0.5) #in the paper [ylim = c(0, 5), xlim = c(0, 8),]
		# mtext(expression(italic(X)), side = 1, line = 2.5)
		# mtext(expression(italic(Y)), side = 2, line = 2.5 , las = 1)
		# QPContour(surface = ex1.global, dens = c(100, 100), x.bound = bounds.x, y.bound = bounds.y, c.parm = 5, xlab = "", ylab = "", lwd = 0.5) #in the paper [ylim = c(0, 5), xlim = c(0, 8),]
		# mtext(expression(italic(X)), side = 1, line = 2.5)

	# figure for the paper
	# print.wd <- "/Users/christophermoore/DropBox/QPRPackage/QPotPaper/Figures/"
	# png(paste(print.wd,"Ex1_QP_contour.png",sep="") , width = 800 , height = 400)
	# par(mfrow=c(1,2) , mar = c(1, 1, 1, 1) , oma = c(3,3,0,0))
	# QPContour(ex1.global,c(1000,1000),c(-0.5,20),c(-0.5,20),c.parm=1,lwd=0.5 , xlab = "" , ylab = "" ,n.contour.lines=25)
	# mtext(expression(italic(y)),2,line=2.5,las=1)
	# mtext(expression(italic(x)),1,line=2.5)
	# QPContour(ex1.global,c(1000,1000),c(-0.5,20),c(-0.5,20),c.parm=4.5,lwd=0.5 , ylab = "" , xlab = "",n.contour.lines=25)
	# mtext(expression(italic(X)),1,line=2.5)
	# dev.off()

##### 0.6 vector field decompisition!!! #####
	# 0.6.0 all fields
	VDAll <- VecDecomAll(surface = ex1.global, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
	VecDecomPlot(field = list(VDAll[,,1], VDAll[,,2]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, x.lim = c(0, 11), y.lim = c(0, 6), arrow.type = "equal", tail.length = 0.25, head.length = 0.025)
	VecDecomPlot(field = list(VDAll[,,3], VDAll[,,4]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.25, head.length = 0.025)
	VecDecomPlot(field = list(VDAll[,,5], VDAll[,,6]), dens = c(25, 25), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.35, head.length = 0.025)


	# 0.6.1 vector field
	VDV <- VecDecomVec(x.num.steps = step.number.x, y.num.steps = step.number.y, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
		VecDecomPlot(field = list(VDV[,,1], VDV[,,2]), dens = c(50, 50), x.bound = bounds.x, y.bound=bounds.y, x.lim = c(0, 11), y.lim=c(0, 6), arrow.type="proportional", tail.length=0.75, head.length=0.03)

	# 0.6.2 gradient field	
	VDG <- VecDecomGrad(ex1.global)
		VecDecomPlot(field = list(VDG[,,1], VDG[,,2]), dens=c(50, 50), x.bound=bounds.x, y.bound = bounds.y, arrow.type = "proportional", head.length = 0.03, tail.length = 0.5)

	# 0.6.3 remainder field
	VDR <- VecDecomRem(surface = ex1.global, x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, y.bound = bounds.y)
		VecDecomPlot(field = list(VDR[,,1], VDR[,,2]), dens = c(50, 50), x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", tail.length = 0.5, head.length = 0.03)

	# plots for the paper figure
		# print.wd <- "/Users/christophermoore/DropBox/QPRPackage/QPotPaper/Figures/"
			# # Ex1_VecDecom.png
		# png(paste(print.wd,"Ex1_VecDecom.png",sep=""), width = 720, height = 720)
		# par(mfrow=c(2,2),mar = c(0.5,0.5,0,0),oma=c(2.5,3.5,2,1))
# VecDecomPlot(field = list(VDAll[,,3], VDAll[,,4]), dens = c(25, 25), x.bound = c(-0.5, 20), y.bound=c(-0.5,20),head.length=0.02,xlab="",ylab="",arrow.type="proportional",tail.length=0.25,xaxt="n",yaxt="n",lwd=0.9)
		# axis(1,labels=F,tck=-0.015)
		# axis(2,tck=-0.015,labels=F,las=1)
		# mtext(1:4,2,at=1:4,line=0.5,las=1,cex=0.8)
		# mtext(expression(italic(y)),2,line=1.75,las=1,cex=0.8)
		# mtext("Proportional arrow lengths",2,line=2.75)
		# mtext(expression(Gradient~field~(-nabla~phi(x,y))),3,line=0.1,las=1)

	# VecDecomPlot(field = list(VDAll[,,5], VDAll[,,6]), dens = c(25, 25), x.bound = c(-0.5, 20), y.bound = c(-0.5, 20), head.length = 0.02, xlab = "", ylab = "", arrow.type = "proportional", tail.length = 0.25, xaxt = "n", yaxt = "n", lwd = 0.9)
		# axis(1,labels=F,tck=-0.015)
		# axis(2,labels=F,tck=-0.015)
		# mtext(expression(Remainder~field~(bold(r)(x,y))),3,line=0.1,las=1)

	# VecDecomPlot(field = list(VDAll[,,3], VDAll[,,4]), dens = c(25, 25), x.bound = c(-0.5,20), y.bound = c(-0.5, 20), head.length = 0.02, xlab = "", ylab = "", arrow.type = "equal", tail.length = 0.2, lwd = 0.9, xaxt = "n", yaxt = "n")
		# axis(1,labels=F,tck=-0.015)
		# axis(2,tck=-0.015,labels=F,las=1)
		# mtext(1:4,2,at=1:4,line=0.5,las=1,cex=0.8)
		# mtext(seq(0,8,2),1,at=seq(0,8,2),line=0.5,las=1,cex=0.8)
		# mtext(expression(italic(y)),2,line=1.75,las=1,cex=0.8)
		# mtext(expression(italic(x)),1,line=1.75,cex=0.8)
		# mtext("Equal arrow lengths",2,line=2.75)

	# VecDecomPlot(field = list(VDAll[,,5], VDAll[,,6]), dens = c(25, 25), x.bound = c(-0.5, 20), y.bound = c(-0.5, 20), head.length = 0.02, xlab="", ylab = "", arrow.type = "equal", tail.length = 0.2, yaxt ="n", xaxt = "n", lwd = 0.9)
		# axis(1,labels=F,tck=-0.015)
		# axis(2,labels=F,tck=-0.015)
		# mtext(seq(0,8,2),1,at=seq(0,8,2),line=0.5,las=1,cex=0.8)
		# mtext(expression(italic(x)),1,line=1.75,cex=0.8)
		# dev.off()
##### 0.7 3D graphs!!! #####
	library(rgl)
	dens.sub <- c(4000,4000)
	global.sub <- ex1.global[round(seq(1,nrow(ex1.global),length.out=dens.sub[1])) , round(seq(1,ncol(ex1.global),length.out=dens.sub[2]))]
	# global.sub[global.sub > 0.02] <- NA # to cut off large values for a better view of the basin
	persp3d(x = global.sub, col = "orange", expand = 1.1, xlim = c(0.05, 0.35), ylim = c(0.1, 0.3), zlim = c(0, 0.01), xlab = "X", ylab = "Y", zlab = intToUtf8(0x03A6))
	# rgl.snapshot("/Users/christophermoore/Dropbox/QPRpackage/QPotPaper/Figures/Ex1_3D.png", fmt = "png", top = TRUE ) #to print view as a png
