# An R package for stochastic differential equation quasi-potential analysis
[![cran version](https://www.r-pkg.org/badges/version/QPot)]( https://CRAN.R-project.org/package=QPot)
[![rstudio mirror per-month downloads](https://cranlogs.r-pkg.org/badges/QPot)](https://github.com/r-hub/cranlogs.app)
[![rstudio mirror total downloads](https://cranlogs.r-pkg.org/badges/grand-total/QPot?color=yellowgreen)](https://github.com/r-hub/cranlogs.app)

### Christopher M. Moore, Christopher R. Stieha, Ben C. Nolting, Maria K. Cameron, and Karen C. Abbott

QPot offers a range of tools to simulate, analyze, and visualize the dynamics of two-dimensional systems of stochastic differential equations.  QPot offers tools to compute the quasi-potential, which is useful when comparing the relative stabilities of equilibria in systems with multiple stable equilibria. 

### Literature: ###

M. K. Cameron. 2012. Finding the quasipotential for nongradient SDEs. *Physica D*, 241(18):1532–1550.

Moore, C.M., Stieha, C.R., Nolting, B.C., Cameron, M.K., and Abbott, K.C. *Resubmitted to The R Journal.* QPot: An R package for stochastic differential equation quasi-potential analysis

B. C. Nolting and K. C. Abbott. 2016. Balls, cups, and quasi-potentials: quantifying stability in stochastic systems. *Ecology,* 97(4):850-864.

### Functions ###

* Model2String()	converts equations in function format to strings
* QPContour()		creates a contour plot of the quasi-potential
* QPGlobal()		creates a global quasi-potential surface
* QPInterp()		evaluates the quasi-potential at point (x,y)
* QPotential()		computes the quasi-potential 
* TSDensity()		creates a density plot of population trajectories from TSTraj()
* TSPlot()			plots population trajectories from TSTraj()
* TSTraj()			simulates a time series based off the equations
* VecDecomAll()		returns three vector fields: VecDecomGrad(), VecDecomRem(), VecDecompVec()
* VecDecomGrad()	vector field of the negative gradient of the quasi-potential
* VecDecomRem()		vector field of the remainder
* VecDecomPlot()	plots the vector field
* VecDecomVec()		vector field of deterministic skeleton

### Notation for mathematical equations ###

*The library expression_parser (https://github.com/jamesgregson/expression_parser) is used under the GPLv2 for non-commerical use.*

Reads strings of mathematical expression in infix notation.  
* standard arithmetic operations (+,-,*,/) with operator precedence
* exponentiation ^ and nested exponentiation
* unary + and -
* expressions enclosed in parentheses ('(',')'), optionally nested
* built-in math functions: pow(x,y), sqrt(x), log(x), exp(x), sin(x), asin(x), cos(x), acos(x), tan(x), atan(x), atan2(y,x), abs(x), fabs(x), floor(x), ceil(x), round(x), with input arguments checked for domain validity, e.g. 'sqrt( -1.0 )' returns an error.

See http://jamesgregson.blogspot.com/2012/06/mathematical-expression-parser-in-c.html.

### License ###
 
*This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 2 of the License.*

### Example 1 from article submitted to R Journal ###

With the equations and parameters, there are two stable equilibria.  This example simulates the time series (0.2.2 Time series), computes the local quasi-potential around the two equilibria (0.3 Compute the local quasi-potentials) and combines the two local quasi-potential to make the global quasi-potential (0.4 global quasi-potential).  These results can be viewed using a variety of functions (section 0.6 to 0.7). 

#### 0.2 Initialization ####
```R
require(QPot)

var.eqn.x <- "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa+(x^2)))"
var.eqn.y <- "((gamma*(x^2)*y)/(kappa+(x^2))) - mu*(y^2)"
model.parms <- c(alpha = 1.54, beta = 10.14, delta = 1, gamma = 0.476, 
	kappa = 1, mu = 0.112509)
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms)
parms.eqn.y <- Model2String(var.eqn.y, parms = model.parms, 
	supress.print = TRUE) # does not print to screen
model.state <- c(x = 1, y = 2)
model.sigma <- 0.05
model.time <- 1000 # we used 12500 in the figures
model.deltat <- 0.025
```
#### 0.2.2 Time series ####
```R
ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, 
     x.rhs = parms.eqn.x, y.rhs = parms.eqn.y, sigma = model.sigma)

# Could also use TSTraj to combine equation strings and parameter values
#ts.ex1 <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, 
#     x.rhs = var.eqn.x, y.rhs = var.eqn.y, parms = model.parms, sigma = model.sigma)

```
#### 0.2.3 Time series plots ####
```R
TSPlot(ts.ex1, deltat = model.deltat)
TSPlot(ts.ex1, deltat = model.deltat, dim = 2)
TSDensity(ts.ex1, dim = 1)
TSDensity(ts.ex1, dim = 2)
```
#### 0.3 Compute the local quasi-potentials ####
```R
parms.eqn.x <- Model2String(var.eqn.x, parms = model.parms)
   # if not done in a previous step
parms.eqn.x <- Model2String(var.eqn.y, parms = model.parms, 
   supress.print = TRUE) # does not print to screen
   
# Could also input the values by hand and use this version
# parms.eqn.x = "1.54*x*(1.0-(x/10.14)) - (y*(x^2))/(1.0+(x^2))"
# parms.eqn.x = "((0.476*(x^2)*y)/(1+(x^2))) - 0.112509*(y^2)"

bounds.x = c(-0.5, 20.0)
bounds.y = c(-0.5, 20.0)
step.number.x = 4100
step.number.y = 4100

eq1.x = 1.40491
eq1.y = 2.80808
eq2.x = 4.9040
eq2.y = 4.06187

eq1.local <- QPotential(x.rhs = equation.x, x.start = eq1.x, x.bound = bounds.x, 
				x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq1.y, 
				y.bound = bounds.y, y.num.steps = step.number.y)

eq2.local <- QPotential(x.rhs = equation.x, x.start = eq2.x, x.bound = bounds.x, 
				x.num.steps = step.number.x, y.rhs = equation.y, y.start = eq2.y, 
				y.bound = bounds.y, y.num.steps = step.number.y)
```

#### 0.4 Global quasi-potential ####
```R
ex1.global <- QPGlobal(local.surfaces = list(eq1.local, eq2.local), 
						unstable.eq.x = c(0, 4.2008), unstable.eq.y = c(0, 4.0039), 
						x.bound = bounds.x, y.bound = bounds.y)
```

#### 0.5 Quasi-potential visualization ####
```R
QPContour(surface = ex1.global, dens = c(1000, 1000), 
			x.bound = bounds.x, y.bound = bounds.y, c.parm = 5)
```

#### 0.6 Vector field decomposition ####
```R
VDAll <- VecDecomAll(surface = ex1.global, x.rhs = parms.eqn.x, 
					y.rhs = parms.eqn.y, x.bound = bounds.x, y.bound = bounds.y)

VecDecomPlot(x.field = VDAll[,,1], y.field = VDAll[,,2], dens = c(25, 25), 
			x.bound = bounds.x, y.bound = bounds.y, xlim = c(0, 11), 
			ylim = c(0, 6), arrow.type = "equal", tail.length = 0.25, 
			head.length = 0.025)

VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25, 25), 
			x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", 
			tail.length = 0.25, head.length = 0.025)

VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25, 25), 
			x.bound = bounds.x, y.bound = bounds.y, arrow.type = "proportional", 
			tail.length = 0.35, head.length = 0.025)
```

#### 0.6.1 vector field ####
```R
VDV <- VecDecomVec(x.num.steps = step.number.x, y.num.steps = step.number.y, 
					x.rhs = equation.x, y.rhs = equation.y, x.bound = bounds.x, 
					y.bound = bounds.y)

VecDecomPlot(field = list(VDV[,,1], VDV[,,2]), dens = c(50, 50), 
				x.bound = bounds.x, y.bound=bounds.y, xlim = c(0, 11), 
				ylim=c(0, 6), arrow.type="proportional", 
				tail.length=0.75, head.length=0.03)
```

#### 0.6.2 gradient field ####
```R
VDG <- VecDecomGrad(ex1.global)

VecDecomPlot(field = list(VDG[,,1], VDG[,,2]), dens=c(50, 50), 
				x.bound=bounds.x, y.bound = bounds.y, arrow.type = "proportional", 
				head.length = 0.03, tail.length = 0.5)
```

#### 0.6.3 remainder field ####
```R
VDR <- VecDecomRem(surface = ex1.global, x.rhs = equation.x, y.rhs = equation.y, 
					x.bound = bounds.x, y.bound = bounds.y)
```

#### 0.7 3D graphs ####
```R
require(rgl)

dens.sub <- c(4000,4000)

global.sub <- ex1.global[round(seq(1,nrow(ex1.global),length.out=dens.sub[1])),
						round(seq(1,ncol(ex1.global),length.out=dens.sub[2]))]

persp3d(x = global.sub, col = "orange", expand = 1.1, xlim = c(0.05, 0.35), 
		ylim = c(0.1, 0.3), zlim = c(0, 0.01), xlab = "X", ylab = "Y", 
		zlab = intToUtf8(0x03A6))
```
 


