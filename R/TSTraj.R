#' Simulate two-dimensional stochastic differential equations
#'
#' This function allows you to simulate two-dimensional stochastic differential equations.
#' @param y0 the a two-element vector of the initial conditions for the state variables.
#' @param time numeric value indicating the total time over which the simulation is to be run.
#' @param deltat numeric value indicating the frequency of stochastic perturbation, as \eqn{\Delta t}.
# @param func function containing deterministic equations formatted as \pkg{deSolve}.
#' @param x.rhs A string containing the right hand side of the equation for x.
#' @param y.rhs A string containing the right hand side of the y equation.
#' @param parms a named vector of paramters and their respective values for the deterministic equations.
#' @param sigma numberic value specifying the noise intensity.
#' @param lower.bound numeric value specifying a lower bound in the simulation.
#' @param upper.bound numeric value specifying an upper bound in the simulation.
#' @keywords Stochastic simulation
#' 
#' @examples
#' # First, the parameter values
#' state <- c(x = 3 , y = 3)
#' model.sigma <- 0.1
#' model.deltat <- 0.005
#'
#' # Second, write out the deterministic skeleton of the equations to be simulated
#' #Example 1 from article
#' equationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
#' equationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
#'
#' # Third, Run it
#' LVModelOut <- TSTraj(y0=state, time=250, deltat=model.deltat, 
#' x.rhs=equationx, y.rhs=equationy, sigma=model.sigma)
#'
#' #Can also input x.rhs and y.rhs as strings that contain parameter names
#' #and include parms with names and values of parameters
#' model.state <- c(x=1 , y=2)
#' model.parms <- c(alpha=1.54, beta=10.14, delta=1, kappa=1, gamma=0.476, mu=0.112509)
#' model.sigma <- 0.05
#' model.time <- 1000
#' model.deltat <- 0.2
#'
#' test.eqn.x = "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))"
#' test.eqn.y = "((gamma*(x^2)*y)/(kappa + (x^2))) - mu*(y^2)"
#'
#' ts.out.ex1 <- TSTraj(y0= model.state, time=model.time, deltat=model.deltat, 
#'   x.rhs=test.eqn.x, y.rhs= test.eqn.y, parms=model.parms, sigma=model.sigma)


TSTraj <- function(y0, time, deltat, x.rhs, y.rhs, parms = NA, sigma, lower.bound = NA, upper.bound = NA) {
	func <- function(t, state, parms) {
		with(as.list(c(state, parms)), {
		dx <- eval(parse(text=x.rhs))
		dy <- eval(parse(text=y.rhs))
		list(c(dx,dy))
		})
	}
	
	time.vals <- seq(from = 0 , to = time , by = deltat)[-1]
	mat <- matrix(data = NA , nrow = length(time.vals) , ncol = 3)
	colnames(mat) <- c("t" , names(y0))
	if(length(sigma) == 1) {
		sigma.a <- sigma
		sigma.b <- sigma
		} else {
			sigma.a <- sigma[1]
			sigma.b <- sigma[2]
		}
	if(length(upper.bound) == 1) {
		upper.bound.a <- upper.bound
		upper.bound.b <- upper.bound
		} else {
			upper.bound.a <- upper.bound[1]
			upper.bound.b <- upper.bound[2]
		}
	if(length(lower.bound) == 1) {
		lower.bound.a <- lower.bound
		lower.bound.b <- lower.bound
		} else {
			lower.bound.a <- lower.bound[1]
			lower.bound.b <- lower.bound[2]
		}

	if (is.numeric(upper.bound) == TRUE) {
			if (is.numeric(lower.bound) == TRUE) {
			mat[1,] <- c(1, y0[1], y0[2])
					for (i in 2:length(time.vals)) {
						mat[i,] <- i
						mat.list <- c(mat[i-1,2] ,mat[i-1,3])
						fX <- unlist(func(1 , mat.list , parms))[1]
						fY <- unlist(func(1 , mat.list , parms))[2]
						X <- mat.list[1] + fX*deltat + sigma.a*rnorm(1,0, sqrt(deltat))
						Y <- mat.list[2] + fY*deltat + sigma.b*rnorm(1,0, sqrt(deltat))
						mat[i,2] <- ifelse(X <= upper.bound.a, ifelse(X >= lower.bound.a,X,lower.bound.a), upper.bound.a)
						mat[i,3] <- ifelse(Y <= upper.bound.b, ifelse(Y >= lower.bound.b,Y,lower.bound.b), upper.bound.b)	
				}
			} else {
			mat[1,] <- c(1, y0[1], y0[2])
					for (i in 2:length(time.vals)) {
						mat[i,] <- i
						mat.list <- c(mat[i-1,2] ,mat[i-1,3])
						fX <- unlist(func(1 , mat.list , parms))[1]
						fY <- unlist(func(1 , mat.list , parms))[2]
						X <- mat.list[1] + fX*deltat + sigma.a*rnorm(1,0, sqrt(deltat))
						Y <- mat.list[2] + fY*deltat + sigma.b*rnorm(1,0, sqrt(deltat))
						mat[i,2] <- ifelse(X <= upper.bound.a, X, upper.bound.a)
						mat[i,3] <- ifelse(Y <= upper.bound.b, Y, upper.bound.b)
					}
				}
	
		} else {
			if (is.numeric(lower.bound) == TRUE) {
			mat[1,] <- c(1, y0[1], y0[2])
				for (i in 2:length(time.vals)) {
					mat[i,] <- i
					mat.list <- c(mat[i-1,2] ,mat[i-1,3])
					fX <- unlist(func(1 , mat.list , parms))[1]
					fY <- unlist(func(1 , mat.list , parms))[2]
					X <- mat.list[1] + fX*deltat + sigma.a*rnorm(1,0, sqrt(deltat))
					Y <- mat.list[2] + fY*deltat + sigma.b*rnorm(1,0, sqrt(deltat))
					mat[i,2] <- ifelse(X >= lower.bound.a, X, lower.bound.a)
					mat[i,3] <- ifelse(Y >= lower.bound.b, Y, lower.bound.b)
		
				}
			} else {
			mat[1,] <- c(1, y0[1], y0[2])
				for (i in 2:length(time.vals)) {
					mat[i,] <- i
					mat.list <- c(mat[i-1,2] ,mat[i-1,3])
					fX <- unlist(func(1 , mat.list , parms))[1]
					fY <- unlist(func(1 , mat.list , parms))[2]
					mat[i,2] <- mat.list[1] + fX*deltat + sigma.a*rnorm(1,0, sqrt(deltat))
					mat[i,3] <- mat.list[2] + fY*deltat + sigma.b*rnorm(1,0, sqrt(deltat))		
				}
			}
		}
		mat
	}
