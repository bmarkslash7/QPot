#' Function that converts differential equations from 
#' function-format to string-format
#' Specifically, it reads in a function, searches for the differential
#' equations within the function, and returns a list of strings containing
#' the differential equations from the function.
#' Will also replace parameters of those equations with numerical values.
#' 
#' @param model.function function containing the differential equations as given to TSTraj()
#' @param parms a named vector of paramters and their respective values for the deterministic equations.
#' @param x.lhs.term string containing the left hand side of the first equation to search for, default is 'dx'
#' @param y.lhs.term string containing the left hand side of the second equation to search for, default is 'dx'
#' @param supress.print Default it FALSE, supress output.  TRUE prints out equations from function
#' @return equations a list with two elements, the first is the x equation, the second is the y equation
#'
# @examples
# test.eqn.x = "(alpha*x)*(1-(x/beta)) - ((delta*(x^2)*y)/(kappa + (x^2)))"
# test.eqn.y = "((gamma*(x^2)*y)/(kappa + (x^2))) - mu*(y^2)"
# # 0.2.1 parameters
# model.state <- c(x=1 , y=2)
# model.parms <- c(alpha=1.54, beta=10.14, delta=1, kappa=1, gamma=0.476, mu=0.112509)
# equations.as.strings <- Model2String(

Model2String <- function(model.function, parms = 'NULL', x.lhs.term = 'dx', y.lhs.term = 'dy', supress.print = FALSE) {
	if (!supress.print) {
		print("Note: This function is supplied as duct tape.  Long equations, equations spanning multiple lines, equations with strange notation, etc, may not work.  Always check the output.")
	}
#	if (parms[1] == 'NULL') {stop("Need to define parms, the names and values of the model parameters")} 
	
	temp <- deparse(model.function, width.cutoff = 500)
#Dump function into a list of character strings
#Go through each string and determine if it contains an equation
	foundx = 0	#flag for making sure dx is only found once
	foundy = 0	#flag for making sure dy is only found once

#remove the lhs and return the rhs
	for (i in 1:length(temp)) {
	#when searching, first look for the lhs defining whether the derivative is for x or y
	#once found, look inside the string and use either '<-' or '=' to separate
	# the lhs from the rhs
		if ((foundx == 0) && isTRUE(grep(pattern = x.lhs.term, x = temp[i]) == 1)) {
			foundx = 1
			if (grep(pattern = '<-', x = temp[i]) == 1) {
				location <- regexpr(pattern = '<-', text = temp[i])
				x.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
			}
			else if (grep(pattern = '=', x = temp[i]) == 1) {
				location <- regexpr(pattern = '=', text = temp[i])
				x.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
			} else {stop("Equation does not contain = or <-")}
		}

		if ((foundy == 0) && isTRUE(grep(pattern = y.lhs.term, x = temp[i]) == 1)) {
			foundy = 1
			if (grep(pattern = '<-', x = temp[i]) == 1) {
				location <- regexpr(pattern = '<-', text = temp[i])
				y.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
			}
			else if (grep(pattern = '=', x = temp[i]) == 1) {
				location <- regexpr(pattern = '=', text = temp[i])
				y.equation = substr(temp[i], start = (location+2), stop = nchar(temp[i]))
			} else {stop("Equation does not contain = or <-")}
		}
	}
	
	equations = c(x.equation, y.equation)
#if parameters are not declared, then we do not have to replace anything
	if (!(parms[1] == 'NULL')){
#replace the parameter names in the equations with their values
		allnames <-names(parms)
		for (i in 1:length(parms)) {
			currname <- allnames[i]
			value = toString(parms[[i]])
			equations <- gsub(pattern = currname, replacement = value, x = equations)
		}
	} #if (!(parms[1] == 'NULL')){
	
	if (!supress.print) {print(paste(equations))}
	return(equations)

}
