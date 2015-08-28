#' Wrapper for call to quasipotential computation using the upwind ordered method
#' 
#' 
#' @param x.rhs A string containing the right hand side of the equation for x.
#' @param x.start The starting value of x, usually the x value of the current equilibrium.
#' @param x.bound The x boundaries denoted at c(minimum, maximum).
#' @param x.num.steps The number of steps between the minimum and maximum x value defined in x range.
#' @param y.rhs A string containing the right hand side of the y equation.
#' @param y.start The starting value of y, usually the y value of the current equilibrium.
#' @param y.bound The y boundaries denoted at c(minimum, maximum).
#' @param y.num.steps The number of steps between the minimum and maximum y value defined in y range.
#' @param save.to.R Output the matrix of results for the upwind-ordered method to the current R session.  The default is not to write the matrix to the R session.  save.to.R=TRUE writes the output matrix to the R session.
#' @param save.to.HD Write the matrix of results for the upwind-ordered method to the hard drive named filename.  Default is TRUE.
#' @param filename String for the name of the file saved to the hard drive.  If filename is left blank, output file saved as defaultname-xX.STARTyY.START.txt, where X.START and Y.START are values in x.start and y.start, respectively. 
#' @param bounce By default, the upwind-ordered method stops when the boundaries are reached.  The bounce parameter allows the (d)efault action, only (p)ositive values to be tested, or reflection near the boundaries (bounce = 'b').
#' @param bounce.edge If bounce = 'b', then to prevent the upwind-ordered method from reaching the boundaries, temporary boundaries are created inside the boundaries defined by x.bound and y.bound.  The boundary edge is bounce.edge of the total range.  Default is 0.01
#' @param verboseR NOT IMPLEMENTED: Flag (default = FALSE) for printing out information in QPotential Rwrapper
#' @param verboseC NOT IMPLEMENTED: Flag (default = FALSE) for printing out useful-for-everyone information in quasipotential.C
#' @param debugC NOT IMPLEMENTED: Flag (default = FALSE) for printing out debugging C code 
#' @return filetoHD If save.to.HD enabled, then saves a file in the current directory as either filename or as defaultname-xXSTARTyYSTART.txt
#' @return filetoR If save.to.R enabled, then the function QPotential returns a matrix containing  the upwind-ordered results to be used for plotting.  Requires a variable to catch the returned matrix, i.e. storage <- QPotential(parameters...)
#'
#' @examples
#' #Example 1, Equilibrium 1 from article
#' #Equations
#' equationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
#' equationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
#' #Boundaries of the x and y state space
#' xbounds = c(-0.5, 20.0)
#' ybounds = c(-0.5, 20.0)
#' # Step number
#' # In the examples, these are 4100
#' # set to a smaller number for this example
#' xstepnumber = 1000
#' ystepnumber = 1000
#' # Equilibrium x and y values for Equilibrium 1 in Example 1
#' xinit1 = 1.40491
#' yinit1 = 2.80808
#' storage.eq1 <- QPotential(x.rhs = equationx, x.start = xinit1, 
#'  x.bound = xbounds, x.num.steps = xstepnumber, 
#'  y.rhs = equationy, y.start = yinit1, 
#'  y.bound = ybounds, y.num.steps = ystepnumber)

QPotential <- function (x.rhs = 'NULL', x.start = 'NULL', x.bound = 'NULL', x.num.steps = 'NULL', y.rhs = 'NULL', y.start = 'NULL', y.bound = 'NULL', y.num.steps = 'NULL', filename = 'NULL', save.to.R = TRUE, save.to.HD = FALSE, bounce = 'd', bounce.edge = 0.01, verboseR = FALSE, verboseC = FALSE, debugC = FALSE)
{
# ----------------------------------------------------------------------
# Break apart function parameters into things that C code will use
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# check if any component of x is missing
# ----------------------------------------------------------------------
if (verboseR) {print("check if any component of x is missing")}
if (x.rhs == 'NULL') {stop("No equation defined for x. Define x.rhs.")} else {equationx <- x.rhs}
if (grepl(pattern="=", x = x.rhs)) {stop("Equals sign (=) found in x.rhs parameter.  Please give only right hand side of the equation.")}
if (x.start == 'NULL') {stop("No starting value for x. Define x.start.")} else {startxval <- x.start}
if (x.num.steps == 'NULL') {stop("Need the number of steps in x range. Define x.num.steps.")} else {numofstepsx <- x.num.steps}
lengthequationx	<- nchar(x.rhs)
if (length(x.bound) < 2) stop("Not enough values for x range in parameter x.bound.")
if (length(x.bound) > 2) stop("Too many values for x range in parameter x.bound.")
lowerboundsx	<- x.bound[1]
upperboundsx	<- x.bound[2]

# ----------------------------------------------------------------------
# check if any component of y is missing
# ----------------------------------------------------------------------
if (verboseR) {print("check if any component of y is missing")}
if (y.rhs == 'NULL') {equationy = '0'} else {equationy <- y.rhs}
if (grepl(pattern="=", x = y.rhs)) {stop("Equals sign (=) found in y.rhs parameter.  Please give only right hand side of the equation.")}
lengthequationy <- nchar(y.rhs)
if ( (length(y.bound) == 1) && (y.bound == 'NULL') ) {lowerboundsy <- 0; upperboundsy <- 0}
if (y.bound[1] != 'NULL') {
	if (length(y.bound) < 2) stop("Not enough values for y range in variable y.bound.")
	if (length(y.bound) > 2) stop("Too many values for y range in variable y.bound.")
	lowerboundsy <- y.bound[1]; upperboundsy	<- y.bound[2]
} # end of if y.bound != 'NULL 
if (y.bound[1] == 'NULL') stop('No minimum and maximum y values.  Parameter y.bound not defined')
#The numofstepsy cannot equal 1, because hy=(LY2-LY1)/(NY-1), where NY is numofsteps
if (y.num.steps == 'NULL') {numofstepsy <- 2} else {numofstepsy <- y.num.steps}
if (y.start == 'NULL') {startyval = 0} else {startyval <- y.start}

#TODO only create this storage array if save.to.R == TRUE
# I think the .C function creates its own storage matrix, despite what 
# I tell it or what it is supposed to do.
# storage <- array(1.0, dim=c(1,(numofstepsx*numofstepsy)))

# ----------------------------------------------------------------------
# Components to enable file saving
# ---------------------------------------------------------------------- 
if (verboseR) {print("components to enable file saving")}
#Save in whatever format the user wants
if (filename != 'NULL') {lengthfilename = nchar(filename)}
else 					{lengthfilename = 0}

#default filename has a restriction on the name size
if ( (filename == 'NULL') && (save.to.HD == TRUE) && ( (abs(x.start) > 99999) || (abs(y.start) > 99999) ) ) { stop('Cannot use default filename because program will crash.  Please supply filename')}

if ((save.to.R == TRUE) && (save.to.HD == TRUE)) 		{datasave = 3}
else if ((save.to.R == TRUE) && (save.to.HD == FALSE)) 	{datasave = 2}
else if (isTRUE(save.to.HD))							{datasave = 1}
else												{datasave = 4}
if (save.to.HD == 'testrun')							{datasave = 4}
if (verboseR) {print(paste("Variable datasave is: ", datasave, sep = ""))}


# ----------------------------------------------------------------------
# Determine what C code does at edges of x.bound and y.bound
# bounce (used to make bounce.style), bounce.edge
# ----------------------------------------------------------------------
if (verboseR) {print("Determine what C code does at edges of x.bound and y.bound")}
if (!is.numeric(bounce.edge)) stop('Parameter bounce.edge must be a number')

if (bounce == FALSE) {bounce.style = 'd'}
else if (bounce == 'd') {bounce.style = 'd'}
else if ((bounce == 'p') || (bounce == 'b')) {bounce.style = bounce}
else 	{stop('Parameter bounce must be left blank, (d)efault), (p)ositivevalues, or (b)ounce')}

# ----------------------------------------------------------------------
# Error checking before C code is called
# ----------------------------------------------------------------------

{ 
# check to make sure:
# check that LX1 < xeq < LX2
# check that LY1 < yeq < LY2
if ( (lowerboundsx > startxval) || (upperboundsx < startxval) ) {stop("Starting x value x.start is outside x.bound range")}
if ( (lowerboundsy > startyval) || (upperboundsy < startyval) ) {stop("Starting y value y.start is outside y.bound range")}
if ( lowerboundsx > upperboundsx ) {stop("In x.bound, upper bound is less than lower bound.  x.bound[2] < x.bound[1]")}
if ( lowerboundsy > upperboundsy ) {stop("In y.bound, upper bound is less than lower bounds. y.bound[2] < y.bound[1]")}
 

# ----------------------------------------------------------------------
# warn about memory size is number of steps is very large
# ----------------------------------------------------------------------
if (numofstepsx*numofstepsy > 7000*7000) {
# 8 bytes per numeric, C code and R code each has an array, total number of elements
	tempsizeofarray = 2*8*numofstepsx*numofstepsy
	if (tempsizeofarray > 1e+06) {tempsizeofarraystring <- paste(tempsizeofarray / 1e+06, "MB", sep = " ")}
	if (tempsizeofarray > 1e+09) {tempsizeofarraystring <- paste(tempsizeofarray / 1e+09, "GB", sep = " ")}
	warningstring <- paste("Program expected to use at least", tempsizeofarraystring, "of memory. \n If program crashes, reduce the number of steps and retry.", sep = " ")
	warning(warningstring)
} #end memory size check

} # end of error checking

# ----------------------------------------------------------------------
# Call to upwind ordered method
# ----------------------------------------------------------------------
# Produce multiple versions: 	
#	1) store output in R and on harddrive
#	2) store output only on harddrive
#	3) store output only in R
#	4) no data saved, testing purposes only
if (datasave == 1) {
	#no R write; HD write
	storage = 0;
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), PACKAGE="QPot")
#	print(ls())
	return(TRUE)
}
else if (datasave == 2) {
	#R write; no HD write
	storage <- array(1.0, dim=c(1,(numofstepsx*numofstepsy)))
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), PACKAGE="QPot")
	storage = out2[[1]]
	storage <- matrix(storage, nrow = x.num.steps, byrow = TRUE)
	#1.0e+6 is the INFTY place holder in the C code
	#it means that no QP value was computed
	tstorage = t(storage)
	tstorage[tstorage > ((1.0e+6) - 1)] = NA 
	rm(storage)
	return(tstorage)
}
else if (datasave == 3) {
	# R write; HD write
	storage <- array(1.0, dim=c(1,(numofstepsx*numofstepsy)))
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), PACKAGE="QPot")
	storage = out2[[1]]
	storage <- matrix(storage, nrow = x.num.steps, byrow = TRUE)
	#1.0e+6 is the INFTY place holder in the C code
	#it means that no QP value was computed
	tstorage = t(storage)
	tstorage[tstorage > ((1.0e+6) - 1)] = NA 
	rm(storage)
	return(tstorage)
}
else if (datasave == 4) {
	# no R write; no HD write
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bounce.style, as.double(bounce.edge), PACKAGE="QPot")
	return(TRUE)
}
else {print("datasave is not a possible number.  How did you get here?")}

} # end of QPotential
