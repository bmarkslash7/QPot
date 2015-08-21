#' Wrapper for call to quasipotential compution using the upwind ordered method
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
#' @param filename String for the name of the file saved to the hard drive.  If filename is left blank, output file saved as defaultname-xX.STARTyY.START.txt, where X.START and Y.START are values in x.start and y.start, respectively. 
#' @param save.to.R Output the matrix of results for the upwind-ordered method to the current R session.  The default is not to write the matrix to the R session.  save.to.R=TRUE writes the output matrix to the R session.
#' @param save.to.HD Write the matrix of results for the upwind-ordered method to the hard drive named filename.  Default is TRUE.
#' @param bounce Dy default, the upwind-ordered method stops when the boundaries are reached.  The bounce parameter allows the (d)efault action, only (p)ositive values to be tested, or reflection near the boundaries (bounce = 'b').
#' @param bounce.edge If bounce = 'b', then to prevent the upwind-ordered method from reaching the boundaries, temporary boundaries are created inside the boundaries defined by x.bound and y.bound.  The boundary edge is bounce.edge of the total range.  Default is 0.01
#' @return filetoHD If save.to.HD enabled, then saves a file in the current directory as either filename or as defaultname-xXSTARTyYSTART.txt
#' @return filetoR If save.to.R enabled, then the function QPotential returns a matrix containing  the upwind-ordered results to be used for plotting.  Requires a variable to catch the returned matrix, i.e. storage <- QPotential(parameters...)

# R CMD SHLIB -I/usr/local/include -L/usr/local/lib -lmatheval upwindorderedMATHEVALv3.c -lm
# -I adds directory to the head of the list of directories containing libraries
# -L add directory to be searched for that contains library listed in -l
# -l search the library
# -lm I don't know why or what.  google-fu failure

# File includes: 
# upwindordered -  the R wrapper to call the C code
# testupwindordered - checks to see if upwindordered throws the correct errors
#
# TODO - make initial checks make sure wrapper is given numerical values when it expects it

#currentupwindordered = "upwindorderedMATHEVALv4"

QPotential <- function (x.rhs = 'NULL', x.start = 'NULL', x.bound = 'NULL', x.num.steps = 'NULL', y.rhs = 'NULL', y.start = 'NULL', y.bound = 'NULL', y.num.steps = 'NULL', filename = 'NULL', save.to.R = 'NULL', save.to.HD = TRUE, bounce = 'd', bounce.edge = 0.01)
{

#currentupwindordered = "upwindorderedMATHEVALv4"
# Break apart function parameters into things that C code will use

# check if any component of x is missing
if (x.rhs == 'NULL') {stop("No equation defined for x. Define x.rhs.")} else {equationx <- x.rhs}
if (x.start == 'NULL') {stop("No starting value for x. Define x.start.")} else {startxval <- x.start}
if (x.num.steps == 'NULL') {stop("Need the number of steps in x range. Define x.num.steps.")} else {numofstepsx <- x.num.steps}
lengthequationx	<- nchar(x.rhs)
if (length(x.bound) < 2) stop("Not enough values for x range in parameter x.bound.")
if (length(x.bound) > 2) stop("Too many values for x range in parameter x.bound.")
lowerboundsx	<- x.bound[1]
upperboundsx	<- x.bound[2]


# check if any component of y is missing
if (y.rhs == 'NULL') {equationy = '0'} else {equationy <- y.rhs}
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
storage <- array(1.0, dim=c(1,(numofstepsx*numofstepsy)))

#Save in whatever format the user wants
if (filename != 'NULL') {lengthfilename = nchar(filename)}
else 					{lengthfilename = 0}

#default filename has a restriction on the name size
if ( (filename == 'NULL') && (save.to.HD == TRUE) && ( (abs(x.start) > 99999) || (abs(y.start) > 99999) ) ) { stop('Cannot use default filename because program will crash.  Please supply filename=')}

if ((save.to.R == TRUE) && (save.to.HD == TRUE)) 		{datasave = 3}
else if ((save.to.R == TRUE) && (save.to.HD == FALSE)) 	{datasave = 2}
else if (isTRUE(save.to.HD))							{datasave = 1}
else												{datasave = 4}
if (save.to.HD == 'testrun')							{datasave = 4}
print(paste("Variable datasave is: ", datasave, sep = ""))

# bounce.edge is defined in function
if (!is.numeric(bounce.edge)) stop('Parameter bounce.edge must be a number')

if (bounce == FALSE) {bouncestyle = 'd'}
else if (bounce == 'd') {bouncestyle = 'd'}
else if ((bounce == 'p') || (bounce == 'b')) {bouncestyle = bounce}
else 	{stop('Parameter bounce must be left blank, (d)efault), (p)ositivevalues, or (b)ounce')}

# ----------------------------------------------------------------------
# Error checking before C code is called
# ----------------------------------------------------------------------

# Check to make sure upwindorderedMATHEVAL.so exists
# if so, then load it
# if not, make it then load it
#try( dyn.load(paste(currentupwindordered, ".so", sep = "")) )
#system("R CMD SHLIB -I/usr/local/include -L/usr/local/lib -lmatheval upwindorderedMATHEVAL.c -lm")

{ # beginning of error checking
# check to make sure:
# check that LX1 < xeq < LX2
# check that LY1 < yeq < LY2
if ( (lowerboundsx > startxval) || (upperboundsx < startxval) ) {stop("Starting x value is outside x range")}
if ( (lowerboundsy > startyval) || (upperboundsy < startyval) ) {stop("Starting y value is outside y range")}
if ( lowerboundsx > upperboundsx ) {stop("Upper bounds of x is less than lower bounds of x")}
if ( lowerboundsy > upperboundsy ) {stop("Upper bounds of y is less than lower bounds of y")}
 


# warn about memory size is number of steps is very large
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
# Produce multiple versions: 	1) store output in R and on harddrive
#								2) store output only on harddrive
#								3) store output only in R

#TODO check that removing out2 <- should prevent R from printing everything out,
#probably not required because of call to storage 
#out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounce.edge))
if (datasave == 1) {
	#no R write; HD write
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounce.edge), PACKAGE="QPot")
	return(TRUE)
}
else if (datasave == 2) {
	#R write; no HD write
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounce.edge), PACKAGE="QPot")
	storage = out2[[1]]
	storage <- matrix(storage, nrow = x.num.steps, byrow = TRUE)
	return(storage)
}
else if (datasave == 3) {
	# R write; HD write
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounce.edge), PACKAGE="QPot")
	storage = out2[[1]]
	storage <- matrix(storage, nrow = x.num.steps, byrow = TRUE)
	return(storage)
}
else if (datasave == 4) {
	# no R write; no HD write
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounce.edge), PACKAGE="QPot")
	return(TRUE)
}
else 					{print("datasave is not a possible number.  How did you get here?")}

#out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounce.edge))
#storage = out2[[1]]
#return(storage)
} # end of upwind ordered function call
