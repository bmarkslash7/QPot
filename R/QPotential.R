#' Wrapper for call to quasipotential compution using the upwind ordered method
#' 
#' @export
#' 
#' @param xrhs xstart xrange xsteps yrhs ystart yrange ysteps filename savetoR savetoHD bounce bounceedge MORE
#' @return filetoHD filetoR


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

QPotential <- function (xrhs = 'NULL', xstart = 'NULL', xrange = 'NULL', xsteps = 'NULL', yrhs = 'NULL', ystart = 'NULL', yrange = 'NULL', ysteps = 'NULL', filename = 'NULL', savetoR = 'NULL', savetoHD = TRUE, bounce = 'd', bounceedge = 0.01)
{

#currentupwindordered = "upwindorderedMATHEVALv4"
# Break apart function parameters into things that C code will use

# check if any component of x is missing
if (xrhs == 'NULL') {stop("No equation defined for x. Define xrhs.")} else {equationx <- xrhs}
if (xstart == 'NULL') {stop("No starting value for x. Define xstart.")} else {startxval <- xstart}
if (xsteps == 'NULL') {stop("Need the number of steps in x range. Define xsteps.")} else {numofstepsx <- xsteps}
lengthequationx	<- nchar(xrhs)
if (length(xrange) < 2) stop("Not enough values for x range in parameter xrange.")
if (length(xrange) > 2) stop("Too many values for x range in parameter xrange.")
lowerboundsx	<- xrange[1]
upperboundsx	<- xrange[2]


# check if any component of y is missing
if (yrhs == 'NULL') {equationy = '0'} else {equationy <- yrhs}
lengthequationy <- nchar(yrhs)
if ( (length(yrange) == 1) && (yrange == 'NULL') ) {lowerboundsy <- 0; upperboundsy <- 0}
if (yrange[1] != 'NULL') {
	if (length(yrange) < 2) stop("Not enough values for y range in variable yrange.")
	if (length(yrange) > 2) stop("Too many values for y range in variable yrange.")
	lowerboundsy <- yrange[1]; upperboundsy	<- yrange[2]
} # end of if yrange != 'NULL 
if (yrange[1] == 'NULL') stop('No minimum and maximum y values.  Parameter yrange not defined')
#The numofstepsy cannot equal 1, because hy=(LY2-LY1)/(NY-1), where NY is numofsteps
if (ysteps == 'NULL') {numofstepsy <- 2} else {numofstepsy <- ysteps}
if (ystart == 'NULL') {startyval = 0} else {startyval <- ystart}

#TODO only create this storage array if savetoR == TRUE
# I think the .C function creates its own storage matrix, despite what 
# I tell it or what it is supposed to do.
storage <- array(1.0, dim=c(1,(numofstepsx*numofstepsy)))

#Save in whatever format the user wants
if (filename != 'NULL') {lengthfilename = nchar(filename)}
else 					{lengthfilename = 0}

if ((savetoR == TRUE) && (savetoHD == TRUE)) 		{datasave = 3}
else if ((savetoR == TRUE) && (savetoHD == FALSE)) 	{datasave = 2}
else if (isTRUE(savetoHD))							{datasave = 1}
else												{datasave = 4}
if (savetoHD == 'testrun')							{datasave = 4}
print(paste("Variable datasave is: ", datasave, sep = ""))

# bounceedge is defined in function
if (!is.numeric(bounceedge)) stop('Parameter bounceedge must be a number')

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
#out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounceedge))
if (datasave == 1) {
	#no R write; HD write
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounceedge))
	return(TRUE)
}
else if (datasave == 2) {
	#R write; no HD write
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounceedge))
	storage = out2[[1]]
	storage <- matrix(storage, nrow = xsteps, byrow = TRUE)
	return(storage)
}
else if (datasave == 3) {
	# R write; HD write
	out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounceedge))
	storage = out2[[1]]
	storage <- matrix(storage, nrow = xsteps, byrow = TRUE)
	return(storage)
}
else if (datasave == 4) {
	# no R write; no HD write
	.C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounceedge))
	return(TRUE)
}
else 					{print("datasave is not a possible number.  How did you get here?")}

#out2 <- .C("quasipotential", as.double(storage), as.double(lowerboundsx), as.double(upperboundsx), as.integer(numofstepsx), as.double(lowerboundsy), as.double(upperboundsy), as.integer(numofstepsy), as.double(startxval), as.double(startyval), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy), filename, as.integer(lengthfilename), as.integer(datasave), bouncestyle, as.double(bounceedge))
#storage = out2[[1]]
#return(storage)
} # end of upwind ordered function call
