#' Test run of QPotential function
#' 
#' 
#' @return filetoHD

defaultTest <- function () {
#Test Case 1 from the QPot dropbox folder
xbounds = c(-0.5, 20.0)
ybounds = c(-0.5, 20.0)
xstepnumber = 4100
ystepnumber = 4100
#xstepnumber = 1000
#ystepnumber = 1000
xinit = 1.40491
yinit = 2.80808
testequationx = "1.54*x*(1.0-(x/10.14))-(y*x*x)/(1.0+x*x)"
testequationy = "((0.476*x*x*y)/(1+x*x))-0.112590*y*y"

choice <- readline(prompt = "This function runs a series of tests on the function QPotential() that takes a lot of time, hard drive space, and RAM.  This could cause your computer to be unhappy.  Do you want to continue (y/N)?")
if ( (choice == "Y") || (choice == "y") ) {print("I warned you.\n")}
else (stop("Exiting defaultTest function."))


print("########################################################")
print("Starting QPotential()")
print("########################################################")
print("DEFAULT: This writes the file defaultname-x1.4049y2.8081.txt to your hard drive")
QPotential(x.rhs = testequationx, x.start = xinit, x.bound = xbounds, x.num.steps = xstepnumber, 
			y.rhs = testequationy, y.start = yinit, y.bound = ybounds, y.num.steps = ystepnumber, 
			filename = 'NULL', save.to.R = 'NULL', save.to.HD = TRUE, 
			bounce = 'd', bounce.edge = 0.01, 
			verboseR = FALSE, verboseC = FALSE, debugC = FALSE)
			
#print("DEFAULT: This writes the file default file with filename CHRISdefaulttest.txt to your hard drive")
QPotential(x.rhs = testequationx, x.start = xinit, x.bound = xbounds, x.num.steps = xstepnumber, 
			y.rhs = testequationy, y.start = yinit, y.bound = ybounds, y.num.steps = ystepnumber, 
			filename = 'CHRISdefaulttest.txt', save.to.R = 'NULL', save.to.HD = TRUE, 
			bounce = 'd', bounce.edge = 0.01, 
			verboseR = FALSE, verboseC = FALSE, debugC = FALSE)

print("writeHDwriteR: This writes the file defaultname-writeHDwriteR.txt to your hard drive and stores the file in your instance of R")
storage <- 
QPotential(x.rhs = testequationx, x.start = xinit, x.bound = xbounds, x.num.steps = xstepnumber, 
			y.rhs = testequationy, y.start = yinit, y.bound = ybounds, y.num.steps = ystepnumber, 
			filename = 'defaultname-writeHDwriteR.txt', save.to.R = TRUE, save.to.HD = TRUE, 
			bounce = 'd', bounce.edge = 0.01, 
			verboseR = FALSE, verboseC = FALSE, debugC = FALSE)
print("Size of writeHDwriteR storage (R matrix):")
print("(May be backwards)")
print(paste("Number of Rows ", nrow(storage), " should be ", xstepnumber, sep=""))
print(paste("Number of Columns ", ncol(storage), " should be ", ystepnumber, sep=""))

print("nowriteHDwriteR: This stores the results of QPotential() in your instance of R only")
storage2 <- 
QPotential(x.rhs = testequationx, x.start = xinit, x.bound = xbounds, x.num.steps = xstepnumber, 
			y.rhs = testequationy, y.start = yinit, y.bound = ybounds, y.num.steps = ystepnumber, 
			save.to.R = TRUE, save.to.HD = FALSE, 
			bounce = 'd', bounce.edge = 0.01, 
			verboseR = FALSE, verboseC = FALSE, debugC = FALSE)
print("Size of nowriteHDwriteR storage (R matrix):")
print("(May be backwards)")
print(paste("Number of Rows ", nrow(storage2), " should be ", xstepnumber, sep=""))
print(paste("Number of Columns ", ncol(storage2), " should be ", ystepnumber, sep=""))


print("########################################################")
print("Reading in matrices on hard drive")
print("########################################################")
TEMP_CHRISdefault <- read.table(file = "CHRISdefaulttest.txt", sep = "\t", header = FALSE) # I ran this after running writeHDwriteR
TEMP_default <- read.table(file = "defaultname-x1.4049y2.8081.txt", sep = "\t", header = FALSE)
TEMP_HD_writeHDwriteR <- read.table(file = "defaultname-writeHDwriteR.txt", sep = "\t", header = FALSE) # this is when writeHDwriteR is run first
#TEMP_HD_writeHDwriteRSECOND <- read.table(file = "TESTdefaultname-writeHDwriteR.txt", sep = "\t", header = FALSE) #this is when writeHDwriteR is run after default

TEMP_CORRECT <- read.table(file = "CORRECTdefaultname-x1.4049y2.8081.txt", sep = "\t", header = FALSE)

print("########################################################")
print("Plotting the Matrices for easy comparison")
print("########################################################")
QPContour(as.matrix(TEMP_CORRECT), c(4100,4100), xbounds, ybounds, c.parm=5)
QPContour(as.matrix(TEMP_CHRISdefault), c(4100,4100), xbounds, ybounds, c.parm=5)
QPContour(as.matrix(TEMP_div2), c(4100,4100), xbounds, ybounds, c.parm=5)
QPContour(as.matrix(TEMP_HD_withHDwriteR), c(4100,4100), xbounds, ybounds, c.parm=5)
QPContour(storage, c(4100,4100), xbounds, ybounds, c.parm=5)

print("########################################################")
print("Testing for equality among matrices")
print("########################################################")
########################################################################
# USE as.data.frame before using as.matrix
########################################################################
if (isTRUE(all.equal(TEMP_default, TEMP_CORRECT, tolerance = 10^-4))) {print("SUCCESS - Code from DEFAULT and TEMP_CORRECT get same answer")} else {print("FAIL - code from DEFAULT and TEMP_CORRECT produced different matrices")} #PASSES

if (isTRUE(all.equal(TEMP_HD_writeHDwriteR, TEMP_CORRECT, tolerance = 10^-4))) {print("SUCCESS - Code from TEMP_HD_writeHDwriteR and TEMP_CORRECT get same answer")} else {print("FAIL - code from TEMP_HD_writeHDwriteR and TEMP_CORRECT produced different matrices")} #FAILS

if (isTRUE(all.equal(TEMP_default, TEMP_HD_writeHDwriteR, tolerance = 10^-4))) {print("SUCCESS - Code from DEFAULT and writeHDwriteR get same answer")} else {print("FAIL - code from DEFAULT and writeHDwriteR produced different matrices")} #FAILS

if (isTRUE(all.equal(TEMP_default, storage, tolerance = 10^-4))) {print("SUCCESS - Code from DEFAULT and storage-writeHDwriteR get same answer")} else {print("FAIL - code from DEFAULT and storage-writeHDwriteR produced different matrices")}

if (isTRUE(all.equal(TEMP_CHRISdefault, TEMP_CORRECT, tolerance = 10^-4))) {print("SUCCESS - Code from CHRIS and CORRECT get same answer")} else {print("FAIL - code from CHRIS and CORRECT produced different matrices")} #FAILS

if (isTRUE(all.equal(TEMP_HD_writeHDwriteR, storage, tolerance = 10^-4))) {print("SUCCESS - Code from TEMP_HD_writeHDwriteR and storage-writeHDwriteR get same answer")} else {print("FAIL - code from TEMP_HD_writeHDwriteR and storage-writeHDwriteR produced different matrices")} #FAIL



}
