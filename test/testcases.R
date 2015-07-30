# R CMD SHLIB upwindorderedMATHEVALv3.c DOES NOT WORK
# R CMD SHLIB upwindorderedMATHEVALv3.c -lm  DOES NOT WORK
# R CMD SHLIB -lmatheval upwindorderedMATHEVALv3.c -lm WORKS, but why?
# R CMD SHLIB -L./libmatheval -lmatheval upwindorderedMATHEVALv3.c -lm WORKS, but is like above, which works

testQPotential <- function() {
#dyn.load(paste(currentupwindordered, ".so", sep=""))
xbounds = c(-5.0, 60.0)
ybounds = c(-5.0, 60.0)
xstepnumber = 1000
ystepnumber = 1000
xinit = 6.60341
yinit = 3.04537
testequationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
testequationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
equationx = testequationx
equationy = testequationy

storage <- array(1.0, dim=c(1,(xstepnumber*ystepnumber)))

#full code
QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber)
print("fullcode completed")
#missing a function parameter
print("no x equation")
try(QPotential(xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=equationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber))
print("no initial x")
try(QPotential(xrhs=equationx,xrange=xbounds, xsteps=xstepnumber, yrhs=equationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber))
print("no x bounds")
try(QPotential(xrhs=equationx, xstart=xinit,xsteps=xstepnumber, yrhs=equationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber))
print("only one x bounds")
try(QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds[1], xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber))
print("more than two x bounds")
try(QPotential(xrhs=testequationx, xstart=xinit, xrange=c(xbounds, 4), xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber))
print("no x step number")
try(QPotential(xrhs=equationx, xstart=xinit, xrange=xbounds, yrhs=equationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber))
print("no y equation")
try(QPotential(xrhs=equationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, ystart=yinit, yrange=ybounds, ysteps=ystepnumber))
print("no initial y")
try(QPotential(xrhs=equationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=equationy, yrange=ybounds, ysteps=ystepnumber))
print("no y bounds")
try(QPotential(xrhs=equationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=equationy, ystart=yinit, ysteps=ystepnumber))
print("no y step number")
try(QPotential(xrhs=equationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=equationy, ystart=yinit, yrange=ybounds))

print("test suite finished")

} #end of QPotential test suite

testdatawrite <- function() {
#dyn.load(paste(currentupwindordered, ".so", sep=""))
xbounds = c(-5.0, 60.0)
ybounds = c(-5.0, 60.0)
xstepnumber = 1000
ystepnumber = 1000
xinit = 6.60341
yinit = 3.04537
testequationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
testequationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
equationx = testequationx
equationy = testequationy

storage <- array(1.0, dim=c(1,(xstepnumber*ystepnumber)))

print("starting to test data writing")
QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetest.txt', savetoR = 'NULL', savetoHD = TRUE, bounce = 'd', bounceedge = 0.01)
Vdefault <-read.table("datawritetest.txt", header=FALSE, sep="\t")

print("Bounce test: Use default.  Should equal datawritetest.txt")
QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestbounce-default.txt', savetoR = 'NULL', savetoHD = TRUE)
V01<-read.table("datawritetestbounce-default.txt", header=FALSE, sep="\t")
if( isTRUE(all.equal(V01, Vdefault)) ) {print("RESULT: datawritetestbounce-default.txt is equal to datawritetest.txt")}
else {print("RESULT: datawritetestbounce-default.txt DOES NOT EQUAL datawritetest.txt")}

print("Bounce test: Use bounce.  Should NOT equal datawritetest.txt")
QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestbounce-bounce.txt', savetoR = 'NULL', savetoHD = TRUE, bounce = 'b')
V01<-read.table("datawritetestbounce-bounce.txt", header=FALSE, sep="\t")
if( isTRUE(all.equal(V01, Vdefault)) ) {print("RESULT: datawritetestbounce-bounce.txt is equal to datawritetest.txt")}
else {print("RESULT: datawritetestbounce-bounce.txt DOES NOT EQUAL datawritetest.txt")}

print("Bounce test: Use positive.  Should NOT equal datawritetest.txt")
QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestbounce-positive.txt', savetoR = 'NULL', savetoHD = TRUE, bounce = 'b')
V01<-read.table("datawritetestbounce-positive.txt", header=FALSE, sep="\t")
if( isTRUE(all.equal(V01, Vdefault)) ) {print("RESULT: datawritetestbounce-positive.txt is equal to datawritetest.txt")}
else {print("RESULT: datawritetestbounce-positive.txt DOES NOT EQUAL datawritetest.txt")}


}

testHDwriteRwrite <- function() {
#TESTED THESE BY HAND ON 23 JULY 2015
#currentupwindordered = "upwindorderedMATHEVALv4"
#dyn.load(paste(currentupwindordered, ".so", sep=""))
xbounds = c(-5.0, 60.0)
ybounds = c(-5.0, 60.0)
xstepnumber = 1000
ystepnumber = 1000
xinit = 6.60341
yinit = 3.04537
testequationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
testequationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
equationx = testequationx
equationy = testequationy

FAIL = 0

storage <- array(1.0, dim=c(1,(xstepnumber*ystepnumber)))
onlyones <- array(1.0, dim=c(1,(xstepnumber*ystepnumber)))

print("starting to test HD writing and R writing")
print("write HD, no write R, default listed explicitly")
#correctly writes to HD, correctly does not write to R  -  24 July 2015 CRS
tempdefault <- QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestALLDEFAULTSlisted.txt', savetoR = 'NULL', savetoHD = TRUE)
#try(Vdefault <-read.table("datawritetestONLYHD.txt", header=FALSE, sep="\t"))

print("write HD, no write R, defaults not listed") 
#correctly writes to HD, correctly does not write to R -  24 July 2015 CRS
temp2 <- QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestALLDEFAULTSnotlisted.txt')

V01<-read.table('datawritetestALLDEFAULTSlisted.txt', header=FALSE, sep="\t")
V02<-read.table('datawritetestALLDEFAULTSnotlisted.txt', header=FALSE, sep="\t")
if (isTRUE(all.equal(V01, V02))) {print("PASSED: HD write: defaults (listed) and (not listed) the same")}
else {print("FAILED: HD write: defaults not the same"); FAIL = FAIL + 1}

if (isTRUE(tempdefault) && isTRUE(temp2) && isTRUE(all.equal(tempdefault, temp2))) {print("PASSED: R write: defaults listed and not listed the same")}
else {print("FAILED: R write: defaults not the same"); FAIL = FAIL + 1}
#rm(temp1)
#rm(temp2)
#gc()

print("********************************************")
print("write to HD, write to R")
rm(temp)
tempfull <- QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestBOTHRandHD.txt', savetoR = TRUE, savetoHD = TRUE)
V01<-read.table('datawritetestALLDEFAULTSlisted.txt', header=FALSE, sep="\t")
V02<-read.table('datawritetestBOTHRandHD.txt', header=FALSE, sep="\t")
write.table(tempfull, file = "datawritetest-datafromtempfull-writeR.txt", row.names=FALSE, col.names=FALSE)
V03<-read.table("datawritetest-datafromtempfull-writeR.txt", header=FALSE)
if (isTRUE(all.equal(V01, V02, tolerance = 1e-4))) {print("PASSED: HD write: defaults and wHDwR the same")}
else {print("FAILED: HD write: defaults and wHDwR not the same"); FAIL = FAIL + 1}
if (isTRUE(all.equal(V01, V03, tolerance = 1e-4))) {print("PASSED: HD write: defaults and wR the same")}
else {print("FAILED: HD write: defaults and wR not the same"); FAIL = FAIL + 1}
if (isTRUE(all.equal(V02, V03, tolerance = 1e-4))) {print("PASSED: HD write: wHDwR and wR the same")}
else {print("FAILED: HD write: wHDwR and wR not the same"); FAIL = FAIL + 1}

if (isTRUE(all.equal(tempfull, tempdefault, tolerance = 1e-4))) {print("FAILED: R write: defaults and wHDwR the same and should not be"); FAIL = FAIL + 1}
else {print("PASSED: HD write: defaults and wHDwR different")}
if (isTRUE(all.equal(tempfull, tempdefault, tolerance = 1e-4))) {print("FAILED: R write: temp2 and wHDwR the same and should not be"); FAIL = FAIL + 1}
else {print("PASSED: HD write: temp2 and wHDwR different")}


print("********************************************")
print("no write HD, no write R NULL")
rm(temp)
temp <- QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestonlyRNULL-SHOULDNTWRITE.txt', savetoR = 'NULL', savetoHD = FALSE)
V01<-read.table('datawritetestALLDEFAULTSlisted.txt', header=FALSE, sep="\t")
rm(V02)
V02='NULL'
try(V02<-read.table('datawritetestonlyRNULL-SHOULDNTWRITE.txt', header=FALSE, sep="\t"))
if (isTRUE(all.equal(V01, V02, tolerance = 1e-4))) {print("FAILED: no HD no R: default and readHD the same and should not be"); FAIL = FAIL + 1}
else {print("PASSED: HD write: default and readHD different")}


#TODO: temp should be empty. temp is full of ones
print("********************************************")
print("write to HD, no write to R FALSE")
rm(temp)
temp <- QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestonlyHD.txt', savetoR = FALSE, savetoHD = TRUE)
V01<-read.table('datawritetestALLDEFAULTSlisted.txt', header=FALSE, sep="\t")
try(V02<-read.table('datawritetestonlyHD.txt', header=FALSE, sep="\t"))
if (isTRUE(all.equal(V01, V02, tolerance = 1e-4))) {print("PASSED: HD write: defaults and write HD, no R write are the same")}
else {print("FAILED: HD write: defaults and write HD, no R write not the same"); FAIL = FAIL + 1}

write.table(temp, file = "data-temp.txt", row.names=FALSE, col.names=FALSE)
V03<-read.table("data-temp.txt", header=FALSE)
if (isTRUE(all.equal(V02, V03, tolerance = 1e-4))) {print("FAILED: HD write, no R write: HD write  and R write the same and should not be"); FAIL = FAIL + 1}
else {print("PASSED: HD write, no R write: HD write and R write different")}
if (isTRUE(all.equal(V01, V03, tolerance = 1e-4))) {print("FAILED: HD write, no R write: default and R write the same and should not be"); FAIL = FAIL + 1}
else {print("PASSED: HD write, no R write: default and R write different")}


## STOPPED HERE WITH TEST CASES
if (FAIL == 0) {print("ALL TESTS PASSED")}
if (FAIL > 0) {print(paste("FAILURES OCCURRED: ",FAIL, " fails", sep=""))}
return(TRUE)
#############################################################3
#############################################################
#################################################################

print("********************************************")
print("no write to HD, no write to R FALSE")
QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestonlyRFALSE-SHOULDNTWRITE.txt', savetoR = FALSE, savetoHD = FALSE)
V01<-read.table('datawritetestALLDEFAULTSlisted.txt', header=FALSE, sep="\t")
try(V02<-read.table('datawritetestonlyRFALSE-SHOULDNTWRITE.txt', header=FALSE, sep="\t"))
try(all.equal(V01, V02))

#This shouldn't write anything anywhere, R matrix is all 1s
#this doesn't write anything to the HD, yeah!
print("********************************************")
print("HD is testrun, no write to R FALSE")
rm(temp)
temp <- QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber, filename = 'datawritetestonlytestrun.txt', savetoR = FALSE, savetoHD = 'testrun')

if (FAIL > 0) {print(paste("FAILURES OCCURRED: ",FAIL, "fails", sep=""))}

}
