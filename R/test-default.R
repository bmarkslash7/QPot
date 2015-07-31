#' Test run of QPotential function
#' 
#' @export
#' 
#' @param NONE
#' @return filetoHD

defaultTest <- function () {
xbounds = c(-5.0, 60.0)
ybounds = c(-5.0, 60.0)
xstepnumber = 1000
ystepnumber = 1000
xinit = 6.60341
yinit = 3.04537
testequationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
testequationy = "-4.0*y+((10.0*x*y)/(18.0+x))"

#full code
print("This writes a file to your hard drive\n")
QPotential(xrhs=testequationx, xstart=xinit, xrange=xbounds, xsteps=xstepnumber, yrhs=testequationy, ystart=yinit, yrange=ybounds, ysteps=ystepnumber)


}
