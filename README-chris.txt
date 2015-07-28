notes for R wrapper:


notes for c code:
Figure out how to compile using "./libmatheval/matheval.h"


Files:
loadall.R					recompiles and loads C code, 
							loads R wrapper, loads testcases
upwindordered.R				R wrapper to call Maria's C code
upwindorderedMATHEVALvX		Maria's C code with Stieha additions
							X = version number
testcases.R					A series of calls to upwindorered code
							to test various scenarios
							failure is checked by hand








#################################################
delete after 1 August 2015
#################################################

R CMD SHLIB -lmatheval upwindorderedMATHEVALv3.c -lm

dyn.load("upwindorderedMATHEVALv3.so")
teststepsize =  1000
storage <- array(1.0, dim=c(1,(teststepsize*teststepsize)))
equationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
lengthequationx <- nchar(equationx)
equationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
lengthequationy <- nchar(equationy)
out2 <- .C("quasipotentialSTRINGSwithstoragePOINTERSreturnvalue", as.double(storage), as.double(-5.0), as.double(60.0), as.integer(teststepsize), as.double(-5.0), as.double(60.0), as.integer(teststepsize), as.double(6.60341), as.double(3.04537), equationx, as.integer(lengthequationx), equationy, as.integer(lengthequationy))
