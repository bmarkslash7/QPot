#loadall.R reads in the three function files
print("Recompiling C code")
#system("R CMD SHLIB -I/usr/local/include -L/usr/local/lib -lmatheval upwindorderedMATHEVALv4.c -lm")
system("R CMD SHLIB ./src/upwindorderedMATHEVALv4.c -lm")	#compiles, but does not find evaluator_evaluate_x_y
system("R CMD SHLIB -Lsrc/libmatheval -lmatheval upwindorderedMATHEVALv4.c -lm") # no rule to make .o needed for .so
system("R CMD SHLIB -Isrc/libmatheval -lmatheval upwindorderedMATHEVALv4.c -lm") #no rule to make .o needed by .so
system("R CMD SHLIB -o src/libmatheval/matheval.o -c src/libmatheval/matheval.c")
system("R CMD SHLIB src/upwindorderedMATHEVALv4.c src/libmatheval/matheval.o -lm") # no rule to make upwind.o needed by upwind.so
print("Reloading C code into R")
currentupwindordered = "src/upwindorderedMATHEVALv4"
try( dyn.load(paste(currentupwindordered, ".so", sep = "")) )
print("Loading functions")
source('./R/upwindordered.R')
source('./test/testcases.R')


gcc -o myprog myprog.c -L/home/newhall/lib -lmine
