# An R package for stochastic differential equation quasi-potential analysis

### Christopher M. Moore, Christopher Stieha, Ben C. Nolting, Maria K. Cameron, and Karen C. Abbott

Literature:


Uses the expression_parser library from http://jamesgregson.blogspot.com/2012/06/mathematical-expression-parser-in-c.html, licensed under GPL v2. Reads strings of mathematical expression in infix notation.  
*standard arithmetic operations (+,-,*,/) with operator precedence
*exponentiation ^ and nested exponentiation
*unary + and -
*expressions enclosed in parentheses ('(',')'), optionally nested
*built-in math functions: pow(x,y), sqrt(x), log(x), exp(x), sin(x), asin(x), cos(x), acos(x), tan(x), atan(x), atan2(y,x), abs(x), fabs(x), floor(x), ceil(x), round(x), with input arguments checked for domain validity, e.g. 'sqrt( -1.0 )' returns an error.
 
File Structure:
```
all of these are not needed;  if empty, remove---or R will do it for you.
https://cran.r-project.org/doc/manuals/R-exts.html#Package-structure
./configure	configure instructions for win, linux, mac
./data		data files
./demo		R scripts that show functionality of package, so examples from journal article(?);
		needs "00Index file with one line for each demo, giving its name and a description separated by a tab or at least three spaces."
./exec		contains executable scripts , usually for Perl, shell, etc.
./inst		"The contents of the inst subdirectory will be copied recursively to the installation directory." 
		include CITATION file if automatic file from DESCRIPTION is not good enough
./man		R documentation in .Rd format; QPot-internal.Rd defines internal objects unavailable to user
./po		used for localization
./R		source code
./src		source and headers for compiled code, Makefile if needed
./test		package-specific test code
./tools		"place for auxiliary files needed during configuration, and also for sources need to re-create scripts"
./vignettes
DESCRIPTION	Basic information on authors, licenses, etc.
NAMESPACE	lists what is imported (from other libraries) and exported (made available to the user)
NEWS		file containing updates and news
README.md	file containing basic information on the package
TODO		list of things to add; can also be incorporated into the issue tracker.
loadall.R	R script to load all files not as a package
READMEN-chris.txt		stieha todo file with notes - delete by 1 sept 2015
```



