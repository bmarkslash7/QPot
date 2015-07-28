An R package for stochastic differential equation quasi-potential analysis

Christopher M. Moore, Christopher Stieha, Ben C. Nolting, Maria K. Cameron, and Karen C. Abbott

Literature:


Uses the matheval library from http://www.gnu.org/software/libmatheval/, licensed under GPL v3. 

```
File Structure:
all of these are not needed;  if empty, remove
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
```



