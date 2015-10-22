# An R package for stochastic differential equation quasi-potential analysis

### Christopher M. Moore, Christopher R. Stieha, Ben C. Nolting, Maria K. Cameron, and Karen C. Abbott

QPot offers a range of tools to simulate, analyze, and visualize the dynamics of two-dimensional systems of stochastic differential equations.  QPot offers tools to compute the quasi-potential, which is useful when comparing the relative stabilities of equilibria in systems with multiple stable equilibria. 

Literature:
M. K. Cameron. 2012. Finding the quasipotential for nongradient SDEs. *Physica D*, 241(18):1532–1550.

Moore, C.M., Stieha, C.R., Nolting, B.C., Cameron, M.K., and Abbott, K.C. *Submitted to The R Journal.* QPot: An R package for stochastic differential equation quasi-potential analysis

B. C. Nolting and K. C. Abbott. *Accepted.* Balls, cups, and quasi-potentials: quantifying stability in stochastic systems. *Ecology.*

Uses the expression_parser library from http://jamesgregson.blogspot.com/2012/06/mathematical-expression-parser-in-c.html, licensed under GPL v2. Reads strings of mathematical expression in infix notation.  
* standard arithmetic operations (+,-,*,/) with operator precedence
* exponentiation ^ and nested exponentiation
* unary + and -
* expressions enclosed in parentheses ('(',')'), optionally nested
* built-in math functions: pow(x,y), sqrt(x), log(x), exp(x), sin(x), asin(x), cos(x), acos(x), tan(x), atan(x), atan2(y,x), abs(x), fabs(x), floor(x), ceil(x), round(x), with input arguments checked for domain validity, e.g. 'sqrt( -1.0 )' returns an error.
 
*This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation version 2 of the License.*
*The library expression_parser (https://github.com/jamesgregson/expression_parser) is used under the GPLv2 for non-commerical use.*

