#' Vector decomposition and remainder fields
#'
#' This function calculates the gradient field.
#' @param surface matrix output from QPGlobal.
#' @keywords vector field decompoosition, gradient field
#' 
#' @examples
#' #Example 1 from article. Equations are not required
#' equationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
#' equationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
########################################################################
# Need to decide how to handle e1.global in examples
########################################################################
# #' e1.global <- matrix()
# #' VDG <- VecDecompGrad(e1.global)
# #' VecDecompPlot(field=list(VDG[,,1],VDG[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))


VecDecompGrad <- function(surface){
	qpr <- nrow(surface)
	qpc <- ncol(surface)

	# column derivative
	r.dc <- (surface[,qpc] - surface[,(qpc-1)])
	l.dc <- (surface[,2] - surface[,1])
	int.dc <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in (qpc-1):2){
		int.dc[,i] <- (surface[,(i+1)]-surface[,(i-1)])/2
	}
	int.dc[,1] <- l.dc
	int.dc[,qpc] <- r.dc
	dc <- -int.dc

	# row derivative
	t.dr <- (surface[1,] - surface[2,])
	b.dr <- (surface[(qpr-1),] - surface[qpr,])
	int.dr <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in 2:(qpr-1)){
		int.dr[i,] <- (surface[(i-1),]-surface[(i+1),])/2
	}
	int.dr[1,] <- t.dr
	int.dr[qpr,] <- b.dr
	dr <- int.dr

	array(data=c(dr,dc),dim=c(qpr,qpc,2))
}
