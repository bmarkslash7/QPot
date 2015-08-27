#' Vector decomposition and remainder fields
#'
#' This function calculates the remainder field.
#' @param surface matrix output from QPGlobal.
#' @param x.rhs COPY FROM QPOTENTIAL().
#' @param y.rhs COPY FROM QPOTENTIAL().
#' @param x.bound two-element vector with respective minimum and maximum x values.
#' @param y.bound two-element vector with respective minimum and maximum y values.
#' @keywords vector field decompoosition, remainder field
#' 
#' @examples
#' #Example 1 from article
#' equationx = "1.5*x*(1.0-(x/45.0))-(y*x*5.0)/(18.0+x)"
#' equationy = "-4.0*y+((10.0*x*y)/(18.0+x))"
#' VDR <- VecDecompRem(surface=e1.global, x.rhs=testequationx, y.rhs=testequationy, x.bound=c(-0.5,20), y.bound=c(-0.5,20))
#' VecDecompPlot(field=list(VDR[,,1],VDR[,,2]), dens=c(50,50), x.bound=c(-0.5,20), y.bound=c(-0.5,20))

VecDecompRem <- function(surface,x.rhs,y.rhs,x.bound,y.bound){
	qpr <- nrow(surface)
	qpc <- ncol(surface)
	equations <- list(x.rhs,y.rhs)

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

	#remainder fields
	x.val <- seq(min(x.bound),max(x.bound),length.out=qpr)
	y.val <- seq(min(y.bound),max(y.bound),length.out=qpc)
	n.eq <- length(equations)
	z.list <- vector(mode="list" , length = n.eq)
	for(i in 1:n.eq){
		z.list[[i]] <- 
		outer(x.val,y.val,function(x,y){eval(parse(text=equations[i]))})
		}


	#remainder field
	rx <- z.list[[1]]+dr
	ry <- z.list[[2]]+dc

	array(data=c(rx,ry),dim=c(qpr,qpc,2))
}
