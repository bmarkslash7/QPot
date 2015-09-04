#' Finding the global quasi-potential
#'
#' This function allows you to find the global quasi-potential values for several local quasi-potential surfaces
#' @param local.surfaces A list of local quasi-potential surfaces, each of which is stored in discretized form as a matrix.
#' @param unstable.eq.x A vector of the x-coordinates of the unstable equilibria.  Must be in the same order as unstable.eq.y.
#' @param unstable.eq.y A vector of the y-coordinates of the unstable equilibria.  Must be in the same order as unstable.eq.x.
#' @param x.bound A two-element vector with the minimum and maximum x values used for computing the quasi-potential.
#' @param y.bound A two-element vector with the minimum and maximum y values used for computing the quasi-potential.
#' @keywords Global quasi-potential
#' 
########################################################################
# COMMENTED OUT FOR PASSING devtools::check() - DECLARE local.1 and local.2
########################################################################
# #' @examples 
# #' QPGlobal(list(local.1,local.2),c(0,4),c(0,4),c(-1,5),c(-1,5))

QPGlobal <- function(local.surfaces , unstable.eq.x , unstable.eq.y , x.bound , y.bound) {
	n.surfaces <- length(local.surfaces)
	n.unstable.eq.x <- length(unstable.eq.x)
	n.unstable.eq.y <- length(unstable.eq.y)
		if(n.unstable.eq.x != n.unstable.eq.y){stop("Unstable x and y points not equal")}
	n.unstable.pts <- length(unstable.eq.x)
	mesh.xy <- dim(local.surfaces[[1]])
	x.range <- max(x.bound)-min(x.bound)
	y.range <- max(y.bound)-min(y.bound)

	unstable.xy <- cbind(unstable.eq.x , unstable.eq.y) #unstable eq. as one object
	unstable.phi.loc <- matrix(data = NA , nrow = n.unstable.pts , ncol = 2 , byrow = T, dimnames=list(paste("unstab.eq",1:n.unstable.pts,sep="") , c("x","y"))) #local indeces

	for (i in 1:n.unstable.pts){ #indexes unstable cooridnates
		x.loc <- round((unstable.xy[i,][1]-min(x.bound))/x.range*mesh.xy[1])
		y.loc <- round((unstable.xy[i,][2]-min(y.bound))/y.range*mesh.xy[2])
		unstable.phi.loc[i,] <- c(x.loc, y.loc)
		}

	unstable.phi <- matrix(data = NA , nrow = n.unstable.pts , ncol = n.surfaces , byrow = F , dimnames=list(paste("unstab.eq",1:n.unstable.pts,sep=""),paste("surface",1:n.surfaces,sep=""))) #local phi values for each unstable equilibrium pair
	for(i in 1:n.unstable.pts) {
		for (j in 1:n.surfaces) {
		unstable.phi[i,j] <- local.surfaces[[j]][unstable.phi.loc[i,1],unstable.phi.loc[i,2]]
		}
	}

	max.phi <- max(unstable.phi,na.rm=T)
	max.phi.arr <- which(unstable.phi == max.phi , arr.ind=T)

	if(nrow(max.phi.arr) != n.unstable.pts){ #if not all max(phi) are the same, then they need to be aligned
		global.max <- unstable.phi[max.phi.arr[1,1],max.phi.arr[1,2]]
		phi.diff <- matrix(data = NA , nrow = n.surfaces , ncol = 1, dimnames=list(paste("surface",1:n.surfaces,sep="")))
		for (i in 1:n.surfaces){
			phi.diff[i,] <- global.max - unstable.phi[max.phi.arr[1,1],i]
		}
	} else {
		phi.diff <- matrix(data = 0 , nrow = n.surfaces , ncol = 1)
	}

	eq.arr.unadj <- array(data = unlist(local.surfaces) , dim = c(mesh.xy[1] , mesh.xy[2] , n.surfaces))
	eq.arr <- array(data = NA , dim = c(mesh.xy[1] , mesh.xy[2] , n.surfaces))
	for(i in 1:n.surfaces){
		eq.arr[,,i] <- eq.arr.unadj[,,i]+phi.diff[i]
	}
	rm(eq.arr.unadj)
	gc()
	global.qp <- apply(eq.arr , c(1:n.surfaces) , min)
	rm(eq.arr)
	gc()
	global.qp
	}