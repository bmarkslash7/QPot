#' Finding the global quasi-potential
#'
#' This function allows you to find the global quasi-potential values for several local quasi-potential surfaces
#' @param local.surfaces A list of local quasi-potential surfaces as matrices.
#' @param unstable.eq.x A vector of unstable equilibria.  Must be in the same order as unstable.eq.y.
#' @param unstable.eq.y A vector of unstable equilibria.  Must be in the same order as unstable.eq.x.
#' @param xlim A two-element vector with the minimum and maximum x values used for computing the quasi-potential.
#' @param ylim A two-element vector with the minimum and maximum y values used for computing the quasi-potential.
#' @keywords Global quasi-potential
#' @export
#' @examples 
#' QPGlobal(list(local.1,local.2),c(0,4),c(0,4),c(-1,5),c(-1,5))

QPGlobal <- function(local.surfaces , unstable.eq.x , unstable.eq.y , xlim , ylim) {
	n.surfaces <- length(local.surfaces)
	n.unstable.eq.x <- length(unstable.eq.x)
	n.unstable.eq.y <- length(unstable.eq.y)
		if(n.unstable.eq.x != n.unstable.eq.y){stop("Unstable x and y points not equal")}
	mesh.xy <- dim(local.surfaces[[1]])
	x.range <- max(xlim)-min(xlim)
	y.range <- max(ylim)-min(ylim)

	unstable.xy <- cbind(unstable.eq.x , unstable.eq.y) #unstable eq. as one object
	unstable.phi.loc <- matrix(data = NA , nrow = n.unstable.eq.x , ncol = 2 , byrow = T) #local indeces
	unstable.phi <- matrix(data = NA , nrow = n.unstable.eq.x , ncol = n.surfaces , byrow = F) #local phi values for each unstable equilibrium pair
for (i in 1:n.unstable.eq.x){
	x.loc <- round((unstable.xy[i,][1]-min(xlim))/x.range*mesh.xy[1])
	y.loc <- round((unstable.xy[i,][2]-min(ylim))/y.range*mesh.xy[2])
	unstable.phi.loc[i,] <- c(x.loc, y.loc)
	}

	for(i in 1: n.unstable.eq.x){
		for (j in 1:n.surfaces)
		unstable.phi[i,j] <- local.surfaces[[j]][unstable.phi.loc[i,1],unstable.phi.loc[i,2]]
	}

	max.phi <- which(unstable.phi == max(unstable.phi,na.rm= T),arr.ind=T)

	if(nrow(max.phi) != n.unstable.eq.x){ #if not all max(phi) are the same, then they need to be aligned
		global.max <- unstable.phi[max.phi[1,1],max.phi[1,2]]
		phi.diff <- matrix(data = NA , nrow = n.surfaces , ncol = 1)
		for (i in 1:n.surfaces){
			phi.diff[i,] <- global.max - unstable.phi[max.phi[1,1],i]
		}
	}

	eq.arr.unadj <- array(data = unlist(local.surfaces) , dim = c(mesh.xy[1] , mesh.xy[2] , n.surfaces))
	eq.arr <- array(data = NA , dim = c(mesh.xy[1] , mesh.xy[2] , n.surfaces))
	for(i in 1:n.surfaces){
		eq.arr[,,i] <- eq.arr.unadj[,,i]+phi.diff[1]
	}
	rm(eq.arr.unadj)
	gc()
	global.qp <- apply(eq.arr , c(1,2) , min)
	rm(eq.arr)
	gc()
	global.qp
}