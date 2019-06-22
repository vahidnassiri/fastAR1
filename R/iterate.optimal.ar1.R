# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' computing iterated optimal weights
#' @param ck number of cluster with each unique cluster size
#' @param nk unique cluster sizes
#' @param mu.split estimated \eqn{\mu} from different splits 
#' @param var.comp.split estimated variance components (rho and sigma^2)
#' from different splits
#' @param tol the tolerance
#' @return a list with computed materials 
#' 
#' @author Vahid Nassiri
#' @export
iterate.optimal.ar1 <- function(ck,nk,mu.split,var.comp.split,tol){
	num.split=length(ck)
	W.size.prop=(nk*ck)/sum(nk*ck)
	rho.hat= sum(var.comp.split[1,]*W.size.prop)
	sigma2.hat=sum(var.comp.split[2,]*W.size.prop)
	diff=10
	var.comp.hat=matrix(c(rho.hat,sigma2.hat),2,1)
	count=0
	while (diff>tol){
		WW=cov.ar1 (ck,nk,var.comp.hat[1,1],var.comp.hat[2,1])
		W.mu=WW$w.mu
		W.varcomp=WW$w.varcomp
		var.comp.hat=0
		for (i in 1:num.split){
			var.comp.hat=var.comp.hat+W.varcomp[,,i]%*%var.comp.split[,i]
		}
		W.mu.old=W.mu
		W.varcomp.old=W.varcomp
		WW=cov.ar1 (ck,nk,var.comp.hat[1,1],var.comp.hat[2,1])
		W.mu=WW$w.mu
		W.varcomp=WW$w.varcomp
		diff1=norm(as.matrix(W.mu-W.mu.old))
		diff2=sum(apply(W.varcomp-W.varcomp.old,3,norm))
		diff=max(c(diff1,diff2))
		count=count+1
	}
	var.comp.hat=0
	for (i in 1:num.split){
		var.comp.hat=var.comp.hat+W.varcomp[,,i]%*%var.comp.split[,i]
	}
	W.total=0
	WW=cov.ar1 (ck,nk,var.comp.hat[1,1],var.comp.hat[2,1])
	for (i in 1:num.split){
		W.total=W.total+solve(WW$var.varcomp[,,i])
	}
	mu.hat=sum(W.mu*mu.split)
	return(list(W.mu=W.mu,W.varcomp=W.varcomp,mu.hat=mu.hat
					,varcomp.hat=var.comp.hat, var.mu.hat=1/sum(1/WW$var.mu),
					var.varcomp.hat=solve(W.total),num.iterate=count))
}


