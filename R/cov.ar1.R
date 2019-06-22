# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' estimates covariace matrix for AR1 parameters
#' @param ck number of cluster with each unique cluster size
#' @param nk unique cluster sizes
#' @param rho estimated rho
#' @param sigma2 estimated sigma2
#' @return a list with computed materials 
#' 
#' @author Vahid Nassiri
#' @export
cov.ar1 <- function(ck,nk,rho,sigma2){
	num.split=length(ck)
	var.mu1=(ck/(sigma2*(1-(rho^2))))* (((nk-2)*rho^2)-(2*((nk-1)*rho))+nk)
	var.mu=1/var.mu1
	w.mu=var.mu1/sum(var.mu1)
# Note that the unbiased version of the covariance is used here
	v22=2*(sigma2^2)*(1+(rho^2))
	v12=2*sigma2*(1-(rho^2))
	v11=(1-rho^2)^2
	var.varcomp1=matrix(c(v11,v12,v12,v22),2,2)
	varcomp.coef=1/(ck*(nk-((nk-2)*(rho^2))))
	var.varcomp=outer(var.varcomp1,varcomp.coef)
	W.total=0
	for (i in 1:num.split){
		W.total=W.total+solve(var.varcomp[,,i])
	}
	w.varcomp=array(0,c(2,2,num.split))
	for (i in 1:num.split){
		w.varcomp[,,i]=solve(W.total)%*%solve(var.varcomp[,,i])
	}
	return(list(var.mu=var.mu,var.varcomp=var.varcomp
					,w.mu=w.mu,w.varcomp=w.varcomp))
}





