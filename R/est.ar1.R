# TODO: Add comment
# 
# Author: Vahid Nassiri
###############################################################################


#' estimates parameters of a model with AR1 covariance structure
#' for balanced data
#' @param C number of clusters
#' @param n cluster sizes
#' @param Y the response
#' @param Plot if 1 (Default) will plot the third degree polynomial which is 
#' solved to find the estimates.
#' @return a list with estimated materials 
#' 
#' @author Vahid Nassiri
#' @export
est.ar1 <- function(C,n,Y,Plot=1){
# making a matrix out of the response vector
	Resp=matrix(Y,n,C)
# Computing cross products
	SS=crossprod(t(Resp))
# Computing S, \tilde{S} and R
	S=sum(diag(SS))
	S.tilde=sum(diag(SS)[2:(n-1)])
	tmp.R=SS
	diag(tmp.R)=NA
	tmp.R2 = (matrix(tmp.R[which(!is.na(tmp.R))],nrow=n,ncol=n-1))
	R=sum(tmp.R2[1,])
# Finding the coefficients of the 3rd degree polynomial and its roots
	P1=(n-1)*S.tilde
	P2=(n-2)*R
	P3=((n*S.tilde)+S)
	P4=n*R
	PP=polynomial(c(P4,-P3,-P2,P1))
	roots=polyroot(PP)
	Roots=Re(roots)[abs(Im(roots)) < 1e-6]
	rho.hat=Roots[abs(Roots)<1]
# Plotting the 3rd degree polynomial if requested
	if (Plot==1){
		plot(PP,xlim=c(-1.5,1.5),xlab="rho",ylab="3rd degree polynomial")
		abline(h=0,col=2)
		abline(v=Roots[abs(Roots)<1],lty=2,lwd=2,col=2)
		abline(v=-1,lty=2)
		abline(v=1,lty=2)
	}
# Estimating \sigma2
	tmp1=1/(C*n)
	tmp2=1/(1-(rho.hat^2))
	tmp3=S+((rho.hat^2)*S.tilde)
	tmp4=2*rho.hat*R
	sigma2.hat=(tmp1*tmp2)*(tmp3-tmp4)
	return(list(rho.hat=rho.hat,sigma2.hat=sigma2.hat))
}

