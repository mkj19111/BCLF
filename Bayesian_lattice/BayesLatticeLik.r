# Search maximum log likelihood values at each lattice stage for TVAR(P) 
#
#	Description:
#
#	Usage:
# 	BayesLatticeLik(signal, P, disSys, disMea)
#
# Arguments:
#   signal: 1 x N dependent variable.
#   P     : maximum order of TVAR models.
#   disSys: a vector of discount factors for system errors (between 0 and 1).
#   disMea: a vector of discount factors for measurement errors (between 0 and 1).  
#
# Values: 
#  para_combin: a P x 4 matrix with columns: order, alpha, disMea, and loglikelihood.  
#
#	Note:
#

BayesLatticeLik = function(signal, P, disSys, disMea){

N <- length(signal)
likp <- array(0, dim = c(length(disSys),length(disMea),P))
para_combin <- matrix(0, nrow = P, ncol = 4)

# initial values at t_0
m0 <- 0.00001 
C0 <- 1
s0 <- var(signal[1:10])
n0 <- 1 

# start the searching
finnov <- signal 
binnov <- signal

for (kk in 1:P) {
	cat("*** Order: ",kk,"\n", sep="")
	flush.console()
	
	maxlike <- -1e300 # give a very small initial value for maxlike
	para_combin[kk,1] <- kk
	
	for (nSys in 1:length(disSys)) {
		for (nMea in  1 : length(disMea)) {
      del <- c(disSys[nSys],disMea[nMea])
			resultf <- dynamicLik(finnov,c(numeric(kk), binnov[1:(N-kk)]), 0,del, m0, C0, s0, n0,1)
			
			likp[nSys,nMea,kk] <- resultf$llike
			if (resultf$llike > maxlike) {
				para_combin[kk,2] <- disSys[nSys]
				para_combin[kk,3] <- disMea[nMea]
				para_combin[kk,4] <- resultf$llike
				maxlike <- resultf$llike;
			}
		}
	}
	del <- c(para_combin[kk,2], para_combin[kk,3])
	resultf <- dynamicLik(finnov,c(numeric(kk), binnov[1:(N-kk)]),0, del, m0, C0, s0, n0, 1)
  
	resultb <- dynamicLik(binnov,c(finnov[(kk+1):N], numeric(kk)), 0, del, m0, C0, s0, n0, 1)
		
	finnov <- resultf$e
	binnov <- resultb$e		
}
	return(para_combin)
}

