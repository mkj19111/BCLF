# Estimate coefficients, innovation variances, and PAROCR coefficients
# of TVAR(P) 
#
#	Description:
# 
#	Usage:
#		BayesLattice(signal, Dfactor)
#
# Arguments:
#   signal :  1 x N vector dependent variable
#   Dfactor:  P x 2 matrix with columns: gamma and delta
#
# Values: 
#  coefs : P x N matrix of estimated coefficients of TVAR models 
#  s2  	 : 1 x N vector of estimated innovation variances  
#  parcor: P x N matrix of estimated PARCOR
#
#	Note:
#

BayesLattice = function(signal, Dfactor){

N <- length(signal)
P <- nrow(Dfactor)

fcoefs <- array(0, dim = c(P,P,N))
fcoefsC <- array(0, dim = c(P,P,N))
bcoefs <- array(0, dim = c(P,P,N))

# initial values at t_0
m0 <- 0.4 
C0 <- 100
s0 <- var(signal[1:10]) 
n0 <- 1 

finnov <- signal 
binnov <- signal
del <- Dfactor[1,]

cat("*** Order: 1\n", sep="")
flush.console()
	
# get the instantaneous forward PARCOR
resultf <- dynamicLik(finnov, c(0, binnov[1:(N-1)]), 0, del, m0, C0, s0, n0, 1)

# get the instantaneous back PARCOR
resultb <- dynamicLik(binnov,c(finnov[2:N], 0), 0, del, m0, C0, s0, n0, 1);

fcoefs[1,1,] <- resultf$m
fcoefsC[1,1,]  <- resultf$CC
bcoefs[1,1,] <- resultb$m 


if (P > 1) {
	finnov <- resultf$e
	binnov <- resultb$e
	
	for (kk in 2:P) {
  	cat("*** Order: ",kk,"\n", sep="")
		flush.console()
	
		del <- Dfactor[kk,]

		# get the instantaneous forward PARCOR
		resultf <- dynamicLik(finnov,c(numeric(kk), binnov[1:(N-kk)]), 0,del, m0, C0, s0, n0,1)
		
		# get the instantaneous back PARCOR
		resultb <- dynamicLik(binnov,c(finnov[(kk+1):N], numeric(kk)), 0,del, m0, C0, s0, n0,1)
		
		fcoefs[kk,kk,] <- resultf$m
            fcoefsC[kk,kk,]  <- resultf$CC
		bcoefs[kk,kk,] <- resultb$m
    
		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			fcoefs[kk,jj,] <- fcoefs[ii,jj,] - fcoefs[kk,kk,] * bcoefs[ii,hh,]
                  fcoefsC[kk,jj,] <- fcoefsC[kk,kk,] * bcoefs[ii,hh,]^2
		}
    
		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			bcoefs[kk,jj,] <- bcoefs[ii,jj,] - bcoefs[kk,kk,] * fcoefs[ii,hh,]
		}
    
		finnov <- resultf$e
		binnov <- resultb$e
  }
}

s2 <- resultf$s
n <- resultf$n
parcor <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    parcor[ii,] <- fcoefs[ii,ii,]
}

coefs <- matrix(0, nrow = P, ncol = N)
coefsC <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    coefs[ii,] <- fcoefs[P,ii,]
    coefsC[ii,] <- fcoefsC[P,ii,]
}

return(list(coefs = coefs, coefsC = coefsC, parcor = parcor, n = n, s2 = s2))
}

###  Poisson TV-AR

BayesLattice_P = function(signal, Dfactor){

N <- length(signal)
P <- nrow(Dfactor)

fcoefs <- array(0, dim = c(P,P,N))
fcoefsC <- array(0, dim = c(P,P,N))
bcoefs <- array(0, dim = c(P,P,N))

# initial values at t_0
m0 <- 0.4 
C0 <- 100
s0 <- var(signal[1:10]) 
n0 <- 1 

finnov <- signal 
binnov <- signal
del <- Dfactor[1,]

cat("*** Order: 1\n", sep="")
flush.console()
	
# get the instantaneous forward PARCOR
resultf <- dynamicLik(finnov, c(0, binnov[1:(N-1)]), 0, del, m0, C0, s0, n0, 1)

# get the instantaneous back PARCOR
resultb <- dynamicLik(binnov,c(finnov[2:N], 0), 0, del, m0, C0, s0, n0, 1);

fcoefs[1,1,] <- resultf$m
fcoefsC[1,1,]  <- resultf$CC
bcoefs[1,1,] <- resultb$m 


if (P > 1) {
	finnov <- resultf$e
	binnov <- resultb$e
	
	for (kk in 2:P) {
  	cat("*** Order: ",kk,"\n", sep="")
		flush.console()
	
		del <- Dfactor[kk,]

		# get the instantaneous forward PARCOR
		resultf <- dynamicLik(finnov,c(numeric(kk), binnov[1:(N-kk)]), 0,del, m0, C0, s0, n0,1)
		
		# get the instantaneous back PARCOR
		resultb <- dynamicLik(binnov,c(finnov[(kk+1):N], numeric(kk)), 0,del, m0, C0, s0, n0,1)
		
		fcoefs[kk,kk,] <- resultf$m
            fcoefsC[kk,kk,]  <- resultf$CC
		bcoefs[kk,kk,] <- resultb$m
                      

		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			fcoefs[kk,jj,] <- fcoefs[ii,jj,] - fcoefs[kk,kk,] * bcoefs[ii,hh,]
                  fcoefsC[kk,jj,] <- fcoefsC[kk,kk,] * bcoefs[ii,hh,]^2
		}
    
		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			bcoefs[kk,jj,] <- bcoefs[ii,jj,] - bcoefs[kk,kk,] * fcoefs[ii,hh,]
		}
    
		finnov <- resultf$e
		binnov <- resultb$e
  }
}

s2 <- resultf$s
n <- resultf$n

parcor <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    parcor[ii,] <- fcoefs[ii,ii,]
}

coefs <- matrix(0, nrow = P, ncol = N)
coefsC <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    coefs[ii,] <- fcoefs[P,ii,]
    coefsC[ii,] <- fcoefsC[P,ii,]
}

return(list(coefs = coefs, coefsC = coefsC, parcor = parcor, n = n, s2 = s2))
}


###  Poisson TV-AR fixed tau2

BayesLattice_P_const = function(signal, Dfactor){

N <- length(signal)
P <- nrow(Dfactor)

fcoefs <- array(0, dim = c(P,P,N))
fcoefsC <- array(0, dim = c(P,P,N))
bcoefs <- array(0, dim = c(P,P,N))

# initial values at t_0
m0 <- 0.4 
C0 <- 100
s0 <- var(signal[1:10]) 
n0 <- 1 

finnov <- signal 
binnov <- signal
del <- Dfactor[1,]

cat("*** Order: 1\n", sep="")
flush.console()
	
# get the instantaneous forward PARCOR
resultf <- dynamicLik_c(finnov, c(0, binnov[1:(N-1)]), 0, del, m0, C0, s0, n0, 1)

# get the instantaneous back PARCOR
resultb <- dynamicLik_c(binnov,c(finnov[2:N], 0), 0, del, m0, C0, s0, n0, 1);

fcoefs[1,1,] <- resultf$m
fcoefsC[1,1,]  <- resultf$CC
bcoefs[1,1,] <- resultb$m 


if (P > 1) {
	finnov <- resultf$e
	binnov <- resultb$e
	
	for (kk in 2:P) {
  	cat("*** Order: ",kk,"\n", sep="")
		flush.console()
	
		del <- Dfactor[kk,]

		# get the instantaneous forward PARCOR
		resultf <- dynamicLik_c(finnov,c(numeric(kk), binnov[1:(N-kk)]), 0,del, m0, C0, s0, n0,1)
		
		# get the instantaneous back PARCOR
		resultb <- dynamicLik_c(binnov,c(finnov[(kk+1):N], numeric(kk)), 0,del, m0, C0, s0, n0,1)
		
		fcoefs[kk,kk,] <- resultf$m
            fcoefsC[kk,kk,]  <- resultf$CC
		bcoefs[kk,kk,] <- resultb$m
                      

		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			fcoefs[kk,jj,] <- fcoefs[ii,jj,] - fcoefs[kk,kk,] * bcoefs[ii,hh,]
                  fcoefsC[kk,jj,] <- fcoefsC[kk,kk,] * bcoefs[ii,hh,]^2
		}
    
		ii <- kk - 1
		for (jj in 1:(kk-1)) {
			hh <- kk - jj
			bcoefs[kk,jj,] <- bcoefs[ii,jj,] - bcoefs[kk,kk,] * fcoefs[ii,hh,]
		}
    
		finnov <- resultf$e
		binnov <- resultb$e
  }
}

s2 <- resultf$s
n <- resultf$n

parcor <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    parcor[ii,] <- fcoefs[ii,ii,]
}

coefs <- matrix(0, nrow = P, ncol = N)
coefsC <- matrix(0, nrow = P, ncol = N)
for (ii in 1:P) {
    coefs[ii,] <- fcoefs[P,ii,]
    coefsC[ii,] <- fcoefsC[P,ii,]
}

return(list(coefs = coefs, coefsC = coefsC, parcor = parcor, n = n, s2 = s2))
}

