# Estimate coefficients, innovation variances, and PAROCR coefficients
# of TVAR(P) 
#
#	Description:
# 
#	Usage:
#		BayesLattice(signal, Dfactor, np ,MS)
#
# Arguments:
#   signal :  1 x KN  dependent variable
#   Dfactor:  P x 2K matrix with columns: gamma and delta
#   np     :  number of paddings
#
# Values: 
#  coefs : P x N matrix of estimated coefficients of TVAR models 
#  s2  	 : 1 x N vector of estimated innovation variances  
#  parcor: P x N matrix of estimated PARCOR
#
#	Note:
#

MBayesLattice = function(signal, K, Pvar, Dfactor, np, nptype, tvcov, sm=0,MC=0,pred=0,ns=10){

N <- length(signal)/K
Ns <- 50

P <- Pvar*K + K -1
if ( K < 2 ) stop('Need more variables.')
#if ( P > nrow(Dfactor) ) stop('Need more tuning parameters in Dfactor.')

dfactor = seq(Dfactor[1],Dfactor[2],by=Dfactor[3])
ndf = length(dfactor)
Dfactors = array(0,c(P,2*K))

llik = rep(0,P)
lkid = 1

fcoefs <- array(0, dim = c(P,K,P,N))   #  stage, variable, lag, time
fcoefsC <- array(0, dim = c(P,K,P,N))
bcoefs <- array(0, dim = c(P,K,P,N))
bcoefsC <- array(0, dim = c(P,K,P,N))
s2 <- array(0, c(P,K,N))
n <- array(0, c(P,K,N))
fcp <- array(0, c(P,K,P,ns))
bcp <- array(0, c(P,K,P,ns))
fcpW <- array(0, c(P,K))
bcpW <- array(0, c(P,K))

# initial values at t_0
m0 <- 0 
C0 <- 1
n0 <- 1

finnov <- signal 
binnov <- signal

tfinnov <- signal 
tbinnov <- signal

for (p in 1:P) {
    cat("*** Order: ",p,"\n", sep="")

    flush.console()
    for (k in 1:K) {
        if( p <= P-(K-k) ){
          idxk = idx(k,K,N)
          ##  forward
          maxlik <- -1e20
          for (i in 1:ndf){
            for (j in 1:ndf){
              del <- c(dfactor[i],dfactor[j])
              s0 <- var(finnov[idxk[1:10]])
              #s0 <- 1
     #         resultfc <- MdynamicLik(finnov[idxk],c( numeric(nz), binnov[idxk[(nz+1):N]-p]), np, nptype, del, m0, C0, s0, n0, tvcov,sm)
     #   old
     #         resultfc <- MdynamicLik(finnov[idxk],c( binnov[idxk[1*(k==1)]], binnov[idxk[(1+1*(k==1)):N]-1]), np, nptype, del, m0, C0, s0, n0, tvcov)
     #  before nz
              resultfc <- MdynamicLik(finnov[idxk],c( numeric(1*(k==1)), binnov[idxk[(1+1*(k==1)):N]-1]), np, nptype, del, m0, C0, s0, n0, tvcov,sm)

              if (resultfc$llike > maxlik){
                maxlik <- resultfc$llike
                resultf <- resultfc
                Dfactors[p,1:2 + 2*(k-1)] <- del
              }
            }
          }
          llik[lkid] = maxlik
          lkid = lkid + 1
          ##  backward
          

          del <- Dfactors[p,1:2 + 2*(k-1)] 
          s0 <- var(binnov[idxk[1:10]])
          #s0 <- 1
    # before nz
    #      resultb <- MdynamicLik(c(numeric(1*(k==1)), binnov[idxk[(1+1*(k==1)):N]-1]),finnov[idxk], np, nptype, del, m0, C0, s0, n0, tvcov,sm)
          resultb <- MdynamicLik(c(numeric(1*(k==1)), binnov[idxk[(1+1*(k==1)):N]-1]),finnov[idxk], np, nptype, del, m0, C0, s0, n0, tvcov,sm)

	  
          fcoefs[p,k,p,] <- resultf$m
          fcoefsC[p,k,p,] <- resultf$CC
          fcpW[p,k] <- resultf$Wt
          bcoefs[p,k,p,] <- resultb$m 
          bcoefsC[p,k,p,] <- resultb$CC
          bcpW[p,k] <- resultb$Wt
          tfinnov[idxk] <- resultf$e
          tbinnov[idxk] <- resultb$e
          s2[p,k,] <- resultf$s
          n[p,k,] <- resultf$n

          if (p >1){
            if (k==1){
              k1 = K
            }else{
              k1 = k-1
            }
            for (jj in 1:(p-1)) {
	          hh <- p - jj
		    if (k==1) {
                    fcoefs[p,k,jj,] <- fcoefs[p-1,k,jj,] - fcoefs[p,k,p,] * c(0,bcoefs[p-1,k1,hh,-N])
                }else{
                    fcoefs[p,k,jj,] <- fcoefs[p-1,k,jj,] - fcoefs[p,k,p,] * bcoefs[p-1,k1,hh,]
                }
		    fcoefsC[p,k,jj,] <- fcoefsC[p,k,p,] * bcoefs[p-1,k1,hh,]^2
	      }
    
  	      for (jj in 1:(p-1)) {
	          hh <- p - jj
		    if (k==1) {
                    bcoefs[p,k,jj,] <- c(0,bcoefs[p-1,k1,jj,-N]) - bcoefs[p,k,p,] * fcoefs[p-1,k,hh,]
                }else{
                    bcoefs[p,k,jj,] <- bcoefs[p-1,k1,jj,] - bcoefs[p,k,p,] * fcoefs[p-1,k,hh,]
                }
	      }
          }

     }  #end if p <= P-(K-k) 
 
   } #end k loop    
	
   finnov <- tfinnov
   binnov <- tbinnov	
}  #end p loop

if (pred == 1){
  for (p in 1:P) {
    for (k in 1:K) {
      fcp[p,k,p,] = rnorm(ns,fcoefs[p,k,p,N], fcoefsC[p,k,p,N]+fcpW[p,k] )
      bcp[p,k,p,] = rnorm(ns,bcoefs[p,k,p,N-1], bcoefsC[p,k,p,N-1]+bcpW[p,k] )
         #################### Durbin-Levinson
         if (p >1){
            if (k==1){
              k1 = K
            }else{
              k1 = k-1
            }
            for (jj in 1:(p-1)) {
	          hh <- p - jj
                fcp[p,k,jj,] <- fcp[p-1,k,jj,] - fcp[p,k,p,] * bcp[p-1,k1,hh,]
	      }
    
  	      for (jj in 1:(p-1)) {
	          hh <- p - jj
                bcp[p,k,jj,] <- bcp[p-1,k1,jj,] - bcp[p,k,p,] * fcp[p-1,k,hh,]
	      }
         } 
         ###############
    }
  }
}
return(list(coefs = fcoefs, coefsC = fcoefsC, fcp = fcp, bcp = bcp, s2 = s2 ,n=n, Dfactors=Dfactors, MC=MC , llik= llik ))
}

#####################################################################

MBayesLattice2 = function(signal, K, Pvar, Dfactor, np, nptype, tvcov, sm=0,MC=0){

N <- length(signal)/K
Ns <- 50
P <- Pvar*K + K -1
if ( K < 2 ) stop('Need more variables.')
#if ( P > nrow(Dfactor) ) stop('Need more tuning parameters in Dfactor.')


fcoefs0 <- array(0, dim = c(K,P,N))
fcoefs1 <- array(0, dim = c(K,P,N))
fcoefsC <- array(0, dim = c(K,P,N))
bcoefs0 <- array(0, dim = c(K,P,N))
bcoefs1 <- array(0, dim = c(K,P,N))
s2 <- array(0, c(K,N))
# initial values at t_0
m0 <- 0 
C0 <- 10
n0 <- 1 

finnov <- signal 
binnov <- signal

tfinnov <- signal 
tbinnov <- signal
for (p in 1:P) {
    cat("*** Order: ",p,"\n", sep="")

    flush.console()
    for (k in 1:K) {
        if( p <= P-(K-k) ){
          idxk = idx(k,K,N)
          ##  forward
          del <- Dfactor
          s0 <- var(finnov[idxk[1:10]])
          resultf <- MdynamicLik(finnov[idxk],c( numeric(1*(k==1)), binnov[idxk[(1+1*(k==1)):N]-1]), np, nptype, del, m0, C0, s0, n0, tvcov,sm)

          ##  backward
          s0 <- var(binnov[idxk[1:10]])
          resultb <- MdynamicLik(c(numeric(1*(k==1)), binnov[idxk[(1+1*(k==1)):N]-1]),finnov[idxk], np, nptype, del, m0, C0, s0, n0, tvcov,sm)


	  
          fcoefs1[k,p,] <- resultf$m
          fcoefsC[k,p,] <- resultf$CC
          bcoefs1[k,p,] <- resultb$m 
          tfinnov[idxk] <- resultf$e
          tbinnov[idxk] <- resultb$e
          s2[k,] <- resultf$s

          if (p >1){
            if (k==1){
              k1 = K
            }else{
              k1 = k-1
            }
            for (jj in 1:(p-1)) {
	          hh <- p - jj
		    fcoefs1[k,jj,] <- fcoefs0[k,jj,] - fcoefs1[k,p,] * bcoefs0[k1,hh,]
		    fcoefsC[k,jj,] <- fcoefsC[k,p,] * bcoefs0[k1,hh,]^2
	      }
    
  	      for (jj in 1:(p-1)) {
	          hh <- p - jj
		    bcoefs1[k,jj,] <- bcoefs0[k1,jj,] - bcoefs1[k,p,] * fcoefs0[k,hh,]
	      }
          }

        
       }  #end if p <= P-(K-k) 
 
   } #end k loop    
	
   fcoefs0 <- fcoefs1
#   bcoefs0 <- bcoefs1
   bcoefs0[,,-1] <- bcoefs1[,,-N]
   bcoefs0[,,1] <- 0

   finnov <- tfinnov
   binnov <- tbinnov	
}  #end p loop


return(list(coefs = fcoefs1, coefsC = fcoefsC, s2 = s2,n=n))
}
##########################################################
########### no updated below  ############################

BayesLattice_MP = function(signal, K, Pvar, Dfactor, np, tvcov, MC=0){

N <- length(signal)/K
Ns <- 50
P <- Pvar*K + K -1
if ( K < 2 ) stop('Need more variables.')
if ( P > nrow(Dfactor) ) stop('Need more tuning parameters in Dfactor.')
aic = rep(0,Pvar)
bic = rep(0,Pvar)

mu <- array(0, dim = c(K,N))
muC <- array(0, dim = c(K,N))
fcoefs0 <- array(0, dim = c(K,P,N))
fcoefs1 <- array(0, dim = c(K,P,N))
fcoefsC <- array(0, dim = c(K,P,N))
bcoefs0 <- array(0, dim = c(K,P,N))
bcoefs1 <- array(0, dim = c(K,P,N))
coefsMS <- array(0, dim = c(Pvar,K,P,N))
s2MS <- array(0, dim = c(Pvar,K,N))
coefsSamples <- array(0, dim = c(Pvar,K,P,N,Ns))
s2Samples <- array(0, dim = c(Pvar,K,N,Ns))

s2 <- array(0, c(K,N))
n <- array(0, c(K,N))
# initial values at t_0
m0 <- 0.4 
C0 <- 10
n0 <- 1 

finnov <- signal 
binnov <- signal

tfinnov <- signal 
tbinnov <- signal
for (p in 1:P) {
    #cat("*** Order: ",p,"\n", sep="")

    flush.console()
    for (k in 1:K) {
        if( p <= P-(K-k) ){
          idxk = idx(k,K,N)
	    del <- Dfactor[p,1:2 + 2*(k-1)]

          if (p<0){
            s0 <- var(finnov[idxk[1:10]])
            resultf <- dynamicLik2(finnov[idxk],c( binnov[1*(k==1)], binnov[idxk[(1+1*(k==1)):N]-1]), np, del, c(mean(signal[1:N+k-1]),m0), diag(C0,2), s0, n0)
            s0 <- var(binnov[idxk[1:10]])         
            resultb <- dynamicLik( c(binnov[idxk[1*(k==1)]], binnov[idxk[(1+1*(k==1)):N]-1]),finnov[idxk], np, del, m0, C0, s0, n0, tvcov)
	  
            #resultb <- dynamicLik(binnov[idxk],c(finnov[idxk[1:(N-offsetb)]+p], numeric(offsetb)), np, del, m0, C0, s0, n0, tvcov)
            mu[k,] = resultf$m[,1]
            muC[k,] = resultf$CC[,1,1]
            fcoefs1[k,p,] <- -resultf$m[,2]
            fcoefsC[k,p,] <-  resultf$CC[,2,2]
            bcoefs1[k,p,] <- -resultb$m 
            tfinnov[idxk] <- resultf$e
            tbinnov[idxk] <- resultb$e
            s2[k,] <- resultf$s

          }else{
            s0 <- var(finnov[idxk[1:10]])
            resultf <- dynamicLik(finnov[idxk],c( binnov[1*(k==1)], binnov[idxk[(1+1*(k==1)):N]-1]), np, del, m0, C0, s0, n0, tvcov)
            s0 <- var(binnov[idxk[1:10]])         
            resultb <- dynamicLik( c(binnov[idxk[1*(k==1)]], binnov[idxk[(1+1*(k==1)):N]-1]),finnov[idxk], np, del, m0, C0, s0, n0, tvcov)
	  
            #resultb <- dynamicLik(binnov[idxk],c(finnov[idxk[1:(N-offsetb)]+p], numeric(offsetb)), np, del, m0, C0, s0, n0, tvcov)
            fcoefs1[k,p,] <- -resultf$m
            fcoefsC[k,p,] <-  resultf$CC
            bcoefs1[k,p,] <- -resultb$m 
            tfinnov[idxk] <- resultf$e
            tbinnov[idxk] <- resultb$e
            s2[k,] <- resultf$s
            n[k,] <- resultf$n
          }



          if (p >1){
            if (k==1){
              k1 = K
            }else{
              k1 = k-1
            }
            for (jj in 1:(p-1)) {
	          hh <- p - jj
		    fcoefs1[k,jj,] <- fcoefs0[k,jj,] + fcoefs1[k,p,] * bcoefs0[k1,hh,]
		    fcoefsC[k,jj,] <- fcoefsC[k,p,] * bcoefs0[k1,hh,]^2
	      }
    
  	      for (jj in 1:(p-1)) {
	          hh <- p - jj
		    bcoefs1[k,jj,] <- bcoefs0[k1,jj,] + bcoefs1[k,p,] * fcoefs0[k,hh,]
	      }
          }

          if ( ((p-k+1)%%K == 0) & (p>=K) ){
              im = (p-k+1)/K
              coefsMS[im,k,,] <- fcoefs1[k,,]   
              s2MS[im,k,] <- s2[k,]  
              if (MC==1){
                coefsSamples[im,k,p,,] <- array(-resultf$m+sqrt(resultf$CC)*rt(Ns*N,df=K),c(N,Ns))
                s2Samples[im,k,,] <- array(1/rgamma(n=Ns*N,resultf$n/2,(resultf$s*resultf$n)/2),c(N,Ns) )   
                
                if (k==1){
                  k1 = K
                }else{
                  k1 = k-1
                }
                for (jj in 1:(p-1)) {
	              hh <- p - jj
                    for (is in 1:Ns){
                        coefsSamples[im,k,jj,,is] <- fcoefs0[k,jj,] + coefsSamples[im,k,p,,is] * bcoefs0[k1,hh,]
                    }
	          }
              }
          }
       }  #end if p <= P-(K-k) 
 
   } #end k loop    
	
   fcoefs0 <- fcoefs1
   bcoefs0 <- bcoefs1
   finnov <- tfinnov
   binnov <- tbinnov	
}  #end p loop


return(list(mu = mu, muC= muC, coefs = fcoefs1, coefsC = fcoefsC, s2 = s2 , n = n, coefsMS = coefsMS, s2MS=s2MS, coefsSamples=coefsSamples,s2Samples=s2Samples,  MC=MC  ))
}



