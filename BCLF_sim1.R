
rm(list = ls())

# load packages
library("fields")
library("mvtnorm")

# load R functions
source("./Bayesian_lattice/utility.r")
source("./Bayesian_lattice/MBayesLattice.r")
source("./Bayesian_lattice/MdynamicLik.r")
source("./Bayesian_lattice/tvar_spectrum.r")
source("./Bayesian_lattice/sim1.r")

############################################
#  run simulation 1 of Sui(2022)
#  sim1(casem, cases, Nsim, Dfactor,flag_ASE, flag_DW) 
#  -----
#  Arguments: 
#  casem     1,2,3  case of autoregressive coefficients.
#  cases     1,2    case of innovation covariance.
#  Nsim      number of simulation datasets.
#  Dfactor   a sequence of values of discount factor.
#  flag_ASE  1 if ASE is computed, 0 otherwise.
#  flag_DW   1 if DIC and WAIC is computed, 0 otherwise.
#  -----
#  See readme for more details
############################################
Dfactor = c(0.990,1,0.001)
set.seed(20000)
result=sim1(1,1,Nsim=500,Dfactor,flag_ASE=1,flag_DW=0)  
result$ASE


###################################
# plots of parameter estimation   #
###################################
A0 = result$A0
S0 = result$S0
K=2
Pvar=2
N = 1024
A_pm = array(0,c(Pvar,K,K,N))
A_p05 = array(0,c(Pvar,K,K,N))
A_p95 = array(0,c(Pvar,K,K,N))
A_sd = array(0,c(Pvar,K,K,N))
S_pm = array(0,c(K,K,N))
S_p05 = array(0,c(K,K,N))
S_p95 = array(0,c(K,K,N))
for (i in 1:K){
  for (j in 1:K){
    for (p in 1:Pvar){
        A_pm[p,i,j,] = apply(result$A[p,i,j,,],1,mean)
        A_p05[p,i,j,] = apply(result$A[p,i,j,,],1,quantile,0.05)
        A_p95[p,i,j,] = apply(result$A[p,i,j,,],1,quantile,0.95)
        A_sd[p,i,j,] = apply(result$A[p,i,j,,],1,sd)
    }
        S_pm[i,j,] = apply(result$S[i,j,,],1,mean)
        S_p05[i,j,] = apply(result$S[i,j,,],1,quantile,0.05)
        S_p95[i,j,] = apply(result$S[i,j,,],1,quantile,0.95)
  }
}


par(mfrow=c(K,K))
for (k1 in 1:K){
    for (k2 in 1:K){
        plot(A_pm[1,k1,k2,],ylim=c(-1-1.5*((k1==1)*(k2==2)),1),type='l',
             xlab='Time',ylab=paste('Phi',k1,k2,sep=''))
        lines(A_p05[1,k1,k2,],col='blue')
        lines(A_p95[1,k1,k2,],col='blue')
        lines(A0[1,k1,k2,],col='red') 
    }
}

par(mfrow=c(K,K))
for (k1 in 1:K){
    for (k2 in 1:K){
        plot(A_pm[2,k1,k2,],ylim=c(-1,1),type='l',xlab='Time',
             ylab=paste('Phi',k1,k2,sep=''))
        lines(A_p05[2,k1,k2,],col='blue')
        lines(A_p95[2,k1,k2,],col='blue')
        lines(A0[2,k1,k2,],col='red') 
    }
}

par(mfrow=c(K,K))
for (k1 in 1:K){
    for (k2 in 1:K){
        plot(S_pm[k1,k2,],ylim=c(-0.2,5),type='l',xlab='Time',
             ylab=paste('S',k1,k2,sep=''))
        lines(S_p05[k1,k2,],col='blue')
        lines(S_p95[k1,k2,],col='blue')
        lines(S0[k1,k2,],col='red') 
    }
}

