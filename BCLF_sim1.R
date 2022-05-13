
rm(list = ls())

# load packages
library("fields")
library("mvtnorm")

setwd('D:/Git/BCLF/')

# load R functions
source("./Bayesian_lattice/utility.r")
source("./Bayesian_lattice/MBayesLattice.r")
source("./Bayesian_lattice/MdynamicLik.r")
source("./Bayesian_lattice/tvar_spec_Guo.r")
source("./Bayesian_lattice/sim1.r")


Dfactor = c(0.990,1,0.001)
set.seed(20000)
result=sim1(1,1,Nsim=5,Dfactor,tvcov=1,flag_ASE=1,flag_DW=0)
result$ASE
#########
#> result$ASE
#             [,1]         [,2]
#[1,] 0.0356719037 0.0193232577
#[2,] 0.0382539704 0.0174810011
#[3,] 0.0009478236 0.0005074364



##########
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

