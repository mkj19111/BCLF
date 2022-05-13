

sim1 = function(casem, cases, Nsim, Dfactor, flag_ASE=1,flag_DW=0)
{
tt = proc.time()

np = 0  # no pre/pro-sequence padding 
tvcov = 1 # always choose time-varying innovation covariance

## MS bic   2  aic 5  995
##### simulate a bivariate series of N
#set.seed(5)
N = 1024
K = 2
p0 = 2
r1 = 0.1*1:N/N + 0.85
r2 = -0.1*1:N/N + 0.95
r3 = 0.2*1:N/N - 0.9
r4 = 0.2*1:N/N + 0.7
A0 = array(0,c(3,K,K,N))
A0[1,1,1,] = r1 * cos(2*pi/(15*1:N/N+5))
A0[1,2,2,] = r2 * cos(2*pi/(-10*1:N/N+15))

A0[2,1,1,] = -r1^2
A0[2,2,2,] = -r2^2
A0[1,1,2,] = 0

#A0[3,1,1,] = -r1^2/3
#A0[3,2,2,] = -r2^2/3

if (casem == 2){
  A0[1,1,2,] = -0.8
}else if (casem == 3){
  A0[1,1,2,] = r3  
  A0[2,1,2,] = r4
}

S0 = array(0,c(K,K,N))
S0[1,1,] = 1
S0[2,2,] = 1 

if (cases == 2){
#  S0[1,1,] =  4*c(1:(N/2),(N/2):1)/N 
#  S0[2,2,] =  4*c(1:(N/2),(N/2):1)/N
  S0[1,1,] = 1 + 1:N/N
  S0[2,2,] = 1 + 1:N/N
}
#############  true spectrum
if (flag_ASE == 1){
  times <- 1:N
 # freqs <- c(0,seq(0,0.49,by=0.01)+0.001,0.5)
  freqs <- seq(0.001,0.499,by=0.01)
  tspec <- tvar_spec_Guo(A0, S0, times, freqs)
  tch = abs(tspec[1,2,,])^2/(abs(tspec[1,1,,])*abs(tspec[2,2,,]))
}

Pvar = 5

P_sel = rep(0,Nsim)
P_DW = array(0,c(Nsim,2))

XX = array(0,c(Nsim,N*2))
ASE1 = rep(0,Nsim)
ASE2 = rep(0,Nsim)
CH   = rep(0,Nsim)
SP1 = array(0,c(50,1020))
SP2 = array(0,c(50,1020))
Ch = array(0,c(50,1020))

A_m = array(0,c(Pvar,K,K,N,Nsim))
S_m = array(0,c(K,K,N,Nsim))

Dfactors = array(0,c(Nsim,Pvar*K + K -1,4))

for (iter in 1:Nsim){   
#set.seed(20+iter)
X = array(0,c(N+2,K))
X[1,] = rnorm(K)
X[2,] = rnorm(K)
for (n in 1:N)
  X[n+2,] = A0[1,,,n]%*%X[n+1,] + A0[2,,,n]%*%X[n,]+ t(chol(S0[,,n]))%*%rnorm(K)

X = X[p0+(1:N),]



signal = 1:(K*N)
signal[idx(1,K,N)] = X[,1]
signal[idx(2,K,N)] = X[,2]

XX[iter,] = signal
}

for (iter in 1:Nsim){ 
signal = XX[iter,]
tvarp <- MBayesLattice(signal, K, Pvar, Dfactor , np,0,tvcov,0,0)
Dfactors[iter,,] = tvarp$Dfactors
#################### Model selection ###############
MSresult = array(0,c(Pvar,3))
for(pvar in 1:Pvar){
  para = cpara(tvarp$coefs,tvarp$s2,pvar)
  MSresult[pvar,] = MS(signal,para$A,para$S)[[1]]
}
P_sel[iter] = which(MSresult[,3] == min(MSresult[,3]))

if(flag_DW == 1){
ns = 20
lD = rep(0,Pvar)
pDIC = rep(0,Pvar)
DIC = rep(0,Pvar)
llikw = rep(0,Pvar)
WAIC = rep(0,Pvar)

likes = array(0,c(Pvar,ns,N))
for (j in 1:ns){

  tvarpdim = dim(tvarp$coefs)
  tvarpN = prod(tvarpdim)
  coefsp = array(rnorm(tvarpN,tvarp$coefs[1:tvarpN], sqrt(tvarp$coefsC[1:tvarpN])) ,tvarpdim)
 # coefsp = array(rnorm(tvarpN,tvarp$coefs[1:tvarpN], apply(rbind(sqrt(tvarp$coefsC[1:tvarpN]),rep(1,tvarpN)),2,min)) ,tvarpdim)
  s2sp = array(1/rgamma(prod(dim(tvarp$s2)),tvarp$n/2,tvarp$s2*tvarp$n/2),dim(tvarp$s2))
  for (pvar in 1:Pvar){
    para = cpara(coefsp,s2sp,pvar)
    ms = MS(signal,para$A,para$S)
    lD[pvar] = lD[pvar] + ms[[1]][1]
    likes[pvar,j,] = likes[pvar,j,] + ms[[2]]
  }
}
lD = lD/ns
likes = likes/ns


for (pvar in 1:Pvar){
  pDIC[pvar] = 2*(MSresult[pvar,1]-lD[pvar])
  DIC[pvar] = -2*MSresult[pvar,1] +2*pDIC[pvar]

  llikw[pvar] = sum(log(apply(likes[pvar,,(1+pvar):(N-pvar)],2,mean)))
  WAIC[pvar] = 2*llikw[pvar] - 2*lD[pvar]

}

P_DW[iter,1] = which(DIC == min(DIC))
P_DW[iter,2] = which(WAIC == min(WAIC))
}


if(flag_ASE == 1){

  para = cpara(tvarp$coefs,tvarp$s2,P_sel[iter])
  A_m[1:P_sel[iter],,,,iter] =  para$A
  S_m[,,,iter] =  para$S

  spec <- tvar_spec_Guo(para$A, para$S, 1:N, freqs)
  TT = 3:1022
  sp1 = log(abs(spec[1,1,,TT]))
  sp2 = log(abs(spec[2,2,,TT]))
  ch = abs(spec[1,2,,TT])^2/(abs(spec[1,1,,TT])*abs(spec[2,2,,TT]))
  SP1 = SP1 + sp1
  SP2 = SP2 + sp2
  Ch = Ch + ch

  ASE1[iter] = mean( (   sp1 -log(abs(tspec[1,1,,TT]))   )^2 )
  ASE2[iter] = mean( (   sp2 -log(abs(tspec[2,2,,TT]))   )^2 )
  CH[iter]   =  mean((ch - tch[,TT])^2)
}
  print(iter)
}
print(proc.time()-tt)
ASEtable = cbind(ASE1,ASE2,CH)
ASE = rbind(c(mean(ASE1),sd(ASE1)),c(mean(ASE2),sd(ASE2)),c(mean(CH),sd(CH)))
return(list('ASE'=ASE,'ASEtable'=ASEtable,'A'=A_m,'S'=S_m,'A0'=A0,'S0'=S0,'Dfactors'=Dfactors,
     'P_sel'=P_sel,'P_DW'=P_DW,'SP1'=SP1/Nsim,'SP2'=SP2/Nsim,'Ch'=Ch/Nsim))
}