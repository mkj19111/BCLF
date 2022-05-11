

setwd('D:/research/Current/Bayesian_lattice')
rm(list = ls())
library("fields")
library(mvtnorm)
idx = function(k,K,N){ return(K*1:N-K+k) }

source("utility.r")
source("BayesLatticeLik_M.r")
source("BayesLattice_M.r")
source("dynamicLik.r")
source("tvar_spec.r")
source("tvar_spec_Guo.r")

source("MBayesLattice.r")
source("MdynamicLik.r")



##########




##########

ifl = read.csv('D:/research/Current/Pois_AR/inflation.csv',
    header=F,sep=',')
ifl=as.matrix(ifl[47:107,1:12])

ifl=as.numeric(t(ifl)[1:(61*12)])
t =  1
x1 = 1:(61*4)
for( i in 1:(61*4)){
    x1[t] = prod(1+ifl[1:3+3*(t-1)]/100)
    t = t + 1
}

unemployment = read.csv('D:/research/Current/Pois_AR/unemployment.csv',
    header=T,sep=',')
x2 = unemployment[,2]

interest = read.csv('D:/research/Current/Pois_AR/interest.csv',
    header=T,sep=',')
#313  1044 1960 - 2020
interest = interest[313:1044,2]
t =  1
x3 = 1:(61*4)
for( i in 1:(61*4)){
    x3[t] = mean(interest[1:3+3*(t-1)])
    t = t + 1
}

Xdata = cbind(x1,x2,x3)
#1:240   1960 - 2019
#1:168    1960 - 2001
#  1975:1  1981:3 1996:1
#   61     87    145
Xdata = Xdata[1:168,]
if (1){
  for (i in 1:3) {
    mx = mean(Xdata[,i])
    Xdata[,i] = (Xdata[,i] - mx )
  }
}
K = 3
N = 168
signal = t(Xdata)[1:(K*N)]


Pvar = 10
tvarp <- MBayesLattice(signal, K, Pvar, c(0.95,0.995,0.01) , 0,0,1,0,0)


MSresult = array(0,c(Pvar,3))
for(pvar in 1:Pvar){
  para = cparas(tvarp$coefs,tvarp$coefsC,tvarp$s2,tvarp$n,pvar)
  MSresult[pvar,] = MS(signal,para$A,para$S)[[1]]
}
P_sel = which(MSresult[,3] == min(MSresult[,3]))
P_sel







xdate = seq(as.Date("1960/1/1"), as.Date("2001/12/1"), by = "quarter")

par(mfrow=c(3,1))
plot(xdate, Xdata[,1],type='l',ylab='Inflation Rate',xlab='Date',cex.lab=1.5,cex.axis=1.5)
plot(xdate, Xdata[,2],type='l',ylab='Unemployment Rate',xlab='Date',cex.lab=1.5,cex.axis=1.5)
plot(xdate, Xdata[,3],type='l',ylab='Short-term Interest Rate',xlab='Date',cex.lab=1.5,cex.axis=1.5)





##########################
para = cparas(tvarp$coefs,tvarp$coefsC,tvarp$s2,tvarp$n,P_sel)
A = para$A
As = para$As
S = para$S 
Ss = para$Ss

xdate = seq(as.Date("1960/2/1"), as.Date("2001/12/1"), by = "quarter")


TT = 3:N
TTci = 30:N
p=1
par(mfrow=c(3,3))
for (j2 in 1:3){
    for (j1 in 1:3){
        plot(xdate[TT],A[p,j2,j1,TT],type='l',xlab='Date',ylab=bquote(Phi[.(j2)][.(j1)]),
           ylim=c(min(A[p,j2,j1,TTci]-1.96*sqrt(As[p,j2,j1,TTci])),max(A[p,j2,j1,TTci]+1.96*sqrt(As[p,j2,j1,TTci]))) )
        polygon(c(xdate[TT],rev(xdate[TT])) ,c(A[p,j2,j1,TT]-1.96*sqrt(As[p,j2,j1,TT]),rev(A[p,j2,j1,TT]+1.96*sqrt(As[p,j2,j1,TT]))),
            col = "grey75", border=F  )  
        lines(xdate[TT],A[p,j2,j1,TT])
       # lines(xdate[TT],A[p,j2,j1,TT]+1.96*sqrt(As[p,j2,j1,TT]),col='blue')
       # abline(v=as.Date("1987/2/1"),lty=4)
    }
}

p=2
par(mfrow=c(3,3))
for (j2 in 1:3){
    for (j1 in 1:3){
        plot(xdate[TT],A[p,j2,j1,TT],type='l',xlab='Date',ylab=bquote(Phi[.(j2)][.(j1)]),
           ylim=c(min(A[p,j2,j1,TTci]-1.96*sqrt(As[p,j2,j1,TTci])),max(A[p,j2,j1,TTci]+1.96*sqrt(As[p,j2,j1,TTci]))) )
        polygon(c(xdate[TT],rev(xdate[TT])) ,c(A[p,j2,j1,TT]-1.96*sqrt(As[p,j2,j1,TT]),rev(A[p,j2,j1,TT]+1.96*sqrt(As[p,j2,j1,TT]))),
            col = "grey75", border=F  )  
        lines(xdate[TT],A[p,j2,j1,TT])
       # lines(xdate[TT],A[p,j2,j1,TT]+1.96*sqrt(As[p,j2,j1,TT]),col='blue')
       # abline(v=as.Date("1987/2/1"),lty=4)
    }
}



par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(4,4.5,1,0.5))
for (j2 in 1:3){
    for (j1 in 1:3){
        plot(xdate[TT],S[j2,j1,TT],type='l',xlab='Date',ylab=bquote(S[.(j2)][.(j1)]),
            ylim=c(min(S[j2,j1,TT]-1.96*sqrt(Ss[j2,j1,TT])),max(S[j2,j1,TT]+1.96*sqrt(Ss[j2,j1,TT]))),
            cex.lab=1.5,cex.axis=1.5 )

       # lines(xdate[TT],S[j2,j1,TT]-1.96*sqrt(Ss[j2,j1,TT]),col='blue')
       # lines(xdate[TT],S[j2,j1,TT]+1.96*sqrt(Ss[j2,j1,TT]),col='blue')
        polygon(c(xdate[TT],rev(xdate[TT])) ,c(S[j2,j1,TT]-1.96*sqrt(Ss[j2,j1,TT]),rev(S[j2,j1,TT]+1.96*sqrt(Ss[j2,j1,TT]))),
            col = "grey75", border=F  )      
        lines(xdate[TT],S[j2,j1,TT])
    }
}

