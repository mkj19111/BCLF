

setwd('D:/research/Current/Bayesian_lattice')
rm(list = ls())
library("fields")
library("expm")
library(mvtnorm)

source("utility.r")
#source("BayesLatticeLik_P.r")
source("MBayesLattice.r")
source("MdynamicLik.r")
source("tvar_spec_Guo.r")


## transformation from direction and speed to X,Y component
trans <- function(sknt, drct){
  n_t <- length(sknt)
  y_1 <- sknt*sin(pi*drct/180) # X component
  y_2 <- sknt*cos(pi*drct/180) # Y component
  return(rbind(y_1, y_2))
}

### import Monterey Bay data 
asos.py <- read.csv("D:/research/Current/asos.py_MRY.txt", stringsAsFactors=FALSE)
dtimes <- asos.py$valid
dtimes <- as.POSIXlt(dtimes)
asos.py$drct <- as.numeric(asos.py$drct)
asos.py$sknt <- as.numeric(asos.py$sknt)
asos.py$sep <- cut(dtimes, breaks = "4 hours")

## the median of every four hours
tmp_sknt <- by(data = asos.py$sknt, INDICES = asos.py$sep, FUN = median, na.rm = TRUE)
tmp_drct <- by(data = asos.py$drct, INDICES = asos.py$sep, FUN = median, na.rm = TRUE)
data_median_MRY <- trans(tmp_sknt, tmp_drct)

### import Salinas data 
asos.py <- read.csv("D:/research/Current/asos.py_SNS.txt", stringsAsFactors=FALSE)
dtimes <- asos.py$valid
dtimes <- as.POSIXlt(dtimes)
asos.py$drct <- as.numeric(asos.py$drct)
asos.py$sknt <- as.numeric(asos.py$sknt)
asos.py$sep <- cut(dtimes, breaks = "4 hours")

## the median of every four hours
tmp_sknt <- by(data = asos.py$sknt, INDICES = asos.py$sep, FUN = median, na.rm = TRUE)
tmp_drct <- by(data = asos.py$drct, INDICES = asos.py$sep, FUN = median, na.rm = TRUE)
data_median_SNS <- trans(tmp_sknt, tmp_drct)

### import Watsonville data 
asos.py <- read.csv("D:/research/Current/asos.py_WVI.txt", stringsAsFactors=FALSE)
dtimes <- asos.py$valid
dtimes <- as.POSIXlt(dtimes)
asos.py$drct <- as.numeric(asos.py$drct)
asos.py$sknt <- as.numeric(asos.py$sknt)
asos.py$sep <- cut(dtimes, breaks = "4 hours")

## the median of every four hours
tmp_sknt <- by(data = asos.py$sknt, INDICES = asos.py$sep, FUN = median, na.rm = TRUE)
tmp_drct <- by(data = asos.py$drct, INDICES = asos.py$sep, FUN = median, na.rm = TRUE)
data_median_WVI <- trans(tmp_sknt, tmp_drct)


data_median <- rbind(data_median_MRY,data_median_SNS,data_median_WVI)
#for (j in 1:6)  data_median[j,] = data_median[j,] - mean(data_median[j,])
    
data_median1 <- data_median[c(1,3,5),]
data_median2 <- data_median[c(2,4,6),]

data_median <- as.matrix(rbind(data_median1,data_median2))


#for (j in 1:6)
#    data_median[j,] = data_median[j,] - mean(data_median[j,])



K = 6
N = 450
signal = data_median[1:(6*450)]


Pvar = 10
Dfactor = c(0.95,0.99,0.005)
tvarp <- MBayesLattice(signal, K, Pvar, Dfactor , 0,0,1,0,0)

#################### Model selection ###############
MSresult = array(0,c(Pvar,3))
for(pvar in 1:Pvar){
  para = cpara(tvarp$coefs,tvarp$s2,pvar)
  MSresult[pvar,] = MS(signal,para$A,para$S)[[1]]
}
P_sel = which(MSresult[,3] == min(MSresult[,3]))

P_sel = 6
para = cpara(tvarp$coefs,tvarp$s2,P_sel)
A = para$A
S = para$S



  times <- 1:N
  freqs <- seq(0.001,0.499,by=0.01)

  spec <- tvar_spec_Guo(para$A, para$S, times, freqs)

  #sp1 = log(abs(spec[1,1,,TT]))
  #sp2 = log(abs(spec[2,2,,TT]))
  #ch = abs(spec[1,2,,TT])^2/(abs(spec[1,1,,TT])*abs(spec[2,2,,TT]))

par(mfrow=c(2,3),oma=c( 0,0,0,2.5),mar=c(4.5,4.5,2,2))
for (j in 1:6){ 
ss = t(log(abs(spec[j,j,,])))
image(times, freqs, ss, main="Estimated Spectrum 1", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
}


par(mfrow=c(3,3),oma=c( 1.2,0,0,2.5),mar=c(4.5,5.5,2,2))
j1 = 1
j2 = 4
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 2
j2 = 5
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 3
j2 = 6
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 1
j2 = 2
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 2
j2 = 3
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 1
j2 = 3
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 4
j2 = 5
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 5
j2 = 6
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)
j1 = 4
j2 = 6
ss = t( abs(spec[j1,j2,,])^2/(abs(spec[j1,j1,,])*abs(spec[j2,j2,,])) )
image(times, freqs, ss, main="Estimated Coherence", xlab = "Time", ylab = "Frequency", col = tim.colors(64))
image.plot(ss, legend.shrink = 1, legend.width = 1,  col = tim.colors(64), legend.only = T, legend.mar = 0.5)




