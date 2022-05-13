# Compute log-likelihood values of dynamic linear models
#
# Description:
#		
# *** Model:
#       obs equation: f_t = phi_t * b_t + e_t,   e_t ~ N(0,sig_t^2)
#
#     	trans equation: phi_t = phi_{t-1} + w_t,   w_t ~ N(0,W_t)
#
# *** Assumption:
#     	evolution error variance: sig_t^2 = (sig_{t-1}^2 * delta)/eta_t,   
#                         
#                              eta_t ~ Beta(delta*k_{t-1}/2, (1-delta)*k_{t-1}/2)
#
#     	evolution sys variance: W_t = C_{t-1}(1-alpha)/alpha
# 
# *** Initial prior:
#     	phi_0 ~ N(m_0, C_0)
#     	sig_0^{-2} ~ Gamma(k_0/2, d_0/2)
#
#	Usage:
#		dynamicLik(y, x, np, del, m0, C0, s0, n0) 
#
# Arguments:
#   y	 : 1 x N dependent variable
#   x  : 1 x N independent variable 
#     np : number of padding  head and tail
# 	del: Discount factors: del[1] = alpha and del[2] = delta  
#  	m0 : initial prior mean for state
#  	C0 :  initial variance for state
#  	s0 :  initial prior estimate of obs var
#  	n0 :  initial prior df
#
# Values: 
#  m    : 1 x N post mean 
#  s    : 1 x N post obs var estimates
#  e    : 1 x N estimated innovations 
#  llike: log-likelihood
#
# Note:
#		

MdynamicLik = function(y, x, np, nptype, del, m0, C0, s0, n0,tvcov,sm=0){


# set output

N <- length(y)
if (np>0){
  if (nptype==0){
    y = c(y[1:np],y,y[1:np + (N-np)])
    x = c(x[1:np],x,x[1:np + (N-np)])
  }else{
    yhead = y[1:np]
    ytail = y[1:np + (N-np)]
    y = c(rnorm(np,mean(yhead),sd(yhead)),y,rnorm(np,mean(ytail),sd(ytail)))
    xhead = x[1:np]
    xtail = x[1:np + (N-np)]
    x = c(rnorm(np,mean(xhead),sd(xhead)),x,rnorm(np,mean(xtail),sd(xtail))) 
  }
}

if (tvcov > 0 )
{
N <- length(y)
m <- numeric(N)  
#mn <- numeric(N) 
C <- numeric(N)  
s <- numeric(N)    
n <- numeric(N)  
d <- del[1] 
b <- del[2]
mt <- m0 
Ct <- C0 
st <- s0 
nt <- n0
llike <- 0
# forward filtering
for (ii in 1:N) {
  G <- x[ii]
  A <- Ct*G/d 
	Q <- G*A + st 
	A <- A/Q 
	e <- y[ii] - G * mt
	nt <- b * nt
	# non-standardized Student's t-distribution
   if ( (ii > np) && (ii< N+np+1) ) llike <- llike + lgamma((nt+1)/2) - lgamma(nt/2) - log(nt*Q)/2 - (nt+1)*log(1+e*e/(Q*nt))/2 
    
	mt <- mt + A*e 
	m[ii] <- mt
	r <- nt + e*e/Q 
	nt <- nt+1 
	r <- r/nt 
	st <- st*r 
	n[ii] <- nt 
	s[ii] <- st   
  Ct <- r * (Ct/d - A*A*Q)
	C[ii] <- Ct 
      Wt <- (1-d)*Ct/d 	
}

# backward smoothing
e <- numeric(N)  
#mn[N] <- m[N]
CC = C
for (ii in (N-1):1) {
  #B = CC[ii]/R[ii+1]
  m[ii] <- (1-d)*m[ii] + d*m[ii+1]
  #mn[ii] <- m[ii] + B*(mn[ii+1]-m[ii+1])
  e[ii] <- y[ii]- m[ii]*x[ii]
     
  n[ii] <- (1-b)*n[ii] + b*n[ii+1]  
  st <- s[ii] 
  s[ii] <- 1 / ((1-b)/st + b/s[ii+1]) 
  CC[ii] <- s[ii] * ((1-d)*CC[ii]/st + d*d*CC[ii+1]/st) 
  #CCn[ii] <- CC[ii] - B^2*(R[ii+1] - CCn[ii+1])
}
	idx = 1:(N-2*np)+np

if (sm==1){
n = 400
xx = m
xx = c(rep(xx[1],n/2),xx,rep(xx[N],n/2))
sa = filter(xx, rep(1 / n, n), sides = 2)
m=sa[n/2+1:N]

xx = CC
xx = c(rep(xx[1],n/2),xx,rep(xx[N],n/2))
sa = filter(xx, rep(1 / n, n), sides = 2)
CC=sa[n/2+1:N]

xx = s
xx = c(rep(xx[1],n/2),xx,rep(xx[N],n/2))
sa = filter(xx, rep(1 / n, n), sides = 2)
s=sa[n/2+1:N]

y = y
x = x
e = y - m*x
}
	return(list(m = m[idx], C = C[idx], CC = CC[idx], Wt = Wt, s = s[idx], e = e[idx], n = n[idx], llike = llike))
}else{

###############################################


N <- length(y)
m <- numeric(N)  
CC <- numeric(N) 
s <- 0   
n <- 0 
 

d <- del[1] 
mt <- m0 
Ct <- C0 
st <- s0 
nt <- 0
llike <- 0
sn = 0
n0 = 1
# forward filtering
for (ii in 1:N) {
  G <- x[ii]
  A <- Ct*G/d 
	Q <- G*A + st 
	A <- A/Q 
	e <- y[ii] - G * mt
	# non-standardized Student's t-distribution
   if ( (ii > np) && (ii< N+np+1) ) llike <- llike  
    
	mt <- mt + A*e 
	m[ii] <- mt
      nt <- nt + 1
      sn = sn + st*e^2/Q
#	r <- 1 + (e^2/Q - 1)/nt 
	st <- (n0*s0+sn)/(n0+ii)   
      s <- st
      Ct <- Ct/d - A*A*Q
#      Ct <- r * (Ct/d - A*A*Q)
	CC[ii] <- Ct 	
}

# backward smoothing
e <- numeric(N)  
#mn[N] <- m[N]
for (ii in (N-1):1) {
  m[ii] <- (1-d)*m[ii] + d*m[ii+1]
  e[ii] <- y[ii]- m[ii]*x[ii]
  CC[ii] <- (1-d)*CC[ii] + d*d*CC[ii+1]
}
	idx = 1:(N-2*np)+np
	return(list(m = m[idx], CC = CC[idx], s = rep(s,N-2*np), e = e[idx], n = rep(nt,N-2*np), llike = llike))
}
}