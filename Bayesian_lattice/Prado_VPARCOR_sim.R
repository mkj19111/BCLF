
Prado_VPARCOR = function(casem, cases, Nsim, Dfactor)
{
tt = proc.time()

N = 1024
K = 2
p0 = 2
r1 = 0.1*1:N/N + 0.85
r2 = -0.1*1:N/N + 0.95
r3 = 0.2*1:N/N - 0.9
r4 = 0.2*1:N/N + 0.7
A0 = array(0,c(2,K,K,N))
A0[1,1,1,] = r1 * cos(2*pi/(15*1:N/N+5))
A0[1,2,2,] = r2 * cos(2*pi/(-10*1:N/N+15))


A0[2,1,1,] = -r1^2
A0[2,2,2,] = -r2^2

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
  S0[1,1,] = 1 + 1:N/N
  S0[2,2,] = 1 + 1:N/N
}

times <- 1:N
#freqs <- c(0,seq(0,0.49,by=0.01)+0.001,0.5)
freqs <- seq(0.001,0.499,by=0.01)
tspec <- tvar_spec_Guo(A0, S0, times, freqs)
tch = abs(tspec[1,2,,])^2/(abs(tspec[1,1,,])*abs(tspec[2,2,,]))
P=5
##### discount factors ######
grid_seq <- seq(Dfactor[1], Dfactor[2], Dfactor[3])
tmp_dim <- dim(as.matrix(expand.grid(grid_seq, grid_seq)))
delta_cpp <- array(dim = c(tmp_dim[1], 4, P))
## discount factor for PARCOR model
for(i in 1:P){
  delta_cpp[, , i] <- cbind(as.matrix(expand.grid(grid_seq, grid_seq)), 
                            as.matrix(expand.grid(grid_seq, grid_seq)))
}



set.seed(20000)
#set.seed(15535)

tt=proc.time()
P = 5
XX = array(0,c(Nsim,2,N))
A_m = array(0,c(Nsim,N,2,2,P))
S_m = array(0,c(Nsim,2,2))
ASE1 = numeric(Nsim)
ASE2 = numeric(Nsim)
CH = numeric(Nsim)
#######################################
#n.sim <- Nsim
#x.sim1 <- rep(list(NA), n.sim)
#n_t <- 1024
#n_pred <- 0
#
#for(i in 1:n.sim){
#  x.sim1[[i]] <- gen.sim.data(I = 2, n_t = n_t+n_pred, 
#                              phi1 = rep(0, n_t + n_pred),
#                              phi2 = rep(0, n_t + n_pred),
#                              phi3 = rep(0, n_t + n_pred)) 
#}
######################################

for (iter in 1:Nsim){   
  X = array(0,c(N+2,K))
  X[1,] = rnorm(K)
  X[2,] = rnorm(K)
  for (n in 1:N)
    X[n+2,] = A0[1,,,n]%*%X[n+1,] + A0[2,,,n]%*%X[n,]+ t(chol(S0[,,n]))%*%rnorm(K)

  X = X[p0+(1:N),]


XX[iter,,] = t(X)
}



for (iter in 1:Nsim){ 
X = XX[iter,,]

n_t <- 1024
n_pred <- 0
n_I <- 2
V1t <- diag(2)
G <- diag(4)
n_0 <- 1
S_0 <- V1t
### bandwidth of frequency ##
w <- seq(0.001, 0.499, 0.01)
n_w <- length(w)
##### prior parameters ######
mk_0 <- matrix(c(0, 0, 0, 0), ncol = 1)
Ck_0 <- diag(4)


result_sim <- run_parcor_parallel(F1 = X, G = G, mk_0 = mk_0, Ck_0 = Ck_0, 
                                          n_0 = n_0, S_0 = S_0, delta = delta_cpp, P = 5,
                                          sample_size = 500, chains = 5, DIC = TRUE)

P_opt <- which.min(round(result_sim$DIC, 3))

A_m[iter,,1,1,1:P_opt] = result_sim$ar_coef[[P_opt]]$forward[1,,]
A_m[iter,,2,1,1:P_opt] = result_sim$ar_coef[[P_opt]]$forward[2,,]
A_m[iter,,1,2,1:P_opt] = result_sim$ar_coef[[P_opt]]$forward[3,,]
A_m[iter,,2,2,1:P_opt] = result_sim$ar_coef[[P_opt]]$forward[4,,]

S_m[iter,,] = result_sim$St_fwd[[P_opt]][,,N-5]



spec <- tvar_spec_Guo_P(array(A_m[iter,,,,1:P_opt],c(1024,2,2,P_opt[iter])), S_m[iter,,], times, freqs)

#TT = 5:1019
TT = 3:1022
ASE1[iter] = mean( log(abs(spec[1,1,,TT])/abs(tspec[1,1,,TT]))^2 )
ASE2[iter] = mean( log(abs(spec[2,2,,TT])/abs(tspec[2,2,,TT]))^2 )
CH[iter]   =  mean((abs(spec[1,2,,TT])^2/(abs(spec[1,1,,TT])*abs(spec[2,2,,TT])) - tch[,TT])^2)
  

  print(iter)
}
#c(mean(ASE1),sd(ASE1))
#c(mean(ASE2),sd(ASE2))
#c(mean(CH),sd(CH))

print(proc.time()-tt)
ASEtable = cbind(ASE1,ASE2,CH)
ASE = rbind(c(mean(ASE1),sd(ASE1)),c(mean(ASE2),sd(ASE2)),c(mean(CH),sd(CH)))
return(list('ASE'=ASE,'ASEtable'=ASEtable,'A'=A_m,'S'=S_m,'A0'=A0,'S0'=S0))
}
