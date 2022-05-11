idx = function(k,K,N){ return(K*1:N-K+k) }

cpara = function(coefs,s2,Pvar){

N = dim(coefs)[4]
K = dim(coefs)[2]
P <- Pvar*K + K -1

A = array(0,c(Pvar,K,K,N))
S = array(0,c(K,K,N))
  for (t in 1:N)
  {
      L = diag(K)
      for (k in 2:K) L[k,1:(k-1)] = -coefs[P-(K-k),k,(k-1):1,t]
      Linv = solve(L)
      Ap = array(0,c(Pvar,K,K))
      for (p in 1:Pvar){
          for (k in 1:K){
              Ap[p,k,] = coefs[P-(K-k),k,(K*p+k-1):(K*(p-1)+k),t]
          }
          A[p,,,t] = Linv%*%Ap[p,,]
      }

      Sd = diag(K)
      for (k in 1:K) Sd[k,k] = s2[P-(K-k),k,t]
      S[,,t] = Linv%*%Sd%*%t(Linv)
  }
return(list(A = A, S = S))
}


cparas = function(coefs,coefsC, s2, n, Pvar){

N = dim(coefs)[4]
K = dim(coefs)[2]
P <- Pvar*K + K -1

A = array(0,c(Pvar,K,K,N))
As = array(0,c(Pvar,K,K,N))
S = array(0,c(K,K,N))
Ss = array(0,c(K,K,N))
  for (t in 1:N)
  {
      L = diag(K)
      for (k in 2:K) L[k,1:(k-1)] = -coefs[P-(K-k),k,(k-1):1,t]
      Linv = solve(L)
      Ap = array(0,c(Pvar,K,K))
      Aps = array(0,c(Pvar,K,K))
      for (p in 1:Pvar){
          for (k in 1:K){
              Ap[p,k,] = coefs[P-(K-k),k,(K*p+k-1):(K*(p-1)+k),t]
              Aps[p,k,] = coefsC[P-(K-k),k,(K*p+k-1):(K*(p-1)+k),t]
          }
          A[p,,,t] = Linv%*%Ap[p,,]
          As[p,,,t] = abs(Linv)%*%Aps[p,,]%*%t(abs(Linv))
      }

      Sd = diag(K)
      Sds = diag(K)
      for (k in 1:K){
          Sd[k,k] = s2[P-(K-k),k,t]
          Sds[k,k] = s2[P-(K-k),k,t]^2/( (n[P-(K-k),k,t]-1)^2*(n[P-(K-k),k,t]-2)  )
      }
      S[,,t] = Linv%*%Sd%*%t(Linv)
      Ss[,,t] = abs(Linv)%*%Sds%*%t(abs(Linv))
  }
return(list(A = A, As = As, S = S, Ss = Ss))
}

cpara_bk = function(coefs,s2,Pvar){
N = dim(coefs)[3]
K = dim(coefs)[1]

A = array(0,c(Pvar,K,K,N))
S = array(0,c(K,K,N))
  for (t in 1:N)
  {
      L = diag(K)
      for (k in 2:K) L[k,1:(k-1)] = -coefs[k,(k-1):1,t]
      Linv = solve(L)
      Ap = array(0,c(Pvar,K,K))
      for (p in 1:Pvar){
          for (k in 1:K){
              Ap[p,k,] = coefs[k,(K*p+k-1):(K*(p-1)+k),t]
          }
          A[p,,,t] = Linv%*%Ap[p,,]
      }

      Sd = diag(K)
      for (k in 1:K) Sd[k,k] = s2[k,t]
      S[,,t] = Linv%*%Sd%*%t(Linv)
  }
return(list(A = A, S = S))
}



MPpara = function(mu,coefs,s2,Pvar){
N = dim(coefs)[3]
K = dim(coefs)[1]

Mu = array(0,c(K,N))
A = array(0,c(Pvar,K,K,N))
S = array(0,c(K,K,N))
  for (t in 1:N)
  {
      L = diag(K)
      for (k in 2:K) L[k,1:(k-1)] = -coefs[k,(k-1):1,t]
      Linv = solve(L)
      Mu[,t] = Linv%*%mu[,t]
      Ap = array(0,c(Pvar,K,K))
      for (p in 1:Pvar){
          for (k in 1:K){
              Ap[p,k,] = -coefs[k,(K*p+k-1):(K*(p-1)+k),t]
          }
          A[p,,,t] = Linv%*%Ap[p,,]
      }

      Sd = diag(K)
      for (k in 1:K) Sd[k,k] = s2[k,t]
      S[,,t] = Linv%*%Sd%*%t(Linv)
  }
return(list(Mu = Mu, A = A, S = S))
}


TVVARMS = function(signal,Pvar,tvarp){
coefsMS = tvarp$coefsMS
s2MS = tvarp$s2MS
coefsSamples = tvarp$coefsSamples
s2Samples = tvarp$s2Samples

if(tvarp$MC == 1 ) Ns = dim(coefsSamples)[5]
N = dim(s2MS)[3]
K = dim(s2MS)[2]
X = rbind(array(0,c(Pvar,K)),t(array(signal,c(K,N))))
order = 1:Pvar
llik = rep(0,Pvar)
aic = rep(0,Pvar)
bic = rep(0,Pvar)
dic = rep(0,Pvar)
dic2 = rep(0,Pvar)
waic = rep(0,Pvar)
Pdic = 0
Pdic2 = 0
Pwaic = 0
for (P in 1:Pvar){
    if(tvarp$MC == 1 ){
        LSLik = rep(0,Ns)
        SLik = array(0,c(N,Ns))
    }

    paraMS = cpara(coefsMS[P,,,],s2MS[P,,],P)
    llike = 0
    for (t in 1:N){
        resid = X[t+Pvar,]
        for (p in 1:P){
            resid = resid - paraMS$A[p,,,t]%*%X[t+Pvar-p,]
        }
        llike = llike + dmvnorm(resid[1:K], mean=rep(0,K), sigma=paraMS$S[,,t], log = TRUE)
    }
    llik[P] = llike

    aic[P] = -2*llike + (P*K + K*(K-1)/2)*2
    bic[P] = -2*llike + (P*K + K*(K-1)/2)*log(N*K)


  if(tvarp$MC == 1 ){
    for (ss in 1:Ns){
        paraMS = cpara(coefsSamples[P,,,,ss],s2Samples[P,,,ss],P)
        for (t in 1:N){
            resid = X[t+Pvar,]
            for (p in 1:P){
                resid = resid - paraMS$A[p,,,t]%*%X[t+Pvar-p,]
            }
            #LSLik[ss] = LSLik[ss] + dmvt(resid[1:K], delta=rep(0,K), sig=paraMS$S[,,t], df=K,log = TRUE)    
            #SLik[t,ss] =  dmvt(resid[1:K], delta=rep(0,K), sig=paraMS$S[,,t], df=K,log = FALSE)             
            LSLik[ss] = LSLik[ss] + dmvnorm(resid[1:K], mean=rep(0,K), sigma=paraMS$S[,,t],log = TRUE)    
            SLik[t,ss] =  dmvnorm(resid[1:K], mean=rep(0,K), sigma=paraMS$S[,,t],log = FALSE)  
        }
    }   
    Pdic = 2*(llike - mean(LSLik)) + Pdic 
    dic[P]= -2*llike + 2*Pdic 
    Pdic2 = 2*var(LSLik) + Pdic2 
    dic2[P] = -2*llike + 2*Pdic2
    Pwaic = 2*( sum(log(apply(SLik,1,mean))) - mean(LSLik) ) + Pwaic
    waic[P] = -2*llike + 2*Pwaic
    
  }
}
######### model selection


return(cbind(order,llik,aic,bic,dic,dic2,waic))
}


########## model selection ##########
MS = function(signal,A,S){
P = dim(A)[1]
N = dim(A)[4]
K = dim(A)[2]
resid = array(0,c(N,K))

X = t(array(signal,c(K,N)))

  llike = 0
  like = rep(0,N)
    for (t in (1+P):(N-P)){
        resid[t,] = X[t,]
        for (p in 1:P){
            resid[t,] = resid[t,] - A[p,,,t]%*%X[t-p,]
        }
        like[t] = dmvnorm(resid[t,1:K], mean=rep(0,K), sigma=S[,,t], log = F)
        llike = llike + log(like[t])
    }


    #AIC = -2*llike + (P+1)*K^2 *2
    #BIC = -2*llike + (P+1)*K^2*log(N*K)
    AIC = -2*llike + (P*K^2 + K*(K-1)/2)*2 *2
    BIC = -2*llike + (P*K^2 + K*(K-1)/2)*2* log(N*K)

return(list(c(llike,AIC , BIC),like))
}



