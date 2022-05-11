# Calculated the time-varying spectrum of a TVAR(p) model
#
#	Description:
# 
#	Usage:
#		tvar_spec(m, s, times, freqs)
#
# Arguments:
#   m    : p x N matrix of posterior means of model
#   s    : 1 x N vector of posterior means of innovation variances
#   times: 1 x N vector of time points selected 
#   freqs: 1 x f vector of frequency selected  (0 to 0.5)
#
# Values: 
#   spec : f x T matrix of spectra
#
#	Note:
#

tvar_spec_Guo = function(m, s, times, freqs) {
nt <- length(times)
nf <- length(freqs)
P <- dim(m)[1]
K <- dim(m)[2]
if (length(dim(s))==3){
    S = s
}else{
    S = array(s,c(dim(s)[1:2],nt))
}

Sp = array(0,c(K,K,nf,nt))
for (t in 1:nt) { 
    tindx <- times[t]
    for (f in 1:nf){
        A = diag(K)
        for (p in 1:P){
            A = A - m[p,,,tindx]*exp(-1i*2*pi*freqs[f]*p)
        }
        A = solve(A)
        Sp[,,f,t] = A%*%S[,,tindx]%*%t(Conj(A))
    }
}
return(Sp)
}


tvar_spec_Guo_P = function(m, s, times, freqs) {
nt <- length(times)
nf <- length(freqs)
P <- dim(m)[4]
K <- dim(m)[2]
if (length(dim(s))==3){
    S = s
}else{
    S = array(s,c(dim(s)[1:2],nt))
}

Sp = array(0,c(K,K,nf,nt))
for (t in 1:nt) { 
    tindx <- times[t]
    for (f in 1:nf){
        A = diag(K)
        for (p in 1:P){
            A = A - m[tindx,,,p]*exp(-1i*2*pi*freqs[f]*p)
        }
        A = solve(A)
        Sp[,,f,t] = A%*%S[,,tindx]%*%t(Conj(A))
    }
}
return(Sp)
}