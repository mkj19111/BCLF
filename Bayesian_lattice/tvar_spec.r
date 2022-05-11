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

tvar_spec = function(m, s, times, freqs) {
nt <- length(times)
p <- nrow(m)
N <- ncol(m)
omega <- 2 * pi * freqs
spec <- matrix(0, nrow = length(omega), ncol = nt)
eom <- exp(-omega * 1i) # complex number

for (tt in 1:nt) { 
 tindx <- times[tt]
 tempx <- c(-m[p:1,tindx], 1)
 sp <- abs(polyval(tempx,eom)) 
 spec[,tt] <- s[tindx] / sp^2
}
return(spec)
}

############ subfunctions ####################
# Polynomial evaluation similar to "polyval" of Matlab.
# The input argument coeff is a vector of length n+1 whose 
# elements are the coefficients in descending powers of 
# the polynomial to be evaluated; that is,
# y = coeff_1 * x_n + coeff_2 * x_n–1 + … + coeff_n * x + coeff_n+1

polyval <- function(coeff,x) {
   n <- length(coeff)
	 y <- rep(coeff[1],length(x))
   for (ii in 2:n) {
     y <- coeff[ii] + x * y
   }
   return(y)
 }
