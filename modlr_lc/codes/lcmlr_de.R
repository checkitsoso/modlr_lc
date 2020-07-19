# Log-concave Modal Linear Regression (LCMLR) using Differential Evoultion (DE)
lcmlr_de = function (y, z, L, U, N, C=0.9, tau=0.8, kmax=100, cl) {
  n = length(y)
  p = ifelse(is.null(nrow(z)), 1, ncol(z))
  if (missing(N)) N = max(50, 20*p)
  
  # Log-likelihood function
  llk = function(beta, y, z) {
    if (p == 1) e = y - z*beta else e = y - z%*%beta
    return(-cnmlcd.mode(e)$ll)
  }
  
  ## DE: Control ##
  if (missing(cl)) {
    control = DEoptim.control(strategy=2, trace=TRUE, 
                              NP=N, CR=C, F=tau, itermax=kmax, packages=c("lsei"))
  } else {
    control = DEoptim.control(strategy=2, trace=TRUE, 
                              NP=N, CR=C, F=tau, itermax=kmax, parallelType=1, cluster=cl, 
                              packages=c("lsei"), parVar=c("y", "z")) 
  }
  
  ## DE: Execution ##
  res = DEoptim(fn=llk, lower=L, upper=U, control=control, y, z)
  return(list(beta=res$optim$bestmem, ll=-res$optim$bestval))
}
