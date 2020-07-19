#### Modal Linear Regression using Log-concave Distributions
lcmlr_aa = function (y, z, beta=NULL, tol=1e-10, maxit=100) {
  z = as.matrix(z)
  n = length(y)
  p = ifelse(is.null(nrow(z)), 1, ncol(z))
  # Initial regression parameter
  if (is.null(beta)) beta = lm(y ~ 0+z)$coefficients
  # Initial log-concave error density
  if (p == 1) e = y - z*beta else e = y - z%*%beta
  res = cnmlcd.mode(e)
  lcd = res$lcd
  ll = res$ll
  
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    beta.old = beta
    lcd.old = lcd
    ll.old = ll
    
    # Step 1: Update regression parameter estimate
    q = length(lcd$theta) + 2
    knot = c(lcd$lower, lcd$theta, lcd$upper)
    phi = dlcd(knot, lcd, log=TRUE)
    slopes = diff(phi) / diff(knot)
    intercepts = -slopes * knot[-q] + phi[-q]
    U = cbind(-kronecker(z, slopes), kronecker(diag(rep(-1, n)), rep(1, q-1)))
    v = -(rep(slopes, n)*rep(y, each=q-1)) - rep(intercepts, n)
    sparse.mat = rbind(cbind(z, Matrix(0, nrow=n, ncol=n, sparse=TRUE)), 
                       cbind(-z, Matrix(0, nrow=n, ncol=n, sparse=TRUE)))
    LP = Rglpk_solve_LP(obj=c(v, y - max(knot), -y + min(knot)),
                        mat=t(rbind(U, sparse.mat)),
                        dir=c(rep("==", n+p)),
                        rhs=c(rep(0, p), rep(-1, n)),
                        bounds=list(lower=list(ind=1:(n*(q-1)), val=rep(0, n*(q-1)))),
                        max=TRUE)
    beta = LP$auxiliary$dual[1:p]
    
    # Step 2: Update log-concave error density estimate
    if (p == 1) e = y - z*beta else e = y - z%*%beta
    res = cnmlcd.mode(e)
    lcd = res$lcd
    ll = res$ll
  }
  list(beta=beta, lcd=lcd, ll=ll, num.iterations=i, convergence=convergence)
}

lcmlr_aa_sg = function (y, z, beta=NULL, tol=1e-10, maxit=100) {
  z = as.matrix(z)
  n = length(y)
  p = ifelse(is.null(nrow(z)), 1, ncol(z))
  # Initial regression parameter
  if (is.null(beta)) beta = lm(y ~ 0+z)$coefficients
  # Initial log-concave error density
  if (p == 1) e = y - z*beta else e = y - z%*%beta
  res = cnmlcd.mode(e)
  lcd = res$lcd
  ll = res$ll
  
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    beta.old = beta
    lcd.old = lcd
    ll.old = ll
    
    # Step 1: Update regression parameter estimate
    q = length(lcd$theta) + 2
    knot = c(lcd$lower, lcd$theta, lcd$upper)
    phi = dlcd(knot, lcd, log=TRUE)
    slopes = diff(phi) / diff(knot)
    intercepts = -slopes * knot[-q] + phi[-q]
    delta = -z %o% slopes
    gamma = y %o% slopes + matrix(1, nrow=n, ncol=1) %o% intercepts
    ll.temp = -Inf
    for (l in 1:100) {
      if (ll <= ll.temp + tol) break
      beta.temp = beta
      ll.temp = ll
      llk = matrix(0, nrow=n, ncol=q-1)
      for (k in 1:(q-1)) llk[, k] = delta[,,k]%*%beta + gamma[,,k]
      idx = apply(llk, 1, which.min)
      subg = apply(-matrix(slopes[idx], nrow=n, ncol=p)*z, 2, sum)
      for (j in 1:100) {
        beta = beta + subg
        for (k in 1:(q-1)) llk[, k] = delta[,,k]%*%beta + gamma[,,k]
        ll = sum(apply(llk, 1, min))
        if (ll <= ll.temp) {
          subg = subg / 2
          beta = beta.temp
          ll = ll.temp
        } else break
      }
    }
    # Step 2: Update log-concave error density estimate
    if (p == 1) e = y - z*beta else e = y - z%*%beta
    res = cnmlcd.mode(e)
    lcd = res$lcd
    ll = res$ll
  }
  list(beta=beta, lcd=lcd, ll=ll, num.iterations=i, convergence=convergence)
}