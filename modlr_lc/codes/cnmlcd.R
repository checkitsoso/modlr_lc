cnmlcd = function (x, lcd, maxit = 100, tol = 1e-06) {
  library(lsei)
  lambda = 1e-15
  n = length(x)
  xw = x.weight(x)
  x = xw$x
  w = xw$w
  nx = length(x)
  lower = x[1]
  upper = x[nx]
  attr(x, "xx") = rev(cumsum(rev(w))) * x - rev(cumsum(rev(x * w)))
  if (missing(lcd)) lcd = new.lcd(alpha = 0, lower = lower, upper = upper)
  ll = logLik.lcd(lcd, x, w)
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    lcd.old = lcd
    ll.old = ll
    g = maxima.gradient(lcd, x, w = w)
    if (length(g$theta) != 0) {
      nsp = g$theta
      nsl = length(nsp)
      if (nsl >= 1) lcd = new.lcd(lcd$alpha, c(lcd$theta, nsp), c(lcd$pi, double(nsl)), lcd$lower, lcd$upper)
    }
    knots = c(lcd$lower, lcd$theta)
    nk = length(knots)
    cpkr = lcd$cpk[, nk] - cbind(0, lcd$cpk[, -nk, drop = FALSE])
    mu = cpkr[2, ] - cpkr[1, ] * knots
    grad = n * mu + attr(x, "xx")[indx(knots, x)]
    mm = cpkr[3, ] - (knots + rep(knots, rep.int(nk, nk))) * 
      cpkr[2, ] + tcrossprod(knots) * cpkr[1, ]
    mm[upper.tri(mm)] = 0
    mm = mm + t(mm)
    diag(mm) = diag(mm)/2
    H = mm - tcrossprod(mu)
    e = eigen(H)
    v2 = sqrt(e$values[e$values >= e$values[1] * lambda])
    kr = length(v2)
    R = t(e$vectors[, 1:kr, drop = FALSE]) * v2
    p = grad/n + drop(c(-lcd$alpha, lcd$pi) %*% H)
    b = drop(crossprod(e$vectors[, 1:kr, drop = FALSE], p))/v2
    r1 = pnnls(R, b, 1)
    lcd1 = lcd
    lcd1$alpha = -r1$x[1]
    lcd1$pi = r1$x[-1]
    r = line.lcd(lcd, lcd1, x, w = w, ll0 = ll.old)
    lcd = r$lcd
    ll = r$ll
    if (any(lcd$pi == 0)) lcd = simplify.lcd(lcd)
  }
  list(lcd = lcd, ll = ll, num.iterations = i, max.gradient = g$gmax, convergence = convergence)
}

cnmlcd.mode = function (x, w, m=0, lcd, maxit=100, tol=1e-10) {
  lambda = 1e-15
  n = length(x)
  if (missing(w)) {
    xw = x.weight(x)
    x = xw$x
    w = xw$w
  }
  if (!(m %in% x)) {
    x = c(x, m)
    w = c(w, 0)
    o = order(x)
    x = x[o]
    w = w[o]
  }
  nx = length(x)
  lower = x[1]
  upper = x[nx]
  if (nx <= 2) {
    lcd = new.lcd(alpha=0, lower=lower, upper=upper)
    return(list(lcd = lcd, ll = logLik.lcd(lcd, x, w), num.iterations = 0, max.gradient = 0, convergence = 0))
  }
  attr(x, "xx") = rev(cumsum(rev(w))) * x - rev(cumsum(rev(x * w)))
  if (missing(lcd)) {
    if (m %in% c(lower, upper)) {
      lcd = new.lcd(alpha=0, lower=lower, upper=upper)
    } else {
      lcd = new.lcd(alpha=0, theta=0, pi=0, lower=lower, upper=upper)
    }
  }
  ll = logLik.lcd(lcd, x, w)
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    lcd.old = lcd
    ll.old = ll
    g = maxima.gradient(lcd, x, w = w)
    if (length(g$theta) != 0) {
      nsp = g$theta
      nsl = length(nsp)
      if (nsl >= 1) {
        lcd = new.lcd(lcd$alpha, c(lcd$theta, nsp), c(lcd$pi, double(nsl)), lcd$lower, lcd$upper)
      }
    }
    knots = c(lcd$lower, lcd$theta)
    nk = length(knots)
    cpkr = lcd$cpk[, nk] - cbind(0, lcd$cpk[, -nk, drop = FALSE])
    mu = cpkr[2, ] - cpkr[1, ] * knots
    grad = n * mu + attr(x, "xx")[lsei::indx(knots, x)]
    mm = cpkr[3, ] - (knots + rep(knots, rep.int(nk, nk))) * cpkr[2, ] + tcrossprod(knots) * cpkr[1, ]
    mm[upper.tri(mm)] = 0
    mm = mm + t(mm)
    diag(mm) = diag(mm)/2
    H = mm - tcrossprod(mu)
    e = eigen(H)
    v2 = sqrt(e$values[e$values >= e$values[1] * lambda])
    kr = length(v2)
    R = t(e$vectors[, 1:kr, drop = FALSE]) * v2
    p = grad/n + drop(c(-lcd$alpha, lcd$pi) %*% H)
    b = drop(crossprod(e$vectors[, 1:kr, drop = FALSE], p))/v2
    g1 = cbind(0, diag(rep(nk-1)))
    knot_m = sum(knots < m)
    if (sum(knots %in% m) == 0) {
      g2 = matrix(0, nrow=1, ncol=nk)
      g2[1, 1:(knot_m-1)] = -1
    } else {
      g2 = matrix(0, nrow=2, ncol=nk)
      g2[1, knots < m] = -1
      g2[2, 1:(knot_m +1)] = 1
    }
    G = rbind(g1, g2)
    H = rep(0, nrow(G))
    lcd1 = lcd
    flag = TRUE
    r1 = tryCatch(limSolve::lsei(A=R, B=b, G=G, H=H, type=2), error = function(e) flag <<- FALSE)
    if (!flag) {r2 = lsei::lsei(a=R, b=b, e=G, f=H)} else {r2 =NULL}
    if (flag) {lcd1$alpha = -r1$X[1];lcd1$pi = r1$X[-1]} else {lcd1$alpha = -r2[1];lcd1$pi = r2[-1]}
    
    r = line.lcd(lcd, lcd1, x, w = w, ll0 = ll.old)
    lcd = r$lcd
    ll = r$ll
    if (any(lcd$pi == 0)) 
      lcd = simplify.lcd(lcd)
  }
  # if (! m %in% c(lower, lcd$theta, upper))
  #   lcd = new.lcd(lcd$alpha, theta = c(lcd$theta, m), pi = c(lcd$pi, 0), lower = lcd$lower, upper = lcd$upper)
  return(list(lcd = lcd, ll = ll, num.iterations = i, max.gradient = g$gmax, convergence = convergence))
}

cnmlcd.symm = function (x, w, m=0, lcd, maxit=100, tol=1e-10) {
  lambda = 1e-15
  n = length(x)
  x = abs(x)
  if (missing(w)) {
    xw = x.weight(x)
    x = xw$x
    w = xw$w
  }
  if (!(m %in% x)) {
    x = c(x, m)
    w = c(w, 0)
    o = order(x)
    x = x[o]
    w = w[o]
  }
  nx = length(x)
  lower = x[1]
  upper = x[nx]
  if (nx <= 2) {
    lcd = new.lcd(alpha=0, lower=lower, upper=upper)
    return(list(lcd = lcd, ll = logLik.lcd(lcd, x, w), num.iterations = 0, max.gradient = 0, convergence = 0))
  }
  attr(x, "xx") = rev(cumsum(rev(w))) * x - rev(cumsum(rev(x * w)))
  if (missing(lcd)) {
    if (m %in% c(lower, upper)) {
      lcd = new.lcd(alpha=0, lower=lower, upper=upper)
    } else {
      lcd = new.lcd(alpha=0, theta=0, pi=0, lower=lower, upper=upper)
    }
  }
  ll = logLik.lcd(lcd, x, w)
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    lcd.old = lcd
    ll.old = ll
    g = maxima.gradient(lcd, x, w = w)
    if (length(g$theta) != 0) {
      nsp = g$theta
      nsl = length(nsp)
      if (nsl >= 1) {
        lcd = new.lcd(lcd$alpha, c(lcd$theta, nsp), c(lcd$pi, double(nsl)), lcd$lower, lcd$upper)
      }
    }
    knots = c(lcd$lower, lcd$theta)
    nk = length(knots)
    cpkr = lcd$cpk[, nk] - cbind(0, lcd$cpk[, -nk, drop = FALSE])
    mu = cpkr[2, ] - cpkr[1, ] * knots
    grad = n * mu + attr(x, "xx")[lsei::indx(knots, x)]
    mm = cpkr[3, ] - (knots + rep(knots, rep.int(nk, nk))) * cpkr[2, ] + tcrossprod(knots) * cpkr[1, ]
    mm[upper.tri(mm)] = 0
    mm = mm + t(mm)
    diag(mm) = diag(mm)/2
    H = mm - tcrossprod(mu)
    e = eigen(H)
    v2 = sqrt(e$values[e$values >= e$values[1] * lambda])
    kr = length(v2)
    R = t(e$vectors[, 1:kr, drop = FALSE]) * v2
    p = grad/n + drop(c(-lcd$alpha, lcd$pi) %*% H)
    b = drop(crossprod(e$vectors[, 1:kr, drop = FALSE], p))/v2
    
    ########
    # r1 = pnnls(R, b, 1)
    C1 = cbind(0, diag(rep(nk-1)))
    C2 = matrix(0, nrow=2, ncol=nk)
    C2[1, knots < m] = -1
    C2[2, knots <= m] = 1
    C = rbind(C1, C2)
    D = rep(0, nrow(C))
    r1 = limSolve::lsei(A=R, B=b, G=C, H=D, type=2)
    ########
    
    lcd1 = lcd
    lcd1$alpha = -r1$X[1]
    lcd1$pi = r1$X[-1]
    r = line.lcd(lcd, lcd1, x, w = w, ll0 = ll.old)
    lcd = r$lcd
    ll = r$ll
    if (any(lcd$pi == 0)) 
      lcd = simplify.lcd(lcd)
  }
  if (! m %in% c(lower, lcd$theta, upper))
    lcd = new.lcd(lcd$alpha, theta = c(lcd$theta, m), pi = c(lcd$pi, 0), lower = lcd$lower, upper = lcd$upper)
  
  # symmetric constraint added
  lcd$alpha = sum(lcd$pi)
  lcd$theta = sort(c(-lcd$theta,lcd$theta))
  lcd$pi    = c(rev(lcd$pi),lcd$pi)
  lcd$C     = 2*lcd$C
  lcd$lower = -lcd$upper
  c0 = c(rev(lcd$coef[1,-1]),lcd$coef[1,])
  c1 = c(-rev(lcd$coef[2,-1]),lcd$coef[2,])
  lcd$coef = matrix(c(c0,c1),nrow=2,ncol=2*dim(lcd$coef)[2]-1,byrow=T,dimnames=dimnames(lcd$coef))
  ll = ll - n*log(2)
  
  return(list(lcd = lcd, ll = ll, num.iterations = i, max.gradient = g$gmax, convergence = convergence))
}
