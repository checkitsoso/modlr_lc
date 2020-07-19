x.weight = function (x) {
  n = length(x)
  if(n == 1) return(list(x=x, w=1))
  x = sort(x)
  i = which(diff(x) > 0)
  i.n = c(i, n)
  list(x=x[i.n], w=i.n-c(0, i))
}

## Log-likelihood function
logLik.lcd = function (object, x, w=NULL, ...) {
  if(is.null(w)) { xw = x.weight(x); x = xw$x; w = xw$w }
  sum(dlcd(x, object, log=TRUE) * w)
}

dlcd = function (x, lcd, log=FALSE) {
  logd = rep(-Inf, length(x))
  j = x >= lcd$lower & x <= lcd$upper
  xj = x[j]
  knots = c(lcd$lower, lcd$theta)
  jk = indx(xj, knots)
  logd[j] = lcd$coef[1,jk] + lcd$coef[2,jk] * xj - log(lcd$C)
  if (log) logd else exp(logd)
}

qlcd = function (p, x, lcd) {
  if (p == 0) {
    q = lcd$lower
  } else if (p == 1) {
    q = lcd$upper
  } else {
    n = length(x)
    phi = log(dlcd(x, lcd))
    Fhat = plcd(x, lcd)
    xj = max(x[Fhat <= p])
    j = length(x[x <= xj])
    a = (p - Fhat[j])/(Fhat[j+1] - Fhat[j])
    b = (x[j+1] - x[j]) * (phi[j+1] - phi[j])
    q = xj + (x[j+1] - x[j]) * qloglin(a, b)
  }
  return(q)
}

rlcd = function (n, x, lcd) {
  x = sort(x)
  u = runif(n)
  x = sapply(u, qlcd, x=x, lcd=lcd)
  return(x)
}

maxima.gradient = function (lcd, x, w, tol = -Inf) {
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  nk = length(knots)
  index = match(knots, x)
  grad1 = gradient.lcd(lcd, x, w, "x")
  grad = round(grad1,6)
  ii = numeric(0)
  for (i in 1:(nk-1)) {
    grad.between = grad[index[i]:index[i+1]]
    if (max(grad.between) > 0) {
      ima = which.max(grad.between)
      ii = c(ii, index[i] + ima - 1)
    }
  }
  theta = x[ii]
  g = grad[ii]
  gmax = ifelse(length(g) > 0, max(g), -Inf)
  j = g > tol & ! theta %in% knots
  list(theta=theta[j], gradient=g[j], gmax=gmax)
}

gradient.lcd = function (lcd, x, w=NULL, theta) {
  if (length(theta) < 1) return(numeric(0))
  if (is.null(w)) { xw = x.weight(x); x = xw$x; w = xw$w }
  knots = c(lcd$lower, lcd$theta)
  nk = length(knots)
  xxt = attr(x, "xx")
  if (is.null(xxt)) xxt = rev(cumsum(rev(w))) * x -  rev(cumsum(rev(x * w)))
  if (!is.numeric(theta) && theta == "knots") {
    xxt = xxt[indx(knots, x)]
    px = cbind(0, lcd$cpk[1:2,-nk,drop=FALSE])
    xxt + sum(w) * (lcd$cpk[2,nk] - px[2,] - (1 - px[1,]) * knots)
  } else {
    if (is.numeric(theta)) xxt = xxt[indx(theta, x)] else theta = x
    cpx = cpx(lcd, theta, order=1)
    return(xxt + sum(w)*(lcd$cpk[2,nk] - cpx[2,] - (1-cpx[1,])*theta))
  }
}

cpx = function (lcd, x, order=0) {
  x[x < lcd$lower] = lcd$lower
  x[x > lcd$upper] = lcd$upper
  knots = c(lcd$lower, lcd$theta)
  nk = length(knots)
  k = indx(x, knots)
  a = knots[k]
  dpx = matrix(0, nrow=order+1, ncol=length(x))
  rownames(dpx) = paste0("x", 0:order)
  coef = lcd$coef[,k,drop=FALSE]
  fkk = lcd$fk[k]
  if(any(j <- a != x & coef[2,] == 0)) {      # horizontal segments
    dpx[1,j] = fkk[j] * (x[j] - a[j])
    if(order >= 1) dpx[2,j] = dpx[1,j] * (x[j] + a[j]) * .5
    if(order >= 2) dpx[3,j] = dpx[1,j] * (x[j]^2 + x[j] * a[j] + a[j]^2) / 3
  }
  if(any(j2 <- a != x & coef[2,] != 0)) {     # non-horizontal segments
    x1 = fkk[j2] / coef[2,j2]
    y1 = exp(coef[1,j2] + coef[2,j2] * x[j2] - log(lcd$C)) / coef[2,j2]
    dpx[1,j2] = y1 - x1
    if(order >= 1) dpx[2,j2] = x[j2] * y1 - a[j2] * x1 - dpx[1,j2] / coef[2,j2]
    if(order >= 2)
      dpx[3,j2] = x[j2]^2 * y1 - a[j2]^2 * x1 - 2 * dpx[2,j2] / coef[2,j2]
  }
  if(nk > 1) cbind(0, lcd$cpk[1:(order+1),,drop=FALSE])[,k,drop=FALSE] + dpx
  else dpx 
}

line.lcd = function (lcd0, lcd1, x, w, ll0, tol=1e-10) {
  llt = function(alpha) {
    m = new.lcd((1 - alpha) * lcd0$alpha + alpha * lcd1$alpha, lcd1$theta,
                (1 - alpha) * lcd0$pi + alpha * lcd1$pi, lcd0$lower, lcd0$upper)
    ll = logLik.lcd(m, x, w)
    list(lcd=m, ll=ll)
  }
  grad = gradient.lcd(lcd0, x, w, "knots")
  grad[1] = - grad[1]
  delta = sum(grad * (c(lcd1$alpha, lcd1$pi) - c(lcd0$alpha, lcd0$pi))) * .333
  convergence = 0
  alpha = 1
  repeat{
    new = llt(alpha)
    if(new$ll >= ll0 + alpha * delta) break 
    if(alpha < tol) {convergence=1; new=list(lcd=lcd0, ll=ll0); break}
    # alpha = alpha * .5
    alpha = alpha * .2
  }
  list(lcd=new$lcd, ll=new$ll, alpha=alpha, convergence=convergence)
}

new.lcd = function (alpha, theta = NULL, pi = NULL, lower, upper) {
  if (length(theta) > 0) {
    k = length(theta)
    if (length(pi) < k) 
      pi = rep(pi, len = k)
    o = order(theta)
    theta = theta[o]
    pi = pi[o]
  }
  else theta = pi = numeric(0)
  c0 = -alpha * lower + c(0, cumsum(pi * theta))
  c1 = alpha - c(0, cumsum(pi))
  c1[abs(c1) <= 1e-6] = 0
  knots1 = c(lower, theta)
  knots2 = c(theta, upper)
  dk = knots2 - knots1
  nk = length(dk)
  fk1 = exp(c0 + c1 * knots1)
  fk2 = c(fk1[-1], exp(c0[nk] + c1[nk] * upper))
  dpk = matrix(0, nrow = 3, ncol = nk)
  rownames(dpk) = paste0("x", 0:2)
  if (any(j <- c1 == 0)) {
    dpk[1, j] = fk1[j] * dk[j]
    dpk[2, j] = dpk[1, j] * (knots2[j] + knots1[j]) * 0.5
    dpk[3, j] = dpk[1, j] * (knots2[j]^2 + knots2[j] * knots1[j] + knots1[j]^2)/3
  }
  if (any(j2 <- !j)) {
    x1 = fk1[j2]/c1[j2]
    y1 = fk2[j2]/c1[j2]
    dpk[1, j2] = y1 - x1
    dpk[2, j2] = knots2[j2] * y1 - knots1[j2] * x1 - dpk[1, j2]/c1[j2]
    dpk[3, j2] = knots2[j2]^2 * y1 - knots1[j2]^2 * x1 - 2 * dpk[2, j2]/c1[j2]
  }
  C = sum(dpk[1, ])
  fk = fk1/C
  cpk = dpk = dpk/C
  for (i in 1:nrow(dpk)) cpk[i, ] = cumsum(cpk[i, ])
  structure(list(alpha = alpha, C = C, theta = theta, pi = pi, lower = lower, upper = upper, 
                 coef = rbind(c0, c1), fk = fk, dpk = dpk, cpk = cpk), class = "lcd")
}

simplify.lcd = function (lcd) {
  if (any(j0 <- lcd$pi == 0)) {
    nk = length(lcd$theta) + 1
    j = which(!j0)
    pi = lcd$pi[j]
    theta = lcd$theta[j]
    j1 = c(1,j+1)
    dpk = cpk = lcd$cpk[, c(j,nk), drop=FALSE]
    for(i in 1:nrow(dpk)) dpk[i,] = c(cpk[i,1], diff(cpk[i,]))
    lcd = structure(list(alpha=lcd$alpha, C=lcd$C, theta=theta, pi=pi,
                         lower=lcd$lower, upper=lcd$upper,
                         coef=lcd$coef[,j1,drop=FALSE],
                         fk=lcd$fk[j1], dpk=dpk, cpk=cpk),
                    class = "lcd")
  }
  lcd
}

qloglin = function (u, t) {
  if (abs(t) > 1e-06) {
    z = log(1 + ((exp(t) - 1) * u))/t
  } else {
    z = u + t * u * (1 - u)/2
  }
  return(z)
}
