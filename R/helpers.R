matrix2df <- function(m, origin = c("bottom-left", "top-left")){
  origin <- match.arg(origin)
  nr <- nrow(m)
  nc <- ncol(m)
  
  n  <- nr*nc
  
  out <- data.frame(x = rep(seq_len(nr), nc),
                    y = rep(seq_len(nc), each = nr),
                    s = NA)
  
  for(i in 1:n){
    x_ind <- out$x[i]
    y_ind <- switch (origin,
      "bottom-left" = out$y[i],
      "top-left"    = nc + 1 - out$y[i]
    )
    
    out$s[i] <- m[x_ind, y_ind]
  }
  
  return(out)
}

rwald <- function(n=1, alpha, nu){
  mu <- alpha/nu
  lambda <- alpha^2
  zeta <- rnorm(n)
  zeta_sq <- zeta^2
  x <- mu + (mu^2*zeta_sq)/(2*lambda) - mu/(2*lambda)*sqrt(4*mu*lambda*zeta_sq + mu^2*zeta_sq^2)
  z <- runif(n)
  
  which_x <- z <= mu / (mu+x)
  
  out <- rep(NA, n)
  out[ which_x] <- x[which_x]
  out[!which_x] <- mu^2/x[!which_x]
  
  return(out)
}

rtnorm <- function(n=1, mean=0, sd=1, lb = 0, ub = Inf){
  pl <- pnorm(lb, mean, sd)
  pu <- pnorm(ub, mean, sd)
  pp  <- runif(n, pl, pu)
  
  qnorm(pp, mean, sd)
}
