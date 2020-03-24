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