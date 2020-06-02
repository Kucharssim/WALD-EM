dwald <- function(x, alpha, nu, log = FALSE) {
  
  lpdf <- log(alpha) - 1.0/2.0 * log(2*pi) - 3.0/2.0*log(x) - (alpha - nu*x)^2/(2*x)
  
  if(log) return(lpdf) else return(exp(lpdf))
}

simulatePath <- function(alpha, nu, freq = 1000){
  evid <- 0
  counter <- 1
  while(evid[counter] < alpha){
    evid <- c(evid, evid[counter] + rnorm(1, nu/freq, sqrt(1/freq)))
    counter <- counter + 1
  }
  evid[counter] <- alpha
  
  return(evid)
}

alpha <- 2
nu <- 5
freq <- 1000
set.seed(2020)
paths <- replicate(500, simulatePath(alpha, nu, freq), simplify = FALSE)
y <- lengths(paths)/freq
tmax <- max(lengths(paths))/freq
ymin <- min(sapply(paths, min))

png(filename = here::here("figures", "wald_distribution.png"), width = 512, height = 512)
par(mfrow = c(2, 1),
    oma = c(5,4,1,1) + 0.1,
    mar = c(0,0,0,0))

hist(y, breaks = 50, freq = FALSE, axes = FALSE, xlab = "", ylab = "", xlim = c(0, tmax), main = "")
rug(y)
curve(dwald(x = x, alpha, nu), from = 0, to = tmax, add = TRUE, lwd = 2)
plot(0, xlim = c(0, tmax), type = "n", ylim = c(ymin, alpha), bty = "n", axes = FALSE, xlab = "", ylab = "")
lapply(paths, function(path) lines(1:length(path)/freq, path, col = adjustcolor("black", alpha = 0.025)))

lines(1:length(paths[[2]])/freq, paths[[2]], lwd = 1.5)
arrows(x0=0, y0=0, x1=alpha/nu, y1 = alpha, lwd = 1.5)
segments(0, 0, alpha/nu, 0, lty = 3)
segments(alpha/nu, 0, alpha/nu, alpha, lty = 3)
text(alpha/nu + 0.1, alpha / 2, expression(nu == alpha/t), cex = 2)
text(alpha/(2*nu), ymin/2, "t", cex = 1.5)
abline(h = alpha, lty = 2)
axis(2, at = c(0, alpha), labels = c(0, expression(alpha)), las = 1, cex.axis = 1.5)
axis(1, at = pretty(c(0, tmax)), cex.axis = 1.25)
mtext("Time (sec)", side = 1, line = 3, cex = 1.75)
dev.off()
