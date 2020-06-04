library(here)

x          <- seq(0, 100, by = 0.1)
lambda     <- 0.1*dnorm(x, 10, sd = 4) + 0.2*dnorm(x, 30, sd = 8) + 0.3*dnorm(x, 50, sd=15) + 0.3*dnorm(x, 80, sd = 5)
lambda     <- lambda/sum(lambda)

png(filename = here::here("figures", "attention_window.png"), width = 1000, height = 600, units = "px", pointsize = 25)
par(mfrow = c(2,3), mar = c(2.5, 1, 2.5, 1))
fix_x      <- 55
att_window <- exp(-(x-fix_x)^2/(2*220))
plot(x, lambda,            type = "l", bty = "n", main = expression(lambda[x](x)),             xlab = "", ylab = "", yaxt = "n", lwd = 3)
plot(x, att_window,        type = "l", bty = "n", main = expression(a[x](x*"|"*s[x]==55)),     xlab = "", ylab = "", lwd = 3)
abline(v = fix_x, lty = 2)
plot(x, lambda*att_window, type = "l", bty = "n", main = expression(omega[x](x*"|"*s[x]==55)), xlab = "", ylab = "", yaxt = "n", lwd = 3)
abline(v = fix_x, lty = 2)

fix_x      <- 20
att_window <- exp(-(x-fix_x)^2/(2*220))
plot(x, lambda,            type = "l", bty = "n", main = expression(lambda[x](x)), xlab = "", ylab = "", yaxt = "n", lwd = 3)
plot(x, att_window,        type = "l", bty = "n", main = expression(a[x](x*"|"*s[x]==20)),     xlab = "", ylab = "", lwd = 3)
abline(v = fix_x, lty = 2)
plot(x, lambda*att_window, type = "l", bty = "n", main = expression(omega[x](x*"|"*s[x]==20)), xlab = "", ylab = "", yaxt = "n", lwd = 3)
abline(v = fix_x, lty = 2)
dev.off()