# Demonstrating that the focus on significant results and/or "file drawer syndrome" causes publication bias in estimates


n <- sort(rep(3:250, times = 10)) # 10 simulations at each sample size


# function to calculate all of the relevant quantities

d_t_mean_p_calc <- function(n, diff = 1.4) {
  
  x <- rnorm(n, mean = 10, sd = 4)
  
  y <- rnorm(n, mean = (10 - diff), sd = 4)
  
  difference <- mean(x) - mean(y)
  
  std_dev_pooled <- sqrt( ( ((n-1)*var(x)) + ((n-1)*var(y)) )   / (2*n - 2) )
  
  d <- (difference/std_dev_pooled) * (1 - (3/(4*(2*n-2) - 1)) ) # corrected for small sample size ala Hedges and Olkin 1985
  
  tt <- t.test(x, y, alternative = "two.sided")
  tstat <- tt$statistic
  pval  <- tt$p.value
  
  mean_x <- mean(x)
  mean_y <- mean(y)
  
  return(c(mean_x = mean_x, 
           mean_y = mean_y, 
           difference = difference,
           d = d, 
           tstat = tstat,
           pval = pval,
           n = n))}
           
           


d_t_mean_p_n <- t(sapply(n, d_t_mean_p_calc)) # generate all of the values under simulation


d_t_mean_p_n <- as.data.frame(d_t_mean_p_n)

par(mfrow = c(2,2))

plot(mean_x ~ n, data = d_t_mean_p_n,
     pch = 20, cex = 1, col = rgb(1,0, 0, alpha = 0.4),
     ylim = c(7, 12),
     ylab = "mean of x (red), y (blue)",
     xlab = "sample size in each group")
abline(h = 10, lwd = 2, lty = 2)


points(y = d_t_mean_p_n$mean_y, x = n,
     pch = 20, cex = 1, col = rgb(0,0,1, alpha = 0.4),
     ylim = c(6, 14),
     ylab = "mean of y",
     xlab = "sample size in each group")
abline(h = 8.6, lwd = 2, lty = 2)

plot(difference ~ n, data = d_t_mean_p_n,
     pch = 20, cex = 1, col = rgb(1, 0.5, 0, alpha = 0.5),
     ylim = c(-4, 6),
     ylab = "difference between means",
     xlab = "sample size in each group")
abline(h = 1.4, lwd = 2, lty = 2)
abline(h = 0, lwd = 2, lty = 3, col = "grey")

plot(d ~ n, data = d_t_mean_p_n,
     pch = 20, cex = 1, col = rgb(0.5,0,1, alpha = 0.5),
     ylim = c(-2, 2),
     ylab = "Cohen's d",
     xlab = "sample size in each group")
abline(h = (1.4/4), lwd = 2, lty = 2)
abline(h = 0, lwd = 2, lty = 3, col = "grey")


plot(tstat.t ~ n, data = d_t_mean_p_n,
     pch = 20, cex = 1, col = "grey", 
     ylim = c(-2, 6),
     ylab = "t statistic",
     xlab = "sample size in each group")
abline(h = 0, lwd = 2, lty = 3, col = "grey")

par(mfrow = c(1,1))



## How about if we exclude non-significant values? i.e simulating a field of discipline where only "significant" effects get published

dim(d_t_mean_p_n)

d_t_mean_significantOnly <- d_t_mean_p_n[d_t_mean_p_n$pval <= 0.05,]


dim(d_t_mean_significantOnly)

par(mfrow = c(2,2))

plot(mean_x ~ n, data = d_t_mean_significantOnly,
     pch = 20, cex = 1, col = "red",
     ylim = c(7, 12),
     ylab = "mean of x (red), y (blue)",
     xlab = "sample size in each group")
abline(h = 10, lwd = 2, lty = 2)

points(y = d_t_mean_significantOnly$mean_y, x = n,
     pch = 20, cex = 1.2, col = "blue",
     ylim = c(6, 14),
     ylab = "mean of y",
     xlab = "sample size in each group")
abline(h = 8.6, lwd = 2, lty = 2)


plot(difference ~ n_, data = d_t_mean_significantOnly,
     pch = 20, cex = 1.2, col = "orange",
     ylim = c(-4, 6),
     ylab = "difference between means",
     xlab = "sample size in each group")
abline(h = 1.4, lwd = 2, lty = 2)
abline(h = 0, lwd = 2, lty = 3, col = "grey")     

plot(d ~ n, data = d_t_mean_significantOnly,
     pch = 20, cex = 1.2, col = "purple",
     ylim = c(-2, 2),
     ylab = "Cohen's d",
     xlab = "sample size in each group")
abline(h = (1.4/4), lwd = 2, lty = 2)
abline(h = 0, lwd = 2, lty = 3, col = "grey")

plot(tstat.t ~ n, data = d_t_mean_significantOnly,
     pch = 20, cex = 1.2, col = "grey",
     ylim = c(-2, 6),
     ylab = "t statistic",
     xlab = "sample size in each group")
abline(h = 0, lwd = 2, lty = 3, col = "grey")

par(mfrow = c(1,1))




## how much bias?

## For sample sizes less than 50 
mean(d_t_mean_significantOnly[n < 50, 4])

mean(d_t_mean_p_n[n < 50, 4])

## the value of d from this is approximately 2 X greater than it should be!

##From 30 to 100, still 50% greater than what it should be for only the significant tests.
mean(d_t_mean_significantOnly[n >= 30 & n <= 100 , 4])

mean(d_t_mean_p_n[n >= 30 & n <= 100 , 4])


##From 101 to 250, still about 10% greater than what it should be for only the significant tests.
mean(d_t_mean_significantOnly[n > 101 & n <= 250 , 4])

mean(d_t_mean_p_n[n > 100 & n <= 250 , 4])

# Even still the estimates are inflated!



## This is a useful way of visualizing the bias (suggested by Ben Bolker)
trans_black <- adjustcolor("black", alpha.f = 0.2)
trans_red <- adjustcolor("red", alpha.f = 0.2)

plot(d_t_mean_p_n$n, abs(d_t_mean_p_n[,"d"]), 
    pch = 16, col = trans_black,
    ylab = "|d|", xlab = "sample size",
    main = "selective publication biases estimates of effects")

lines(lowess(d_t_mean_p_n$n, 
    abs(d_t_mean_p_n[,"d"])), lwd = 2 )

points(d_t_mean_significantOnly$n, 
    abs(d_t_mean_significantOnly[,"d"]),
       col = trans_red,)
       
lines(lowess(d_t_mean_significantOnly$n, 
    abs(d_t_mean_significantOnly[,"d"])), 
    lwd = 2, col = "red")
      
      
legend("topright", 
    legend = c("all", " p < 0.05"),
    pch = c(20, 1), bty = "n", cex = 1.5,
    col = c(trans_black, trans_red))

