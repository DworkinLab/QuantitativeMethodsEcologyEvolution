# Why t is not a good measure of effect.


n <- seq(3, 500, by = 1)


t_calc <- function(n, diff = 2) {
  
  x <- rnorm(n, mean = 10, sd = 3)
  y <- rnorm(n, mean = (10 - diff), sd =3)
  t.test(x, y, alternative = "two.sided")$statistic
}

tvals_n <- sapply(n, t_calc)

plot(tvals_n ~ n,
     pch = 20, cex = 1.2, col = "purple",
     ylab = "t statistic",
     xlab = "sample size in each group")


d_calc <- function(n, diff = 2) {
  
  x <- rnorm(n, mean = 10, sd = 3)
  y <- rnorm(n, mean = (10 - diff), sd =3)
  diff <- mean(x) - mean(y)
  std_dev_pooled <- sqrt( ( ((n-1)*var(x)) + ((n-1)*var(y)) )   / (2*n - 2) )
  d <- diff/std_dev_pooled
  return(d)
}

dvals_n <- sapply(n, d_calc)


plot(dvals_n ~ n,
     pch = 20, cex = 1.2, col = "purple",
     ylab = "Cohen's  d",
     xlab = "sample size in each group")


d_t_mean_calc <- function(n, diff = 2) {
  
  x <- rnorm(n, mean = 10, sd = 3)
  y <- rnorm(n, mean = (10 - diff), sd =3)
  difference <- mean(x) - mean(y)
  std_dev_pooled <- sqrt( ( ((n-1)*var(x)) + ((n-1)*var(x)) )   / (2*n - 2) )
  d <- difference/std_dev_pooled
  tstat <- t.test(x, y, alternative = "two.sided")$statistic
  mean_x <- mean(x)
  mean_y <- mean(y)
  return(c(mean_x = mean_x, mean_y = mean_y, 
           difference = difference,
           d = d, tstat = tstat))
}


d_t_mean_n <- t(sapply(n, d_t_mean_calc))

par(mfrow = c(2,2))

plot(mean_x ~ n, data = d_t_mean_n,
     pch = 20, cex = 1.2, col = "red",
     ylim = c(7, 12),
     ylab = "mean of x (red), y (blue)",
     xlab = "sample size in each group")


points(y = d_t_mean_n[,2], x= n,
     pch = 20, cex = 1.2, col = "blue",
     ylim = c(6, 14),
     ylab = "mean of y",
     xlab = "sample size in each group")

plot(difference ~ n, data = d_t_mean_n,
     pch = 20, cex = 1.2, col = "orange",
     ylab = "difference between means",
     xlab = "sample size in each group")

plot(d ~ n, data = d_t_mean_n,
     pch = 20, cex = 1.2, col = "purple",
     ylab = "Cohen's d",
     xlab = "sample size in each group")


plot(tstat.t ~ n, data = d_t_mean_n,
     pch = 20, cex = 1.2, col = "grey",
     ylab = "t statistic",
     xlab = "sample size in each group")

par(mfrow = c(1,1))




## The rubber ruler problem




d_rubber_calc <- function(n, diff = 2) {
  
  x <- rnorm(n, mean = 10, sd = 3)
  y <- rnorm(n, mean = (10 - diff), sd =3)
  difference <- mean(x) - mean(y)

  std_dev_pooled <- sqrt( ( ((n-1)*var(x)) + ((n-1)*var(x)) )   / (2*n - 2) )
  d <- difference/std_dev_pooled
  
  d_no_rubber <- difference/3
  
  d_rubber <- diff/std_dev_pooled
  
  tstat <- t.test(x, y, alternative = "two.sided")$statistic
  mean_x <- mean(x)
  mean_y <- mean(y)
  return(c(mean_x = mean_x, mean_y = mean_y, 
           estimated_difference = difference, 
           actual_difference = diff,
           d_no_rubber = d_no_rubber,
           d_rubber = d_rubber,
           d = d))}


rubber_time <- t(sapply(n, d_rubber_calc))

par(mfrow = c(2,2))

plot(mean_x ~ n, data = rubber_time,
     pch = 20, cex = 1.2, col = "red",
     ylim = c(7, 12),
     ylab = "mean of x (red), y (blue)",
     xlab = "sample size in each group")

points(y = rubber_time[,2], x= n,
     pch = 20, cex = 1.2, col = "blue",
     ylim = c(6, 14),
     ylab = "mean of y",
     xlab = "sample size in each group")

plot(estimated_difference ~ n, data = rubber_time,
     pch = 20, cex = 1.2, col = "orange",
     ylab = "difference between means",
     xlab = "sample size in each group")

plot(d ~ n, data = rubber_time,
     pch = 20, cex = 1.2, col = "purple",
     ylab = "Cohen's d",
     xlab = "sample size in each group")



par(mfrow = c(1,1))



