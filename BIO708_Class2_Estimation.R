## ----setup, include=FALSE---------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(show.signif.stars=F)
options(digits=3)


## ----librariesToUse---------------------------------------------------------------------------
library(ggplot2)


## ---------------------------------------------------------------------------------------------
sct_data  <- read.csv("http://beaconcourse.pbworks.com/f/dll.csv",
                       h = T, stringsAsFactors = TRUE)


## ---------------------------------------------------------------------------------------------
sct_data <- na.omit(sct_data) 


## ---- echo=F----------------------------------------------------------------------------------
ggplot(sct_data, aes(y = SCT, x = tarsus)) +
     geom_jitter(width = 0, height = 0.3, alpha = 0.25)


## ---- echo=F----------------------------------------------------------------------------------
ggplot(sct_data, aes(y = SCT, x = tarsus, color = genotype)) +
        geom_jitter(width = 0, height = 0.25, alpha = 0.2) +
        geom_smooth( method = loess)


## ---- echo=F----------------------------------------------------------------------------------
ggplot(sct_data, aes(y = SCT, x = tarsus, color = genotype)) +
        geom_jitter(width = 0, height = 0.25, alpha = 0.2) +
        geom_smooth( method = lm)


## ---- echo=F----------------------------------------------------------------------------------
ggplot(sct_data, aes(y = SCT, x = tarsus, color = genotype)) +
        geom_jitter(width = 0, height = 0.25, alpha = 0.2) +
        geom_smooth( method = lm)


## ---- echo = F, include = F-------------------------------------------------------------------
x <- seq(1, 10, by = 0.4)

y <- rnorm(length(x), 
           mean = 6 + 2*x - 1.5*(x^2) + 0.3*(x^3),  sd = 2*x)


## ---------------------------------------------------------------------------------------------
plot(y ~ x, pch = 20,
    ylab = "response variable", xlab = "predictor variable", 
    main = "What is the underlying process generating this data?")


## ---------------------------------------------------------------------------------------------
plot(y ~ x, pch = 20,
    ylab = "response variable", xlab = "predictor variable", 
    main = "What is the underlying process generating this data?")

model_1 <- lm(y ~ poly(x,3))
param <- coef(model_1)
lines(x = x, y = fitted(model_1), col = "purple", lwd=2)


## ---------------------------------------------------------------------------------------------
plot(y ~ x, pch = 20,
    ylab = "response variable", xlab = "predictor variable", 
    main = "What is the underlying process generating this data?")

model_2 <- lm( y ~ x)
lines(x=x, y = fitted(model_2), col="purple", lwd = 2, lty = 2)
text(x=8,y=0, "Parametric models", col = "black")

model_3 <- lm( y ~ poly(x,2))
lines(x = x, y = fitted(model_3), col="blue", lwd = 2)

model_4 <- lm(y~ poly(x, 3))
lines(x = x, y = fitted(model_4), col="red", lwd=2, lty=2)

model_22 <- lm(y~ poly(x, 21))
lines(x = x, y = fitted(model_22), col="grey", lwd = 2)

legend("topleft", bty = "n",
       col = c("purple", "blue", "red", "grey"), lty = 1, 
       legend = c("linear", "quadratic", "cubic", "23rd order polynomial"))


## ---------------------------------------------------------------------------------------------
plot(y ~ x, pch = 20,
    ylab = "response variable", xlab = "predictor variable", 
    main = "a polynomial of order 3 (cubic)")

lines(x=x, y = fitted(model_2), col="purple", lwd = 2, lty = 2)

lines(x = x, y = fitted(model_3), col="blue", lwd = 2)

lines(x = x, y = fitted(model_4), col="red", lwd=2, lty=2)

lines(x = x, y = fitted(model_22), col="grey", lwd = 2)

text(x=4.25, y=(max(y) -5), cex = 1,
    expression(paste("~N( ", mu, " = ", 6 + 2*x -1.5*x^2 +0.3*x^3, ", ",sigma, " = 2x)")) )


## ---------------------------------------------------------------------------------------------
y2 <- rnorm(length(x), 
            mean = 6 + 2*x - 1.5*(x^2) + 0.3*(x^3),  
            sd=2*x)

plot(y2 ~ x, pch = 20, 
     ylab = "response", xlab = "predictor",
     main="What happens to the model fit with new data (but same process)")

lines(x=x, y = fitted(model_2), col="purple", lwd=2, lty=2)
lines(x=x, y = fitted(model_3), col="blue", lwd =2, lty=1)
lines(x=x, y = fitted(model_4), col="red", lwd=2)
lines(x=x, y = fitted(model_22), col="grey", lwd=2)


## ---------------------------------------------------------------------------------------------
y3 <- rnorm(length(x), 
            mean = 6 + 2*x,  
            sd = 1.5)

plot(y3 ~ x, 
     ylab = "response", xlab = "predictor", pch = 20, cex = 2,
     main = "how do we get 'best' estimates of the slope and intercept?")



## ---------------------------------------------------------------------------------------------
plot(y3 ~ x, 
     ylab = "response", xlab = "predictor", pch = 20, cex = 2,
     main = "how do we get 'best' estimates of the slope and intercept?")

abline(a = 16.7, b = 0, col = "grey", lwd = 1.2, lty = 3)
abline(a = 10, b = 1, col = "blue", lwd = 0.8, lty = 2)
abline(a = 3, b = 3, col = "red", lwd = 0.8, lty = 2)
abline(lm(y3 ~ x), col = "purple", lwd = 1, lty = 1)


## ---------------------------------------------------------------------------------------------
plot(y3 ~ x, 
     ylab = "response", xlab = "predictor", pch = 20, cex = 2,
     main = "how do we get 'best' estimates of the slope and intercept?")

abline(a = 16.7, b = 0, col = "grey", lwd = 1.2, lty = 3)
abline(a = 10, b = 1, col = "blue", lwd = 0.8, lty = 2)
abline(a = 3, b = 3, col = "red", lwd = 0.8, lty = 2)
abline(lm(y3 ~ x), col = "purple", lwd = 1, lty = 1)


## ---- echo =T---------------------------------------------------------------------------------
ResSumSq <- function(slope) { sum((y -slope*x )^2) } 


## ---- echo = TRUE-----------------------------------------------------------------------------
x <- seq(from = 0, to = 10, by = 0.25)


## ---- echo = TRUE-----------------------------------------------------------------------------
y <- rnorm(x, mean = 0 + 1.25*x, sd = 2)


## ---------------------------------------------------------------------------------------------
plot(y ~ x, pch = 20, col = "blue")


## ---- echo = TRUE-----------------------------------------------------------------------------
mod1 <- lm(y ~ 0 + x)
summary(mod1)


## ---------------------------------------------------------------------------------------------
plot(y ~ x, pch = 20)
abline(a = 0, b = 1.25, col = "red")
abline(mod1, col = "blue")


## ---- echo = T--------------------------------------------------------------------------------
sum(resid(lm(y ~ 0 + x))^2)


## ---- echo = T--------------------------------------------------------------------------------
ResSumSq <- function(slope) { sum((y - slope*x )^2) } 


## ---- echo = T--------------------------------------------------------------------------------
ResSumSq(slope = 0)


## ---- echo = T--------------------------------------------------------------------------------
ResSumSq(slope = 0.25)


## ---- echo = T--------------------------------------------------------------------------------
ResSumSq(slope = - 1.25)


## ---- echo = T--------------------------------------------------------------------------------
ResSumSq(slope = 1)


## ---- echo=TRUE-------------------------------------------------------------------------------
slope <- seq(-1.5, 3, by = 0.01)

head(slope)


## ---- echo = T--------------------------------------------------------------------------------
rss_est <- sapply(slope, ResSumSq )

head(rss_est)


## ---- echo = TRUE-----------------------------------------------------------------------------
together <- cbind(slope, rss_est)


## ---------------------------------------------------------------------------------------------
plot(rss_est ~ slope, type = "p", pch = 20, cex = 0.5,
     ylim = c(0, max(rss_est)),
     ylab = "Sum of Squares", xlab = "Candidate slope")


## ---- echo = TRUE-----------------------------------------------------------------------------
brute_force_estimate <- together[which.min(together[,2]),] 

brute_force_estimate


## ---------------------------------------------------------------------------------------------
plot(rss_est ~ slope, type = "p", pch = 20, cex = 0.5,
     ylim = c(0, max(rss_est)),
     ylab = "Sum of Squares", xlab = "Candidate slope")

segments(x0 = brute_force_estimate[1], x1 = brute_force_estimate[1],
         y0 = max(rss_est), y1 = brute_force_estimate[2], 
         lwd = 2, col = "purple", lty = 3)

segments(x0 = min(slope), x1 = brute_force_estimate[1],
         y0 = brute_force_estimate[2], y1 = brute_force_estimate[2], 
         lwd = 2, col = "blue", lty = 2)


## ---- echo = T--------------------------------------------------------------------------------
optim_min <- optim(par = 3, ResSumSq , method="BFGS")

optim_min$par
optim_min$value


## ---- echo =T---------------------------------------------------------------------------------
sum(resid(lm(y ~ 0 + x))^2)


## ---------------------------------------------------------------------------------------------
coef(mod1)


## ---- echo =TRUE------------------------------------------------------------------------------
dll_data <- read.csv("http://beaconcourse.pbworks.com/w/file/fetch/35183279/dll.csv", header=TRUE)   #data frame input

dll_data <- na.omit(dll_data) # removing missing values


## ---- echo =TRUE------------------------------------------------------------------------------
dll_lm_1 <- lm(SCT ~ 1 + tarsus, data=dll_data)


## ---------------------------------------------------------------------------------------------
coef(dll_lm_1)


## ---------------------------------------------------------------------------------------------
sum(resid(dll_lm_1)^2)


## ---- echo = T--------------------------------------------------------------------------------
y <- dll_data$SCT
x <- dll_data$tarsus


## ---- echo = TRUE-----------------------------------------------------------------------------
MinSumSq2 <- function(b) { 
	intercept <- b[1] # First value in b, corresponding to the intercept
	slope     <- b[2] # second value in b, corresponding to the slope of the model 
	sum((y -intercept -slope*x)^2)  } 


## ---- echo = TRUE-----------------------------------------------------------------------------
optim_min_2 <- optim(par=c(0, 0), 
                    MinSumSq2, method="BFGS")


## ---- echo=TRUE-------------------------------------------------------------------------------
optim_min_2$par
optim_min_2$value


## ---- echo = TRUE-----------------------------------------------------------------------------
MinSumSq2(b = optim_min_2$par)


## ---- echo = T--------------------------------------------------------------------------------
LeastAbsDev <- function(b) { 
	intercept <- b[1] 
	slope     <- b[2] 
	sum(abs(y -intercept -slope*x))  } 


## ---- echo = T--------------------------------------------------------------------------------
optim_min_3 <- optim(par=c(0.1, 0.1), 
                     LeastAbsDev, method="BFGS")

optim_min_2$par # LSE
optim_min_3$par # LAD


## ---- echo = TRUE-----------------------------------------------------------------------------

fakeRidge <- function(b, lambda = 0.2) { 
	intercept <- b[1] 
	slope     <- b[2] 
	(sum((y -intercept -slope*x)^2) + lambda * sum(b^2))}

fakeLasso <- function(b, lambda = 0.2) { 
	intercept <- b[1] 
	slope     <- b[2] 
	(sum((y -intercept -slope*x)^2) + lambda * sum(abs(b)))}


## ---- echo = TRUE-----------------------------------------------------------------------------
optim_min_4 <- optim(par=c(0.1, 0.1), 
                     fakeRidge, method="BFGS")

optim_min_5 <- optim(par=c(0.1, 0.1), 
                     fakeLasso, method="BFGS")
optim_min_2$par #LSE
optim_min_4$par # Ridge (L2)
optim_min_5$par # lasso (L1)

