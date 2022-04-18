# bootstrap CI - pretty

### Data in

dll.data = read.csv("http://beaconcourse.pbworks.com/f/dll.csv", 
    header = TRUE, stringsAsFactors = TRUE)   #data frame input




dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
dll.data$genotype <- relevel(dll.data$genotype, "wt")  # Setting the wild-type (wt) as reference
dll.data <- na.omit(dll.data)  # I am only doing this to look at some diagnostic plots, which get all bothered by the missing data.
dll.data$tarsus.scaled <- as.numeric(scale(dll.data$tarsus))

################


regression.1 <- lm(SCT ~ tarsus.scaled, data = dll.data)

par(mfrow=c(1,1))
plot(jitter(SCT, factor = 0.75) ~ tarsus.scaled, data = dll.data, 
    xlim = range(dll.data$tarsus.scaled), 
    pch = 20, col = "red")
abline(regression.1, lwd = 2)


# bootstrap for pretty CI curves
# values to predict on

new_tarsus <- data.frame(tarsus.scaled = seq(from = 
    range(dll.data$tarsus.scaled)[1],
        to = range(dll.data$tarsus.scaled)[2], 
    length = 1000))   
    # Makes a new data frame for the x values to predict on

nitt = 1000

predicted.boot <- matrix(NA, nrow = nitt, ncol = 1000)

for(i in 1:nitt) {
	dat.boot <- dll.data[sample(nrow(dll.data), nrow(dll.data), replace = T),]
	model.boot <- lm(SCT ~ tarsus.scaled, data = dat.boot)
	predicted.boot[i,] <- predict(model.boot, new_tarsus, interval = "none")
}



# quick check that we get the lines out. We add these onto the plot
for(i in 1:25) {
	lines(predicted.boot[i,] ~ new_tarsus$tarsus.scaled, col="grey")
}


# but what we want is the quantiles
lower <- apply(X = predicted.boot, MARGIN = 2, FUN = quantile, probs = c(0.025))
upper <- apply(X = predicted.boot, MARGIN = 2, FUN = quantile, probs = c(0.975))


par(mfrow=c(1,1))
plot(jitter(SCT, factor = 0.75) ~ tarsus.scaled, data = dll.data, 
    xlim = range(dll.data$tarsus.scaled), 
    pch = 20, col = "purple",
    ylab = "Sex Comb Teeth", xlab = " Tarsus (Centered Values)")
abline(regression.1, lwd = 2)

# plot the lines from the quantiles of the bootstrap.
lines(x = new_tarsus[,1], y = lower, lwd = 4, lty = 2, col = "blue")
lines(x = new_tarsus[,1], y = upper, lwd = 4, lty = 2, col = "blue")


# let's add on classic  parametric CI's
predicted.model.1  <- predict(regression.1, new_tarsus, interval="confidence")
#lines(x = new_tarsus[,1], y = predicted.model.1[,1], lwd=4) # this would be the same as abline(model.1)
lines(x = new_tarsus[,1], y = predicted.model.1[,3], lwd = 3, lty = 3, col = "red")
lines(x = new_tarsus[,1], y = predicted.model.1[,2], lwd = 3, lty = 3, col = "red")


# One nice thing to do is draw the polygon for the Confidence band, and use some of R's pretty features to make it transparent.
par(mfrow=c(1,1))
plot(jitter(SCT, factor = 0.75) ~ tarsus.scaled, data = dll.data, 
    xlim = range(dll.data$tarsus.scaled), 
    pch = 20, col = "purple",
    ylab = "Sex Comb Teeth", xlab = "Tarsus (Centered Values)")
abline(regression.1, lwd = 2)


# The trick with polygon is remember that for the X values you have to in one direction along the line (for the upper CI), and then go in the reverse order for the lower CI (along the X )
polygon( x = c(new_tarsus[,1], rev(new_tarsus[,1])),  
         # for the x values we first go forward, and then reverse along the vector
         y = c(upper, rev(lower)), 
         # We go in the "forward" direction for the Upper CI, and the reverse direction for lower CI
         col = "blue", border = "blue" )


# The problem with this is you can not see the points. You could just re-plot the points using points. But we can add some transparency


# First we see what colour red is represented as in rgb
col2rgb("red") # 0 0 255

# Now we convert that to hexmode 
as.hexmode(255) ## This is ff

# This means that for the blue  & green channel we use 0 (in fact two zero's each for hexbin)

# The final number in hexmode is for degree of transparency. So if we want it ~ 20% transparent
as.hexmode(80)  # 50
  # So the final colour we want is col="#ff000050"
  
# can also just use rgb(1, 0, 0, 0.2) 
 
# Replot  
par(mfrow=c(1,1))
plot(jitter(SCT) ~ tarsus.scaled, data = dll.data, 
  xlim = range(dll.data$tarsus.scaled), ylab = "# SCT", xlab = "centered tarsus length",
  col = densCols(y = SCT, x = tarsus.scaled), pch = 16)
lines(x = new_tarsus[,1], y = predicted.model.1[,1], lwd = 3, col = "red") # this would be the same as abline(model.1), but with controlled x and y limits


# The trick with polygon is remember that for the X values you have to in one direction along the line (for the upper CI), and then go in the reverse order for the lower CI (along the X )

polygon( x = c(new_tarsus[,1], rev(new_tarsus[,1])),  # for the x values first forward, then reverse along vector
         y = c(upper, rev(lower)), # We go in "forward" direction for Upper CI, reverse direction for lower CI
         col = "#ff000050", border = NA)
  