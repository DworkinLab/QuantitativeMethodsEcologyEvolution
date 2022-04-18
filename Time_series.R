
# How to deal with temporal or spatial correlation (time series data). That is lack of independence among observations (after controlling for predictors)
#EXAMPLE FROM ZUUR ET AL. 2009

library(nlme)

Hawaii <- read.table("./Hawaii.txt", h=T)
#data(Hawaii)
#Abundace of three species of birds collected from 1956-2003.
Hawaii <- na.omit(Hawaii)

Hawaii$Birds <- sqrt(Hawaii$Moorhen.Kauai) 
# for variance stabilization. Of course we could just use a different variance structure, but let's keep this simple

plot(Hawaii$Year, Hawaii$Birds, 
     xlab = "Year", ylab = "sqrt(abundance)",type = "b",
     main = "abundance of Kauai Moorhen",
     pch = 20, col = "red")


# What if we wanted to understand the relationship between annual rainfall and abundance, how do we take into account the increasing abundance over time?

M.0 <- lm(Birds ~ Year + Rainfall,  
          data = Hawaii)

summary(M.0)

abline(M.0, lty = 2, col = rgb(0, 0, 1, 0.5))

# Rainfall seems to have no effect but is that only because year so large an effect? We have not really accounted for the fact that the abundance in one year could be correlated to the abundance in the next year.

par(mfrow=c(2,1))
plot(resid(M.0) ~ Hawaii$Year, type = "b", 
     pch = 20, col = "grey",
     ylab = "model residuals", xlab = "Year")

# We can see this
acf(resid(M.0), main = "autocorrelation of residuals")

# Clearly we have violated the assumption of independence. More importantly we are not accounting for a biologically important source of information!!!!


# Compound symmetry structure in R using gls()

M.1 <- gls(Birds ~ Rainfall + Year, 
           data = Hawaii, 
           correlation = corCompSymm( form= ~Year))

# corAR1. This is the one that we likely want:
M.2 <- gls(Birds ~ Rainfall + Year, 
           data = Hawaii, 
           correlation = corAR1( form = ~ Year))

summary(M.1)
summary(M.2)


# Note below AR1 is a relatively better fit. Compound symmetry is worse than ignoring lack of independence!  

bbmle::AICctab(M.0, M.1, M.2, 
               nobs = nrow(Hawaii), 
               weights = TRUE, base = TRUE)


# For more sense of the various correlation structures
?corClasses
