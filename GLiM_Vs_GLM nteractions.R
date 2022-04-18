# updated march 28th, 2022 - Ian Dworkin


# just because it is count data, does not mean that a poisson or negative binomial give the "best" fit.

require(bbmle)
require(MASS)
require(car)
require(emmeans)
require(ggplot2)
require(glmmTMB)
library(dotwhisker)


dll.data = read.csv("https://raw.githubusercontent.com/DworkinLab/DworkinLab.github.io/master/dataSets/Dworkin2005_ED/dll.csv", 
                    header = TRUE, stringsAsFactors = TRUE)

dll.data = na.omit(dll.data)

dll.data$temp <- factor(dll.data$temp)
dll.data$replicate <- factor(dll.data$replicate)
str(dll.data) 
summary(dll.data)

dll.data$genotype <- relevel(dll.data$genotype, ref = "wt")
dll.data$tarsus_Z <- scale(dll.data$tarsus, center = TRUE, scale = T)
str(dll.data)



ggplot(dll.data, aes(x = SCT, col = genotype:temp)) +
  geom_histogram(binwidth = 1)



#If you remember from our earlier discussions of probability distributions, it appeared as if the normal distribution was a better fit for the # of Sex Comb teeth over either poisson or negative binomial, despite the fact that it represents count data constrained 

fit.normal <- glm(SCT ~ 1 + (genotype + temp +tarsus_Z)^2, 
                  family = "gaussian",
                  data = dll.data)

fit.poisson <- glm(SCT ~ 1 + (genotype + temp +tarsus_Z)^2, 
                   family = "poisson",
                   data = dll.data)

fit.nb0 <- glm.nb(SCT ~ 1 + (genotype + temp +tarsus_Z)^2,
                       data = dll.data)

# the same

fit.nb1 <- glmmTMB(SCT ~ 1 + (genotype + temp +tarsus_Z)^2,
                   family = nbinom1(link = "log"),
                   data = dll.data)

fit.nb2 <- glmmTMB(SCT ~ 1 + (genotype + temp +tarsus_Z)^2,
                   family = nbinom2(link = "log"),
                   data = dll.data)


dwplot(fit.normal) + geom_vline(xintercept = 0, lty = 2)

dwplot(list(fit.poisson, fit.nb0, fit.nb1)) + 
  geom_vline(xintercept = 0, lty = 2)

AICctab(fit.normal, fit.poisson, fit.nb0, fit.nb1, fit.nb2, nobs = 1918, base = T)

# but does this mean we should be using the normal distribution to model the residual variation?

# Let's look at the full models we are interested in, including poisson, NB and Gaussian.

glmm_fit.poisson <- glmmTMB(SCT ~ 1 + (genotype + temp +tarsus_Z)^2 
                        + rr(0 + genotype + temp | line/replicate),
                        family = poisson(link = "log"),
                        data = dll.data)

glmm_fit.nb2 <- glmmTMB(SCT ~ 1 + (genotype + temp +tarsus_Z)^2 
                        + rr(0 + genotype + temp | line/replicate),
                   family = nbinom2(link = "log"),
                   control = glmmTMBControl(optimizer=optim,
                                          optArgs=list(method="BFGS")),
                   data = dll.data)

glmm_fit.nb3 <- glmmTMB(SCT ~ 1 + (genotype + temp + tarsus_Z)^2 
                        + rr(0 + genotype + temp | line)
                        + (1|line:replicate),
                        family = nbinom2(link = "log"),
                        data = dll.data)

glmm_fit.gs <- glmmTMB(SCT ~ 1 + (genotype + temp +tarsus_Z)^2 
                       + rr(0 + genotype + temp | line/replicate),
                       data = dll.data)

glmm_fit.gs2 <- glmmTMB(SCT ~ 1 + (genotype + temp +tarsus_Z)^2 
                        + rr(0 + genotype + temp | line)
                        + (1|line:replicate),
                        data = dll.data)


AICtab(glmm_fit.poisson, glmm_fit.nb2, glmm_fit.nb3,
       glmm_fit.gs, glmm_fit.gs2, 
       weights = T, base = TRUE )

# normal distribution still vastly outperforms

summary(glmm_fit.gs)
summary(glmm_fit.gs2)



# Again you can compare the estimated variance from the data (or residual variance from model) and the estimates from nb or Gaussian to see why.
