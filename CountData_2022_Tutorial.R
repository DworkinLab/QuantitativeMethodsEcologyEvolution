# Count data
library(MASS) # for glm.nb
library(bbmle)
library(emdbook) # for the dzin function
library(glmmTMB) # many options for both nb and zero inflation
library(ggplot2)
library(ggbeeswarm)


ggplot(quine, aes(y = Days, x = Sex, col = Sex)) + geom_quasirandom()

ggplot(quine, aes(x = Days, col = Sex)) + geom_dotplot(stackgroups = TRUE, method = "histodot", binwidth = 1.0, binpositions = "all")  +
  scale_y_continuous(NULL, breaks = NULL)

hist(quine$Days, breaks = max(quine$Days))

mean(quine$Days)
median(quine$Days)


quine.pois1 <- glm(Days ~ Sex, 
                   data = quine, 
                   family = "poisson")


summary(quine.pois1)

-2*logLik(quine.pois1)
plot(quine.pois1)

# Should we worry about over-dispersion?

# For the likelihood, we would write it like

#design matrix

X <- model.matrix(~ Sex, data = quine)

quine.pois.NLL <- function(a, b, y = quine$Days) {
	 det.model <- exp(a + b*X[,2])
	 -sum(dpois(y, lambda = det.model, log = T)) }

quine.pois.MLE <- mle2(quine.pois.NLL, 
    start = list(a = 2, b = 0.16))

summary(quine.pois.MLE)

prof.quine.pois <- profile(quine.pois.MLE)

plot(prof.quine.pois)
confint(prof.quine.pois)



# Should we worry about over-dispersion?
quine.pois0 <- glm(Days ~ 1, 
                   data = quine, 
                   family = "poisson")

quine.pois0 # deviance is 10X the df!

var(quine$Days)
exp(coef(quine.pois0))

sim_poisson <- rpois(146, exp(coef(quine.pois0)))
mean(sim_poisson)
var(sim_poisson)

hist(sim_poisson, breaks = max(sim_poisson))


# If we wanted to use a quasi-likelihood, we could do this with glm
quine.QuassiPois1 <- glm(Days ~ Sex, data = quine, family="quasipoisson")
summary(quine.QuassiPois1)


# But we will not really discuss this here....



# Instead let us utilize a Negative binomial errors approach, to allow for over-dispersion

# Here we assume that y is distributed according to a poisson, where lambda itself varies according to a gamma distribution. Using the pre-built function in MASS, glm.nb

quine.nb1 <- glm.nb(Days ~ Sex, data = quine)
summary(quine.nb1)

#note deviance

#car::Anova(quine.nb1)

quine.nb1_v2 <- glmmTMB(Days ~ Sex, 
                        data = quine, 
                        family = nbinom2(link = "log"))

summary(quine.nb1_v2)

fixef(quine.nb1_v2)

# Fitting the negative binomial model as a full MLE

# we can use dnbinom() , using the "mu"  and "size" (over-dispersion) parameterization (see Bolker page 124 )

quine.nbinom.NLL <- function(a, b, dispersion, y = quine$Days){
	 det.model <- exp(a + b*X[,2])
	 -sum(dnbinom(y, mu = det.model, size = dispersion, log = T)) }

#Confirm we get get the same deviance as glm.nb
-2*quine.nbinom.NLL(a=2.7229, b=0.1649, dispersion=1.0741)


# Now we fit the model using mle2, and compare to glm.nb
quine.nb.MLE <- mle2(quine.nbinom.NLL , start=list(a=2, b=0.16, dispersion=1))

summary(quine.nb.MLE)

summary(quine.nb1)

coef(quine.nb1)

prof.quine.nb <- profile(quine.nb.MLE)

plot(prof.quine.nb)
confint(prof.quine.nb)



## Linear parameterization of the nbinom is a bit different (glmmTMB)
quine.nb2_v2 <- glmmTMB(Days ~ Sex, 
                        data = quine, 
                        family = nbinom1(link = "log"))

quine.nb2_v2

############
# While it will just repeat what we have just done, we can also do this by using the poisson and gamma distributions (page 124 in Bolker).
# Again, we can get a negative binomial as a result of poisson sampling, where lambda varies according to a gamma distribution
# lambda is ~gamma with shape  = dispersion parameter (k or size), and mean = mu
# your random variable (in this example quine$Days) is then poisson distributed with mean lambda
################



# Zero inflated distributions (binomial, poisson, neg. binom)
# One other way of considering this data is as a mix of a two distributions, one all zero's, and one with some positive mean, but including some zeroes. IN other words the negative binomial distribution is inflated with zeroes. 


# Do we have more zeroes then we would expect given either the poisson or nb?

hist(quine$Days, breaks = max(quine$Days))

# we could do a simulation, here without any "structural" zeroes, so all zeroes come from the negative binomial
sample1 <- rzinbinom(146, mu = exp(2.76), size = 1.07412, zprob = 0)
hist(sample1, breaks = max(sample1))


# Just to get a sense of what it means to have more zero inflation, 0.1 probability an observation is a "structural" zero.

sample1 <- rzinbinom(146, mu = exp(2.76), size = 1.07412, zprob = 0.1)
hist(sample1, breaks = max(sample1))


#glmmTMB. We can first just ask whether there is any evidence for zero inflation, before even worrying about how our predictors may influence zi.

quine.zinb1_v1 <- glmmTMB(Days ~ Sex, 
                        ziformula = ~ 1, 
                        data = quine, 
                        family = nbinom2(link = "log"))

quine.zinb1_v1
summary(quine.zinb1_v1)

# the zero inflation is on the logit scale, so to interpret it more easily
plogis(fixef(quine.zinb1_v1)$zi)


# We can deal with this using a zero inflated distribution (page 285)
# Ben Bolker has already done the work for us (in the emdbook package)

emdbook::dzinbinom

# So let us modify our support function to allow for this. We will add a parameter zprop_a for the zero inflation.

quine.zinbinom.NLL <- function(a, b, dispersion, zprob_a, y = quine$Days){
	 det.model <- exp(a + b*X[,2])
	 zif <- zprob_a
	 -sum(dzinbinom(y, mu=det.model, size = dispersion, zprob = zif, log = T)) 
}

# First we double check that we get the same result as the negative binomial when zprob=0
-2*quine.zinbinom.NLL(a=2.7229, b=0.1649, dispersion=1.0741, zprob_a = 0)


# Now we can fit the model. 

# Since I do not know what a good starting point for zprob would be I compute the proportion of zeroes

with(quine, length(Days[Days==0])/length(Days))

# And use this as a starting value

# I will also use the MLE from the pure negative binomial model as starting values for the other parameters

quine.zinb.MLE <- mle2(quine.zinbinom.NLL, 
                       start=list(a = 2.72, b = 0.1649, dispersion = 1.07, zprob_a = 0.06))

summary(quine.zinb.MLE)
coef(quine.zinb.MLE)


# from glmmTMB
fixef(quine.zinb1_v1)

# to compare dzinbinom results transform from probability to logit scale for the zinf terms, so let's transform 
qlogis(coef(quine.zinb.MLE)["zprob_a"]) # takes probability, returns logit

plogis(fixef(quine.zinb1_v1)$zi) # takes logit, returns probability

prof.quine.zinb <- profile(quine.zinb.MLE)
plot(prof.quine.zinb)
confint(prof.quine.zinb)


# How do these models compare to each other?
AICctab(quine.zinb.MLE,quine.nb.MLE, quine.pois.MLE, 
        nobs=length(quine$Days), weights=T, base=T )


# Is the differences in AIC just a function of the additional parameters (i.e no change in the deviance)
deviance(quine.zinb.MLE)
deviance(quine.nb.MLE)

# Of course you may want to examine predicted vs observed to get a sense of whether there is a good global fit!

# There is also an R library which has many useful functions called pscl that includes various zero inflated distributions. I have not used it, but it looks very promising. 
# Also see ZIGP, gamlss, and flexmix. MCMCglmm and glmmTMB for mixed models.


# we could also compare to Gaussian

Gaussian_MLE <- function(a,b, sigma, y=quine$Days){
	 det.model <-a + b*X[,2]
	 -sum(dnorm(y, mean=det.model, sd=sigma, log=T)) }



normal_MLE_fit <- mle2(Gaussian_MLE, 
                       start=list(a = mean(quine$Days), b = 0, sigma = sd(quine$Days)))

normal_MLE_fit

AICctab(quine.zinb.MLE,quine.nb.MLE, quine.pois.MLE, normal_MLE_fit,
        nobs = length(quine$Days), weights = T, base = T )


