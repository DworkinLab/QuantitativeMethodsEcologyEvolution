# Written by Ian Dworkin - last modified on March 7th, 2022

# Chunk 1: libraries and read in data as usual
# Chunk 2: Overview of computation for resampling.
# Chunk 3: Fitting the basic model with lm() before attempting the resampling
# Chunk 4: A basic permutation test using a for() loop
# Chunk 5: A basic permutation test using an explicit functional call, and replicate() 
# Chunk 6: Distribution of p values under permutation
# Chunk 7: using RRPP for more complex models
# Chunk 8: Some notes on the computation of permuted data sets.



# Chunk 1: Read in, and set up data as usual.

library(RRPP)

dll.data = read.csv("https://raw.githubusercontent.com/DworkinLab/DworkinLab.github.io/master/dataSets/Dworkin2005_ED/dll.csv", header=TRUE)
dll.data <- na.omit(dll.data)  
dll.data$temp <- as.factor(dll.data$temp)
dll.data$replicate <- as.factor(dll.data$replicate)


### Chunk 2: Overview of computation for resampling  

 # For performing the randomization/permutation test (or bootstrap) there are two important pieces of code to know. 1) the sample() function and 2) how to perform an iterative loop (with a for loop or replicate).
 
 #  the function that performs the resampling is sample()
 #  sample(x, size, replace=False)
 #  where x is the vector to be resampled from, size is the length of the vector to be generated, and
 #  replace is whether or not sampling with replacement is to take place
 
 # the default value of size is length(x) which is what we want for both permutations and bootstrapping.

   # If all of this code gives you a headache, just follow below, and remember that all you need to do is change N for number of resampling events, and change your model and test statistic.
########



####### Chunk 3: Examining the basic model we are fitting using lm (prior to getting to the resampling). Note since, genotype only has two levels, this is just a t-test.

dll.anova = lm(SCT ~ genotype, data=dll.data) 
summary(dll.anova) # look at model summary statistics
 
dll.anova.fstat  = summary(dll.anova)$fstat # just pulls out the F statistic
dll.anova.fstat   
dll.anova.fstat1 = dll.anova.fstat[1]

anova(dll.anova)

#### Chunk 4:  A basic permutation test
# and here is the actual code for the permutation test

# while I am doing this with an F statistic below (to compare against theoretical quantities), there is no need to do this. I could just have extracted the coefficient I cared about directly.

N = 999 # # of resampling iterations


PermuteFunction <- function(y = dll.data$SCT, x = dll.data$genotype){
	model.resample = lm(sample(y, replace = F) ~ x) # the important bit! shuffle y relative to x.
	#permutes response, then runs model
	
   fstats = summary(model.resample)$fstat[1] # place all of the fstats in this vector 
   
   return(fstats)}


# Try calling it a few times.
PermuteFunction()

# To perform the iterations of the permutations we can use either replicate()
permute.N <- replicate(N, PermuteFunction())
hist(permute.N)

# or a for loop (useful for big data or lots of iterations)
fstats = numeric(N) # initializes a numeric vector to input the resampled F stats

for (i in 1:N) {
	fstats[i] <- PermuteFunction()}

hist(fstats) 

#let's remind ourselves again about what F really is.

# look at the distribution of the F statistics from permutation.

par(mfrow=c(2,1))
hist(fstats, 
     xlim = c(0, 65), ylim = c(0,1), 
     xlab = "F Statistic", 
     main = "Empirically derived via randomization",
     freq = F)

arrows(dll.anova.fstat1, 0.6, dll.anova.fstat1, 0, lwd = 3)


# plot of F values based on F distribution

f.theor <- rf(1000, 1, 1946) # f distribution
hist(f.theor, 
     xlim = c(0, 65), 
     ylim = c(0, 1),
     xlab = "F Statistic", freq = F,
     main = " F values from F-distribution (1, 1946)")
curve(df(x, 1,1946), 0, 60, add=T)


# Computing the p-value. Note the "1" added both to the numerator (for the observed value, which is equal to itself) and the denominator (as we are accounting for the original estimated value along with our N iterations of shuffling)

(length(fstats[fstats >= dll.anova.fstat1]) + 1)/(N + 1)

 # gives the empirical P value, 
 # in this case none of the permuted values is greater than the observed
 # so the empirical p-value = 0.001 (for 1000 permutations)


#### Chunk 6:  Distribution of p values under permutation

# We can also use this to look at the distribution of p values under this model. Make sure you understand this and why it is uniform.

PermutePvals <- function(){
	dll.resample = lm(sample(SCT, replace = F) ~ genotype, data=dll.data) 
	anova(dll.resample)[1,5]}

permute.N.p <- replicate(N, PermutePvals())

hist(permute.N.p)


# Chunk 7: Using RRPP for more complex  models

## For more complex models the problem is that a simple "shuffling" of the data may break all associations and not just for the particular predictor of interest.

# One alternative is to use a residual based approach for resampling:
# 1. i.e. extract residuals from "null model", 
# 2. shuffle these and add back on to fitted values of null model as our new observations
# 3. Then use these new "observations" and fit under full model
# 4. repeat 2 & 3many times.
# 5. Compare the observed value to those from this permutation as we did above.


mod_GxT <- lm(SCT ~ genotype*temp, data = dll.data)

rrpp_out <- lm.rrpp(mod_GxT, iter = 999, 
                    SS.type = "II",
                    data = dll.data)

# what are the reduced ("null")  model this uses for each term? 
reveal.model.designs(rrpp_out)

# in this case a model of just SCT ~ genotype + temp is used as the "reduced" model to evaluate the interaction term.
# A model of just SCT ~ temp is used as reduced model to evaluate the "genotype" term and vice versa for the "temp" term.

# anova table for permutation
anova(rrpp_out)

# permutation for coefficients (most useful)
coef(rrpp_out, test = TRUE)


#############
# Chunk 8: Some notes on the computation of permuted data sets.


# For univariate tests use sample on the dependent variable only. For models with a single explanatory variable it will not matter, but for models with a number of explanatory variables if you use sample() on the explanatory variables it will break up the intra-observation correlation structure, which is not generally something you want. However you can shuffle the relevant columns of the design matrix.

#  For multivariate tests, you can not use sample() quite as easily. For a simple model where you have multiple dependent variables, and a single explanatory variable use sample() on the explanatory variable. For complex models:

# the basic idea is for matrix "a" with m rows
# a[sample(m,m),]  # allows indexing

# Note that RRPP does handle multivariate models just fine
#####################
