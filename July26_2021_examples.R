
# Simple regression

size <- c(3, 3.1, 2, 1.5, 2.3, 3.5, 3.6, 2.1)
protein <- c(1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 1.85)

plot(y = size, x = protein)

mod1 <- lm(size ~ protein)
abline(mod1)
summary(mod1) # estimated intercept, and slope

model.matrix(mod1) # design matrix for the simple regression

model.matrix(~protein) # the design matrix ONLY depends on the predictors.

## second example with two continuous predictors.

# now we have carbohydrates as well.

carbs <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 0.85)

mod2 <- lm( size ~ 1 + protein + carbs)

model.matrix(mod2)


# third example with a single discrete predictor (to demonstrate indicator/dummy variable coding)

# our predictor
water <- c("lake", "lake", "lake", "pond", "pond", "pond")

# our response variables (we model them seperately)

Nitrogen <- c(2,2.5, 2.1, 3.1, 3.1, 2.8 )

Phos <- c(2,2.5, 2.1, 1.6, 1.1, 1.8 )

mod_water <- lm(Nitrogen ~ 1 + water) # remember the "1" is there by default, so even Nitrogen ~ water gives same model
summary(mod_water) 

model.matrix(mod_water) # treatment contrast coding.

#means of each group of water
tapply(Nitrogen, water, mean)


## same idea with a second response, but same predictor (design matrix does not change)
mod_water_P <- lm(Phos ~ 1 + water)
summary(mod_water_P)
tapply(Phos, water, mean)
model.matrix(mod_water_P)

# how about cell means model.
# just need to suppress the global intercept.
# note, this is not the real way of doing cell means, but it is fine for now.

mod_water_cell <- lm(Nitrogen ~ 0 + water)
summary(mod_water_cell)
model.matrix(mod_water_cell)
