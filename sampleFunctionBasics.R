# BIO708 -Old Screencast is available (from a previous course), resampling in R
# If you want to "watch along"
#  https://www.youtube.com/watch?v=xE3KGVT6VLE

# How the sample() function works.

# In R to both sampling with or without replacement we use the sample() function

# The basic function call is like

#sample(x=your.vector, size=length(your.vector), replace=F)

# So if we have a vector of numbers called y
vals <- 1:10

# This will sample without replacement for y (as we will need for permutation tests)
sample(x = vals, size = length(vals), replace = F)

# Which just returns all of the numbers, but shuffled with respect to order

# By default the size argument will be the length of your input vector so you do not need to specify it if you want if you want to sample the same number of observations as in your original data vector.

sample(x = vals, replace = F) 
sample(x = vals, replace = F) # different ordering each time


# Since this is the same set of observations, with just the order shuffled the mean, sd etc will be the same

mean(vals)
mean(sample(x = vals, replace = F))

sd(vals)
sd(sample(x = vals, replace = F))

# If we want to do sampling with replacement (like with the non parametric bootstrap)
sample(x = vals, replace = T) # Note that sometimes you can draw the same observation multiple times for the original vector

sample(x = vals, replace = T)

# Now the mean and sd may well be different.

mean(vals)
mean(sample(x = vals, replace = T))
sd(vals)
sd(sample(x = vals, replace = T))

# We can (using the "size=" argument) also use this to "roll dice"

y2 <- 1:6 # Six sided die

sample(x = y2, size = 1, replace = T) # single roll of the die
sample(x = y2, size = 10, replace = T)  # ten independent rolls of the 6 sided die


# Sometimes we will want to sample rows (or columns of a matrix), not the individual elements of a vector (or individual elements of the matrix).

# To do this we use an R "trick", namely to use the index 

y3 <- 11:20
Y <- cbind(vals,y3) # make a matrix
Y

# We want to sample rows of this matrix with replacement. 
# Note the index [1,], [2,] on the left.

# We will sample the elements of the index (specifying rows)
# To do this we will specify the number of rows as the first argument in the sample()
# This may seem weird the argument seems like it should be a vector.
# It actually is. But if you specify a single number (like 500) it actually generates a sequence of integers from 1 to 500 (1:500)

# So the sample call will look like this
Y[sample(x=nrow(Y), size = nrow(Y), replace=T), ] # sampling along rows

# Again you can skip the size argument and it will default to the length of x, which in this case in the number of rows of our matrix Y

Y[sample(x=nrow(Y), replace = T), ]
