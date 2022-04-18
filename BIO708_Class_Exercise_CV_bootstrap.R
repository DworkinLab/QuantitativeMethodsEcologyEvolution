
# Change this to effect sizes.

# What can the CV's tell us about the variation between mutants and wild-types?

dll.data = read.csv("https://raw.githubusercontent.com/DworkinLab/DworkinLab.github.io/master/dataSets/Dworkin2005_ED/dll.csv", 
                    header=TRUE, stringsAsFactors = TRUE)

dll.data <- na.omit(dll.data)

CVfunction <- function(x){
	return(sd(x)/mean(x))
}

dll.SCT = dll.data[dll.data$genotype=="Dll", "SCT"]   # notice " == "
wt.SCT = dll.data[dll.data$genotype=="wt", "SCT"]
CV.dll <- CVfunction(dll.SCT)
CV.wt <- CVfunction(wt.SCT)

CV.diff <- CV.dll - CV.wt


#Writing a function to resample the data, calculate CVs, and SE and the difference between CV

bootstrapKitchenSink <- function(x1 = wt.SCT, x2 = dll.SCT) {
	ResampleWildType <- sample(x1, length(x1), replace = T)
	ResampleMutant   <- sample(x2, length(x2), replace = T)
	CV.mut <- CVfunction(ResampleMutant)
	CV.wt  <- CVfunction(ResampleWildType)
	CV.diff.boot <- CV.mut - CV.wt
	return(c(DifferenceCV = CV.diff.boot, CV_mutant = CV.mut, CV_wildtype = CV.wt))
}


R <- 10000


CV.boot.out <- t(replicate(R, bootstrapKitchenSink()))


dim(CV.boot.out)

# Calculate the SE of the difference
apply(X = CV.boot.out, MARGIN = 2, FUN = sd)
hist(CV.boot.out[,1])

# percentile Confidence intervals for the mutant, the wild type and the difference.
apply(X = CV.boot.out, 
      MARGIN = 2, 
      FUN = quantile, probs=c(0.025, 0.975) )

# Bootstrap bias
bootstrap_means <- 
  apply(X = CV.boot.out, MARGIN = 2, FUN = mean)

CV.dll - bootstrap_means["CV_mutant"]

CV.wt - bootstrap_means["CV_wildtype"]
