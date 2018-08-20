
library(rcarbon)
library(ggplot2)   


getLMParams <- function(model)
{
  l <- list(intercept = as.numeric(coef(model))[1], slope=as.numeric(coef(model))[2], r2 = summary(model)$r.squared)
}

estimateParams<- function(sites, excavatedPerc)
{
    # generate a sample of excavated sites
    numExcavatedSites = as.integer(nrow(sites)*excavatedPerc)
    excavatedRows <- sample(nrow(sites), numExcavatedSites)
    excavatedSites <- sites[excavatedRows,]

    # direct frequentist model
    mDirect <- lm(arrival~distance, excavatedSites)
    mDirectParams <- getLMParams(mDirect)

    # backcalibration
    uncalibrated <- uncalibrate(excavatedSites$arrival, calCurves = "intcal13")
    arrivalC14 = round(uncalibrated$ccCRA)
    arrivalError = sample(c(10,20,50,100),size=numExcavatedSites,replace=TRUE)
    # calibration 
    arrivalCalibrated = calibrate(arrivalC14, arrivalError)

    # C14 frequentist model
    mWithC14 = lm(arrivalC14~distance, excavatedSites)
    mWithC14Params <- getLMParams(mWithC14)

    # basesian model
    post.samples=jagsSpeed(x=arrivalC14,error=arrivalError,distance=excavatedSites$distance)
    post.samples=as.data.frame(post.samples[[1]])

    # final dataset
    models <- data.frame(name=c("direct", "withC14", "bayesian"), intercept=c(mDirectParams$intercept, mWithC14Params$intercept, median(post.samples$alpha)), slope=c(mDirectParams$slope, mWithC14Params$slope, median(post.samples$beta)), percentage=excavatedPerc)
    return(models)
}

#### PARAMS

# number of sites
N=10000
# range of excavated fractions 
#excavatedPerc = c(0.01,0.02,0.03,0.04,0.05,0.1,0.25,0.5)
excavatedPerc = c(0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
iterations = 10    

alpha=8000
# rate of expansion
betaM = -0.8
betaSD = 0.1
beta=rnorm(N, betaM, betaSD)


#distance from origin
sites <- data.frame(distance=runif(N,min=0,max=5000))
# date of earliest farming at the origin
# generation of possible dates of arrival
sites$arrival=alpha+beta*sites$distance

estimates <- data.frame(name=character(), intercept=double(), slope=double, percentage=double)  
for(j in 1:iterations)
{
    cat(sprintf("#################### starting iteration %d of %d\n", j, iterations))
    for(i in excavatedPerc)
    {
        cat(sprintf("\t################ excavated fraction %f\n", i))
        estimates <- rbind(estimates, estimateParams(sites, i))
    }
}



ggplot(estimates, aes(x=percentage, y=slope, col=name)) + geom_point() + geom_hline(yintercept=betaM)

# ggplot(estimates, aes(x=percentage, y=slope, col=name)) + geom_jitter(width=0.005, height=0) + geom_hline(yintercept=betaM) + xlim(0,0.1)
#ggplot(sites, aes(x=distance, y=arrival)) + geom_point(data=sites, alpha=0.1) + geom_point(data=excavatedSites, size=1) +geom_abline(data=models, aes(slope=slope, intercept=intercept, col=name), size=1, linetype="dashed")

