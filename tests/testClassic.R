
library(rcarbon)
library(ggplot2)   
library(gridExtra)    
source("../jagsCode.R")

getLMParams <- function(model)
{
  l <- list(intercept = as.numeric(coef(model))[1], slope=as.numeric(coef(model))[2], r2 = summary(model)$r.squared)
}

estimateParams<- function(sites)
{
    # direct frequentist model
    mDirect <- lm(arrival~distance, sites)
    mDirectParams <- getLMParams(mDirect)

    # backcalibration
    uncalibrated <- uncalibrate(sites$arrival, calCurves = "intcal13")
    arrivalError = sample(c(10,20,50,100),size=nrow(sites),replace=TRUE)
    # calibration 
    calibratedProbs = calibrate(round(uncalibrated$ccCRA), arrivalError)
    sites$calibrated <- medCal(calibratedProbs)

    # C14 frequentist model
    mWithC14 = lm(calibrated~distance, sites)
    mWithC14Params <- getLMParams(mWithC14)

    # basesian model
    post.samples=jagsSpeed(x=round(uncalibrated$ccCRA),error=arrivalError,distance=sites$distance)
    post.samples=as.data.frame(post.samples[[1]])

    # final dataset
    models <- data.frame(name=c("direct", "withC14", "bayesian"), intercept=c(mDirectParams$intercept, mWithC14Params$intercept, median(post.samples$alpha)), slope=c(mDirectParams$slope, mWithC14Params$slope, median(post.samples$beta)))
    return(models)
}

#### PARAMS

# number of sites
numSites = c(20,40,60,80,100)
iterations = 100

alpha=7000
maxDistance = 2000
# rate of expansion
beta = -0.8
sigma <- c(1,10,25,50,75,100)

allSites <- data.frame(run=integer(), numSites=integer(), sigma=integer(), iteration=integer(), distance=integer(), arrival=integer(), calibrated=integer())  
allEstimates <- data.frame(run=integer(), numSites=integer(), sigma=integer(), iteration=integer(), name=character(), intercept=double(), slope=double())  

run = 0    
totalRuns = iterations*length(numSites)*length(sigma)

for(j in 1:iterations)
{
    for(i in numSites)
    {
        for(z in sigma)
        {
            cat(sprintf("################ starting run %d of %d - iteration %d num sites %d sigma %d\n", run, totalRuns, j, i, z))
            sites <- data.frame(distance=runif(i,min=0,max=maxDistance))
            sites$arrival <- as.integer(rnorm(nrow(sites), mean=alpha+beta*sites$distance, sd=z))
            sites$iteration <- j
            sites$numSites <- i
            sites$sigma <- z
            sites$run <- run

            allSites <- rbind(allSites, sites)

            estimates <- estimateParams(sites)
            estimates$iteration <- j
            estimates$numSites <- i
            estimates$sigma <- z
            estimates$run <- run

            allEstimates <- rbind(allEstimates,estimates)

            run <- run+1
        }

    }
}


# example iteration
ggplot(subset(allSites,iteration==1), aes(x=distance, y=arrival)) + geom_point(alpha=0.3) + facet_grid(numSites~sigma)
# all dataset
ggplot(allSites, aes(x=distance, y=arrival)) + geom_point(alpha=0.01, size=1) + facet_grid(numSites~sigma)
# estimates
ggplot(allEstimates, aes(x=as.factor(numSites), y=slope, fill=name)) + geom_boxplot() + facet_wrap(~sigma, ncol=2)
ggplot(allEstimates, aes(x=as.factor(numSites), y=intercept, fill=name)) + geom_boxplot() + facet_wrap(~sigma, ncol=2)


# I run this 4 times in 4 different consoles to use all CPUs
# these 2 commands should be called on each console (sites1, sites2, sites3 and sites 4).
save(allSites, file="sites1.Rda")
save(allEstimates, file="estimates1.Rda")

# once all 4 runs are finished 
load("sites1.Rda")
allSites1 <- allSites
load("sites2.Rda")
allSites2 <- allSites
load("sites3.Rda")
allSites3 <- allSites
load("sites4.Rda")
allSites4 <- allSites

sites <- rbind(allSites1, allSites2)
sites <- rbind(sites, allSites3)
sites <- rbind(sites, allSites4)


load("estimates1.Rda")
allEstimates1 <- allEstimates
load("estimates2.Rda")
allEstimates2 <- allEstimates
load("estimates3.Rda")
allEstimates3 <- allEstimates
load("estimates4.Rda")
allEstimates4 <- allEstimates

estimates <- rbind(allEstimates1, allEstimates2)
estimates <- rbind(estimates, allEstimates3)
estimates <- rbind(estimates, allEstimates4)

# estimates
ggplot(estimates, aes(x=as.factor(numSites), y=slope, fill=name)) + geom_boxplot() + facet_wrap(~sigma, ncol=2)
ggplot(estimates, aes(x=as.factor(numSites), y=intercept, fill=name)) + geom_boxplot() + facet_wrap(~sigma, ncol=2)


pdf("estimates_slope_base.pdf", width=12, height=18)
ggplot(estimates, aes(x=as.factor(numSites), y=slope, fill=name)) + geom_boxplot() + facet_wrap(~sigma, ncol=2)
dev.off()

pdf("estimates_intercept_base.pdf", width=12, height=18)
ggplot(estimates, aes(x=as.factor(numSites), y=intercept, fill=name)) + geom_boxplot() + facet_wrap(~sigma, ncol=2)
dev.off()

