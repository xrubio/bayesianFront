
library(rcarbon)
library(ggplot2)   
library(gridExtra)    

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
    arrivalError = sample(c(10,20,50,100),size=numExcavatedSites,replace=TRUE)
    # calibration 
    calibratedProbs = calibrate(round(uncalibrated$ccCRA), arrivalError)
    excavatedSites$calibrated <- medCal(calibratedProbs)

    # C14 frequentist model
    mWithC14 = lm(calibrated~distance, excavatedSites)
    mWithC14Params <- getLMParams(mWithC14)

    # basesian model
    post.samples=jagsSpeed(x=round(uncalibrated$ccCRA),error=arrivalError,distance=excavatedSites$distance)
    post.samples=as.data.frame(post.samples[[1]])

    # final dataset
    models <- data.frame(name=c("direct", "withC14", "bayesian"), intercept=c(mDirectParams$intercept, mWithC14Params$intercept, median(post.samples$alpha)), slope=c(mDirectParams$slope, mWithC14Params$slope, median(post.samples$beta)), percentage=excavatedPerc)
    return(models)
}

#### PARAMS

# number of sites
N=1000
# range of excavated fractions 
#excavatedPerc = c(0.01,0.02,0.03,0.04,0.05,0.1)
excavatedPerc = c(0.01,0.03,0.05,0.1)
iterations = 100

alpha=3000
# rate of expansion
betaM = -0.8
#betaSD = 0.05
#beta=rnorm(N, betaM, betaSD)
sigma <- 10

#distance from origin
sites <- data.frame(distance=runif(N,min=0,max=2000))
# date of earliest farming at the origin
# generation of possible dates of arrival
sites$arrival <- as.integer(rnorm(nrow(sites), mean=alpha+betaM*sites$distance, sd=sigma))

ggplot(sites, aes(x=distance, y=arrival)) + geom_point()


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




ggplot(estimates, aes(x=as.factor(percentage), y=slope, fill=name)) + geom_boxplot() + geom_hline(yintercept=betaM)

ggplot(estimates, aes(x=as.factor(percentage), y=intercept, fill=name)) + geom_boxplot() + geom_hline(yintercept=alpha)

pdf("sigma1b.pdf", width=10, height=25)
g1 <- ggplot(sites, aes(x=distance, y=arrival)) + geom_point(alpha=0.25)
g2 <- ggplot(estimates, aes(x=as.factor(percentage), y=slope, fill=name)) + geom_boxplot() + geom_hline(yintercept=betaM)
g3 <- ggplot(estimates, aes(x=as.factor(percentage), y=intercept, fill=name)) + geom_boxplot() + geom_hline(yintercept=alpha)
grid.arrange(g1,g2,g3)
dev.off()


