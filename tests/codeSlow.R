
source("../jagsCode.R")

# 1- linear model in calendar time

## params

#number of samples
N=500
stopPoint = 2000

#distance from origin
d=runif(N,min=0,max=5000)

# divide dataset before and after stopping point
dBeforeStop <- d[d<stopPoint]
dAfterStop <- d[d>=stopPoint]

## diffusion before stopPoint
#intercept (date of earliest farming at the origin)  
alphaBeforeStop = 8000
#beta coefficient (i.e. rate of expansion)
betaBeforeStop=-0.8
#standard deviation (uncertainty in the arrival date at given distance)

## diffusion after stopPoint
alphaAfterStop = alphaBeforeStop+betaBeforeStop*stopPoint
betaAfterStop = -0.4

sigma=200

# generation of possible dates of arrival
arrivalBeforeStop <- rnorm(length(dBeforeStop),mean=alphaBeforeStop+betaBeforeStop*dBeforeStop,sd=sigma)
arrivalAfterStop <- rnorm(length(dAfterStop),mean=alphaAfterStop+betaAfterStop*(dAfterStop-stopPoint),sd=sigma)

# get the sorted list of arrivals and distances    
arrival <- c(arrivalBeforeStop, arrivalAfterStop)
d <- c(dBeforeStop, dAfterStop)
plot(d,arrival,xlab="distance from origin",ylab="calendar Year",pch=20)

# theorethical model
abline(a=alpha,b=beta,lty=2)
logModel <- lm(arrival~d)
# regression based on data
abline(logModel,lty=3)

fitted.alpha <- coefficients(logModel)[1]
fitted.beta <- coefficients(logModel)[2]


# 2- calibrate calendar to C14 dates
library(rcarbon)
uncalibrated <- uncalibrate(arrival, calCurves = "intcal13")

arrivalC14 = round(uncalibrated$ccCRA)
arrivalError = sample(c(20,50,80,100),size=N,replace=TRUE)

# 3 - speed estimation based on the model
post.samples=jagsSpeed(x=arrivalC14,error=arrivalError,distance=d)
post.samples=as.data.frame(post.samples[[1]])

# 4 - plot results for bayesian stuff

dates=apply(post.samples,2,quantile,prob=c(0.025,0.975))[,-c(1:3)]
middates=apply(post.samples,2,median)[-c(1:3)]
plot(d,middates,ylim=range(dates),pch=20,xlab="Distance from Origin",ylab="Calibrated Dates")
for (i in 1:N)
{
	lines(rep(d[i],2),dates[,i])
}
dd=seq(0,10000,200)
obsMean=matrix(NA,50000,length(dd))

for (x in 1:length(dd))
{
    obsMean[,x]=post.samples$alpha+post.samples$beta*dd[x]
}
mu=apply(obsMean,2,mean)
ci=apply(obsMean,2,PI)
lines(dd,mu)
shade(ci,dd)


# 4 - 5tandard approach based on meadian
arrivalCalibrated = calibrate(arrivalC14, arrivalError)
arrivalMedian = medCal(arrivalCalibrated)

mMedians = lm(arrivalMedian~d)
summary(mMedians)

mMedians.alpha <- coefficients(mMedians)[1]
mMedians.beta <- coefficients(mMedians)[2]


# 6 - plot comparison
par(mfrow=c(1,3))
hist(post.samples$alpha,border=NA,col="grey",main="Posterior of Alpha")
abline(v=alphaBeforeStop,lty=2,col="indianred",lwd=2)
abline(v=fitted.alpha,lty=3,col="orange",lwd=2)
abline(v=mMedians.alpha,lty=3,col="blue",lwd=2)
legend("topleft",legend=c("true","fitted", "classic"),lty=c(2,3),col=c("indianred","orange", "blue"))

hist(post.samples$beta,border=NA,col="grey",main="Posterior of Beta",xlab="beta",xlim=c(-0.3, -0.9))
abline(v=betaBeforeStop,lty=2,col="indianred",lwd=2)
abline(v=betaAfterStop,lty=2,col="indianred",lwd=2)
abline(v=fitted.beta,lty=3,col="orange",lwd=2)
abline(v=mMedians.beta,lty=3,col="blue",lwd=2)

hist(sqrt(post.samples$gamma^-1),border=NA,col="grey",main="Posterior of Sigma",xlab="sigma")
abline(v=sigma,lty=2,col="indianred")


