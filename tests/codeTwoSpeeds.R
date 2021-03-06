
source("../jagsCode.R")

# 1- linear model in calendar time

## params

#number of samples
N=500

# percentages for speed1 and speed 2
percSp1 = 0.6
percSp2 = 1-percSp1

threshold <- as.integer(N*percSp1)

#distance from origin
d=runif(N,min=0,max=5000)

# divide dataset before and after stopping point
dSp1 <- d[1:threshold]
dSp2 <- d[(threshold+1):length(d)]

#intercept (date of earliest farming at the origin)  
alpha= 8000
# speeds
betaSp1 = -0.8
betaSp2 = -0.6

sigmaSp1=200
sigmaSp2=300

# generation of possible dates of arrival
arrivalSp1 <- rnorm(length(dSp1),mean=alpha+betaSp1*dSp1,sd=sigmaSp1)
arrivalSp2 <- rnorm(length(dSp2),mean=alpha+betaSp2*dSp2,sd=sigmaSp2)

# get the sorted list of arrivals and distances    
arrival <- c(arrivalSp1, arrivalSp2)
plot(d,arrival,xlab="distance from origin",ylab="calendar Year",pch=20)

# theorethical model
abline(a=alpha,b=betaSp1,lty=2)
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

hist(post.samples$beta,border=NA,col="grey",main="Posterior of Beta",xlab="beta")
abline(v=betaBeforeStop,lty=2,col="indianred",lwd=2)
abline(v=betaAfterStop,lty=2,col="indianred",lwd=2)
abline(v=fitted.beta,lty=3,col="orange",lwd=2)
abline(v=mMedians.beta,lty=3,col="blue",lwd=2)

hist(sqrt(post.samples$gamma^-1),border=NA,col="grey",main="Posterior of Sigma",xlab="sigma")
abline(v=sigma,lty=2,col="indianred")


