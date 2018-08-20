
source("../jagsCode.R")

library(rethinking)
library(rcarbon)

# 0 example test of dating a C14 sample
# par(mfrow=c(2,1))
# plot(calibrate(4500,30)) #using rcarbon
# samples=jagsCalibrate(4500,0.0001,raw=TRUE,plot=FALSE)
# hist(as.numeric(samples[[1]]),xlim=rev(range(samples[[1]])),border=NA,col="grey",breaks=200,main="")



# 1- linear model in calendar time

## params

#number of samples
N=100

#distance from origin
d=runif(N,min=0,max=5000)

#intercept (date of earliest farming at the origin)  
alpha=8000
#beta coefficient (i.e. rate of expansion)

betaM = -0.8
betaSD = 0.05
beta=rnorm(N, betaM, betaSD)
#standard deviation (uncertainty in the arrival date at given distance)
sigma=200

# generation of possible dates of arrival
#arrival <- rnorm(N,mean=alpha+beta*d,sd=sigma)
arrival <- alpha+beta*d

plot(d,arrival,xlab="distance from origin",ylab="calendar Year",pch=20)

# theorethical model
abline(a=alpha,b=beta,lty=2)
logModel <- lm(arrival~d)
# regression based on data
abline(logModel,lty=3)

fitted.alpha <- coefficients(logModel)[1]
fitted.beta <- coefficients(logModel)[2]

# 2- calibrate calendar to C14 dates
uncalibrated <- uncalibrate(arrival, calCurves = "intcal13")

arrivalC14 = round(uncalibrated$ccCRA)
arrivalError = sample(c(10,20,50,80,100),size=N,replace=TRUE)

# 3 - standard approach based on median
arrivalCalibrated = calibrate(arrivalC14, arrivalError)
arrivalMedian = medCal(arrivalCalibrated)

mMedians = lm(arrivalMedian~d)
summary(mMedians)

mMedians.alpha <- coefficients(mMedians)[1]
mMedians.beta <- coefficients(mMedians)[2]

fitted.alpha
fitted.beta

mMedians.alpha
mMedians.beta 

# 4 - speed estimation based on the model
post.samples=jagsSpeed(x=arrivalC14,error=arrivalError,distance=d)
post.samples=as.data.frame(post.samples[[1]])

# 5 - plot results for bayesian stuff

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


# 6 - plot comparison
par(mfrow=c(1,3))
hist(post.samples$alpha,border=NA,col="grey",main="Posterior of Alpha",xlab="alpha")
abline(v=alpha,lty=1,col="indianred",lwd=2)
abline(v=fitted.alpha,lty=1,col="orange",lwd=2)
abline(v=mMedians.alpha,lty=1,col="blue",lwd=2)
abline(v=median(post.samples$alpha),lty=1,col="palegreen4",lwd=2)
legend("topleft",legend=c("true","fitted", "classic", "bayesian"),lty=1,col=c("indianred","orange", "blue", "palegreen4"))

hist(post.samples$beta,border=NA,col="grey",main="Posterior of Beta",xlab="beta")
abline(v=beta,lty=1,col="indianred",lwd=2)
abline(v=fitted.beta,lty=1,col="orange",lwd=2)
abline(v=mMedians.beta,lty=1,col="blue",lwd=2)
abline(v=median(post.samples$beta),lty=1,col="palegreen4",lwd=2)

hist(sqrt(post.samples$gamma^-1),border=NA,col="grey",main="Posterior of Sigma",xlab="sigma")
abline(v=sigma,lty=2,col="indianred")
abline(v=median(sqrt(post.samples$gamma^-1)),lty=1,col="palegreen4")

# show median values for bayesian posterior
mMedians.alpha
mMedians.beta 


median(post.samples$alpha)
median(post.samples$beta)
median(sqrt(post.samples$gamma^-1))


