jagsCalibrate<-function(x,error,calCurves='intcal13',init=NA,raw=FALSE,plot=TRUE)
{

	require(rcarbon)
	calCurveFile <- paste(system.file("extdata", package="rcarbon"), "/", calCurves,".14c", sep="")
	options(warn=-1)
	cctmp <- readLines(calCurveFile, encoding="UTF-8")
	cctmp <- cctmp[!grepl("[#]",cctmp)]
	cctmp <- as.matrix(read.csv(textConnection(cctmp), header=FALSE, stringsAsFactors=FALSE))[,1:3]
	options(warn=0)
	colnames(cctmp) <- c("CALBP","C14BP","Error")


	require(rjags) 
	calBP=cctmp[,1]
	C14BP=cctmp[,2]
	C14err=cctmp[,3]
	dataList=list(nDate=length(x), X=x, sigma=error, calBP=rev(calBP), C14BP=rev(C14BP),C14err=rev(C14err))



	if (is.na(init))
	{
		for (i in 1:length(x))
		{init[i]=calBP[which(abs(C14BP-x[i])==min(abs(C14BP-x[i])))[1]]}
	}
	##Specify JAGS Model
	modelString="
	model{
		for (i in 1:nDate) {
			theta[i] ~ dunif(0,49980)
			mu[i] <- interp.lin(theta[i], calBP[], C14BP[])
			sigmaCurve[i] <- interp.lin(theta[i], calBP[], C14err[])
			tau[i] <- 1/(pow(sigma[i],2)+pow(sigmaCurve[i],2))
			X[i] ~ dnorm(mu[i],tau[i])
			twenty.year[i] <- 20*round(theta[i]/20)
			ten.year[i] <- 10*round(theta[i]/10)
			one.year[i] <- round(theta[i])
		}

	}"

	initsList=list(theta=init)
	jagsModel=jags.model(file=textConnection(modelString),data=dataList,inits=initsList,n.chains=1)
	update(jagsModel,n.iter=3000)
	codaSamples=coda.samples(jagsModel,variable.names=c("one.year"),n.iter=50000)
	if (length(x)>1){plot=FALSE}
	if (plot==TRUE)
		{plot(density(codaSamples[[1]],bw=1),xlim=c(max(codaSamples[[1]]),min(codaSamples[[1]])),xlab="cal BP",ylab="density",main="")}
	if(raw==TRUE){return(codaSamples)}
}


jagsSpeed<-function(x,error,distance,calCurves='intcal13',init=NA,raw=FALSE,plot=TRUE)
{

	require(rcarbon)
	calCurveFile <- paste(system.file("extdata", package="rcarbon"), "/", calCurves,".14c", sep="")
	options(warn=-1)
	cctmp <- readLines(calCurveFile, encoding="UTF-8")
	cctmp <- cctmp[!grepl("[#]",cctmp)]
	cctmp <- as.matrix(read.csv(textConnection(cctmp), header=FALSE, stringsAsFactors=FALSE))[,1:3]
	options(warn=0)
	colnames(cctmp) <- c("CALBP","C14BP","Error")


	require(rjags) 
	calBP=cctmp[,1]
	C14BP=cctmp[,2]
	C14err=cctmp[,3]
	dataList=list(nDate=length(x), X=x, d=distance, sigma=error, calBP=rev(calBP), C14BP=rev(C14BP),C14err=rev(C14err))



	if (is.na(init))
	{
		for (i in 1:length(x))
		{init[i]=calBP[which(abs(C14BP-x[i])==min(abs(C14BP-x[i])))[1]]}
	}
	##Specify JAGS Model
	modelString="
	model{
		for (i in 1:nDate) {
			theta[i] ~ dnorm(mean[i],gamma)	
		        mean[i] <- alpha + beta * d[i]			
			mu[i] <- interp.lin(theta[i], calBP[], C14BP[])
			sigmaCurve[i] <- interp.lin(theta[i], calBP[], C14err[])
			tau[i] <- 1/(pow(sigma[i],2)+pow(sigmaCurve[i],2))
			X[i] ~ dnorm(mu[i],tau[i])
			twenty.year[i] <- 20*round(theta[i]/20)
			ten.year[i] <- 10*round(theta[i]/10)
			one.year[i] <- round(theta[i])
		}
                alpha ~ dunif(0,50000)
		beta ~ dnorm(0,10)
		gamma <- pow(pregamma, -2)
                pregamma ~ dunif(0,500)

	}"

	initsList=list(theta=init)
	jagsModel=jags.model(file=textConnection(modelString),data=dataList,inits=initsList,n.chains=1)
	update(jagsModel,n.iter=3000)
	codaSamples=coda.samples(jagsModel,variable.names=c("one.year","beta","alpha","gamma"),n.iter=50000)
        return(codaSamples)
}

