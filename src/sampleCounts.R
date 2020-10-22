sampleCounts = function(x)
{
	require(rcarbon)
	x = as.matrix(x)
	x = round(x)
	nsim = nrow(x)
	nsites = ncol(x)/2
	
	# Define timeRange of analyses as 90% quantile range of start and end dates
	timeRange=quantile(x,c(0.05,0.95))
	# Convert earlier and later dates to timeRange boundaries
	x[which(x<timeRange[1])]=timeRange[1]
	x[which(x>timeRange[2])]=timeRange[2]

	# Convert to BP
	timeRange=BCADtoBP(timeRange)
	x = apply(x,2,BCADtoBP)

	# Create placeholder matrix
	timeRangeDuration = abs(diff(timeRange)) + 1
	cntMat = matrix(NA,nrow=nsim,ncol=timeRangeDuration)

	# Main Loop (could be improved in performance but this will do)
	index = seq(from=1,by=2,length.out=nsites)
	calBP = timeRange[1]:timeRange[2]
	colnames(cntMat)=calBP
	pb <- txtProgressBar(min=0, max=nsim, style=3)
	for (s in 1:nsim)
	{
		setTxtProgressBar(pb,s)
		tmpMat=matrix(0,nrow=nsites,ncol=timeRangeDuration)
		for (i in 1:nsites)
		{
			st = x[s,index[i]]
			en = x[s,index[i]+1]
			tmpMat[i,which(calBP<=st&calBP>=en)]=1
		}
		cntMat[s,]=apply(tmpMat,2,sum)
	}
	close(pb)
	return(cntMat)
}


plotCntMat = function(x)
{
	calBP = as.numeric(colnames(x))
	avgCnts = apply(x,2,mean)
	lo = apply(x,2,quantile,0.025)
	hi = apply(x,2,quantile,0.975)
	plot(0,0,xlim=rev(range(calBP)),ylim=c(0,max(hi)),xlab='Cal BP',ylab='Counts')
	polygon(c(calBP,rev(calBP)),c(lo,rev(hi)),col='lightblue',border=NA)
	lines(calBP,avgCnts,col='darkblue')
}







