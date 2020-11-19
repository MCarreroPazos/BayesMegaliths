markTest = function(x,marks,nsim=20)
{
  require(rcarbon)
  x = as.matrix(x)
  x = round(x)
  mcmcnsim = nrow(x)
  nsites = ncol(x)/2
  
  # Define timeRange of analyses as 90% quantile range of start and end dates
  timeRange=quantile(x,c(0.05,0.95))
  # Convert earlier and later dates to timeRange boundaries
  x[which(x<timeRange[1])]=timeRange[1]
  x[which(x>timeRange[2])]=timeRange[2]
  
  # Convert to BP
  timeRange=BCADtoBP(timeRange)
  x = apply(x,2,BCADtoBP)
  
  # Create placeholder array
  timeRangeDuration = abs(diff(timeRange)) + 1
  calBP = timeRange[1]:timeRange[2]
  obsArray = array(NA,dim=c(length(unique(marks)),mcmcnsim,timeRangeDuration))
  simArray = array(NA,dim=c(nsim,length(unique(marks)),timeRangeDuration))
  markNames = unique(marks)
  markIndex = sapply(markNames, function(x,y){which(y==x)},y=marks)
 
  # Main Loop (could be improved in performance but this will do)
  index = seq(from=1,by=2,length.out=nsites)
  print('Compute Observed Time Series...')
  pb <- txtProgressBar(min=0, max=mcmcnsim, style=3)
  
  # Compute observe time-series for all regions
  for (s in 1:mcmcnsim)
  {
    setTxtProgressBar(pb,s)
    tmpMat=matrix(0,nrow=nsites,ncol=timeRangeDuration)
    for (i in 1:nsites)
    {
      st = x[s,index[i]]
      en = x[s,index[i]+1]
      tmpMat[i,which(calBP<=st&calBP>=en)]=1
    }
    obsArray[,s,] = t(sapply(markIndex,function(x,mat){apply(mat[x,],2,sum)},mat=tmpMat))
  }
  # Compute the average time-series across the posterior for all regions
  obsMean = apply(obsArray,c(1,3),mean)
  close(pb)
  
  
  # The bits below is for the permuted set, it would be considerably slow considering that the above takes several minutes and this is basically repeating the procedure nsim times.
  print('Running permutation test...')
  pb <- txtProgressBar(min=0, max=nsim, style=3)
  for (ss in 1:nsim)
  {
    setTxtProgressBar(pb,ss)
    markIndexRandom= sapply(markNames, function(x,y){which(y==x)},y=sample(marks))
    tmpArray = array(NA,dim=c(length(unique(marks)),mcmcnsim,timeRangeDuration))
    for (s in 1:mcmcnsim)
    {
      tmpMat=matrix(0,nrow=nsites,ncol=timeRangeDuration)
      for (i in 1:nsites)
      {
        st = x[s,index[i]]
        en = x[s,index[i]+1]
        tmpMat[i,which(calBP<=st&calBP>=en)]=1
      }
      tmpArray[,s,] = t(sapply(markIndexRandom,function(x,mat){apply(mat[x,],2,sum)},mat=tmpMat))
    }
    simArray[ss,,] = apply(tmpArray,c(1,3),mean)
  }
close(pb)
  
  pValueList = numeric(length=length(markNames))
  resList = vector('list',length=length(markNames))
  for (i in 1:length(markNames))
  {
    obs = obsMean[i,]
    sim = t(simArray[,i,])
    sim.mean = apply(sim,1,mean)
    sim.sd = apply(sim,1,sd)
    sim.z = (sim-sim.mean)/sim.sd
    obs.z = (obs-sim.mean)/sim.sd
  
    hi = apply(sim.z,1,quantile,0.975)
    lo = apply(sim.z,1,quantile,0.025)
    
    obs.stat = sum(lo[which(obs.z<=lo)]-obs.z[which(obs.z<=lo)]) + sum(obs.z[which(obs.z>=hi)]-hi[which(obs.z>=hi)])
    sim.stat=apply(sim.z,2,function(x,lo,hi){return(sum(lo[which(x<=lo)]-x[which(x<=lo)]) + sum(x[which(x>=hi)]-hi[which(x>=hi)]))},lo=lo,hi=hi)
    r=sum(obs.stat<= sim.stat) 
    pValueList[i]=(r+1)/(nsim+1)
    resList[[i]]=data.frame(calBP=calBP,obs=obs,hi=apply(sim,1,quantile,0.975),lo=apply(sim,1,quantile,0.025))
  }

 return(list(pValueList=pValueList,resList=resList))
}



plotMarkTest = function(x,index,...)
{
  x = x$resList[[index]]
  

  
  # Boom and Bust Handling ####
  booms <- which(x$obs>x$hi)
  busts <- which(x$obs<x$lo)
  baseline <- rep(NA,nrow(x))
  colpts = rep('grey',nrow(x))
  colpts[booms] = 'red'
  colpts[busts] = 'blue'
  
  boomPlot <- baseline
  if (length(booms)>0){ boomPlot[booms]=x$obs[booms] }
  bustPlot <- baseline
  if (length(busts)>0){ bustPlot[busts]=x$obs[busts] }           
  
  boomBlocks <- vector("list")
  counter <- 0
  state <- "off"
  for (i in 1:length(boomPlot)){
    if (!is.na(boomPlot[i])&state=="off"){
      counter <- counter+1
      boomBlocks <- c(boomBlocks,vector("list",1))
      boomBlocks[[counter]] <- vector("list",2)
      boomBlocks[[counter]][[1]] <- boomPlot[i]
      boomBlocks[[counter]][[2]] <- x[i,"calBP"]
      state <- "on"
    }
    if (state=="on"){
      if (!is.na(boomPlot[i])){
        boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[i])
        boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],x[i,"calBP"])
      }
      if (is.na(boomPlot[i])){
        state <- "off"
      }
    }   
  }
  
  bustBlocks <- vector("list")
  counter <- 0
  state <- "off"
  for (i in 1:length(bustPlot)){
    if (!is.na(bustPlot[i])&state=="off"){
      counter <- counter+1
      bustBlocks <- c(bustBlocks,vector("list",1))
      bustBlocks[[counter]] <- vector("list",2)
      bustBlocks[[counter]][[1]] <- bustPlot[i]
      bustBlocks[[counter]][[2]] <- x[i,"calBP"]
      state <- "on"
    }
    if (state=="on"){
      if (!is.na(bustPlot[i])){
        bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[i])
        bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],x[i,"calBP"])
      }
      if (is.na(bustPlot[i])){
        state <- "off"
      }
    }   
  }
  
  
  
  
  
  # Actual Plot ####
  plot(0,0,type='n',xlim=rev(range(x$calBP)),ylim=range(c(x$obs,x$lo,x$hi),na.rm=T),xlab='cal BP',ylab='Frequency',...)
  polygon(c(x$calBP,rev(x$calBP)),c(x$lo,rev(x$hi)),col='lightgrey',border=NA)
  
  if (length(booms)>0){
    for (i in 1:length(boomBlocks)){
      bbb = unique(boomBlocks[[i]][[2]])
      index = which(x$calBP%in%bbb)
      polygon(c(bbb,rev(bbb)),c(x$obs[index],rev(x$hi[index])),border=NA,col=rgb(0.8,0.36,0.36,0.5))
    }  
  }
  
  if (length(busts)>0){
    for (i in 1:length(bustBlocks)){
      bbb = unique(bustBlocks[[i]][[2]])
      index = which(x$calBP%in%bbb)
      polygon(c(bbb,rev(bbb)),c(x$obs[index],rev(x$lo[index])),border=NA,col=rgb(0.25,0.41,0.88,0.5))
    }  
  }
  lines(x$calBP,x$obs,lwd=1.5)
}



