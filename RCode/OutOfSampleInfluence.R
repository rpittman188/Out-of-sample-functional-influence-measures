#Out of sample influence Code

library(fdaconcur) #needed to write in the datasets
library(fda)
library(pracma)
library(lubridate)
library(MASS) #for mvrnorm
library(corpcor)


#Section 5 Code
fulloutput<-PredictFRegressNormTestUpdatedForRUpdate(FinalXtL1star, FinalYtL1star, Oct15CongHt$CongHt, nBasis = 11,Plot = FALSE)
par(mfrow=c(1,1))

basisList<-11 #going to keep it consistent
EachLeftOut<-Leave1Outconcurrent(PredictorMat=FinalXtL1star, ResponseMat=FinalYtL1star, PredictorVec=Oct15CongHt$CongHt, BasisSelection = "L2", basislist = basisList, plot = FALSE)
par(mfrow=c(1,1))
for(i in 1:10){
  plot(EachLeftOut$`Predicted Responses Without i`[,i], lwd=2,ylab="River Height (ft)", col = "red",type ="l", ylim = c(0,20),main = "Reconstructed 2015 Cedar Creek Stage" )
  lines(fulloutput$PredictedResponse,lwd=2)
  legend( x="topright",
          legend=c("With all 10 Events", paste("With Event", i, "left out")), col = c(1,2), lty = c(1,1),lwd=c(2,2))
}

#Figure 9:
plot(EachLeftOut$`Predicted Responses Without i`[,10], lwd=2, lty = 2,ylab="River Height (ft)", col = "red",type ="l", ylim = c(0,20),main = "Reconstructed 2015 Cedar Creek Stage" )
lines(fulloutput$PredictedResponse,lwd=2)
legend( x="topright",
        legend=c("With all 10 Observations", paste("With Feb. 20 Observation left out")), col = c(1,2), lty = c(1,2),lwd=c(2,2))


#now want L2 distance between fulloutput and 1 lefto out
ri<-c()

squaredifference<-matrix(NA,nrow = nrow(EachLeftOut$`Predicted Responses Without i`),ncol = ncol(EachLeftOut$`Predicted Responses Without i`))
for(i in 1:10){
  squaredifference[,i]<-(fulloutput$PredictedResponse - EachLeftOut$`Predicted Responses Without i`[,i])^2
  ri[i]<-sqrt(trapz(x=1:length(squaredifference[,i]), y = squaredifference[,i])) #need
}

deltaPfunFour<-function(xData,yData,predvec,BasisNum=11, lambda=10^(-1))
{
  FullResults<-PredictFRegressNormTestUpdatedForRUpdate(PredictorMat =xData,ResponseMat=yData,predVec=predvec,nBasis=BasisNum, lambda=10^(-1))$PredictedResponse
  i_leftout<-matrix(NA,nrow = nrow(xData),ncol = ncol(xData))
  for(i in 1:ncol(i_leftout)){
    i_leftout[,i]<-PredictFRegressNormTestUpdatedForRUpdate(PredictorMat =xData[,-i],ResponseMat=yData[,-i],predVec=predvec,nBasis=BasisNum, lambda=10^(-1))$PredictedResponse #calculate the fitted vallues for predictor vec with each left out
  }
  ri<-c()
  squaredifference<-matrix(NA,nrow = nrow(i_leftout),ncol = ncol(i_leftout))
  for(i in 1:ncol(i_leftout)){
    squaredifference[,i]<-(FullResults - i_leftout[,i])^2
    ri[i]<-sqrt(trapz(x=1:length(squaredifference[,i]), y = squaredifference[,i]))
  }
  out<-ri
  return(out)
}

##Code needed to Table 2
alpha = 0
theta<-as.vector((1/ri)^alpha/sum((1/ri)^alpha))
N=ncol(FinalXtL1star)
nboot = 10
boot_deltaP<-c()
for(i in 1:nboot){
  index_samp<-sample(1:N,N,replace=TRUE,prob = theta) #resamples the needed data with the given probability
  temp_deltaP<-deltaPfunFour(xData = FinalXtL1star[,index_samp],yData = FinalYtL1star[,index_samp],BasisNum = 11, predvec = Oct15CongHt$CongHt)
  boot_deltaP<-c(boot_deltaP,temp_deltaP)
}
out0<-c(quantile(boot_deltaP,0.9),quantile(boot_deltaP,0.95),quantile(boot_deltaP,0.99))


#Now doing phi study
absdifference<-matrix(NA,nrow = nrow(EachLeftOut$`Predicted Responses Without i`),ncol = ncol(EachLeftOut$`Predicted Responses Without i`))
for(i in 1:10){
  absdifference[,i]<-abs(fulloutput$PredictedResponse - EachLeftOut$`Predicted Responses Without i`[,i])
}

#Figure 10
plot(absdifference[,1],type="l",lwd = 3, ylab = "Absolute Difference (ft)",main = "Absolute Differences", ylim = c(0,1.5))
for(i in 2:10){
  lines(absdifference[,i],col=i,lwd=3)
}

#Figure 4.11. include "plot" and "lines" below
par(mfrow=c(1,1))
phi=seq(.01,1,by=.01)
#plot(x=phi, y = quantile(absdifference[,3],phi),lwd=2,type="l",xlab="Percentile",ylab="Absolute Difference (ft)",main = "Absolute Difference Percentiles")
hold<-c()
for(i in 1:10) {
  #  lines(phi,quantile(absdifference[,i],phi),lwd=2,col=i)
  hold<-c(hold,trapz(x=1:length(phi),y=quantile(absdifference[,i],phi)))
}

#Table 3 results
area_abs_diff=hold

##Lets put the bootstrap phi right here
##Bootstrapping River stage AIP
##No perturbations needed. Four means it uses Fourier basis functions
area_abs_diffFour<-function(xData,yData,predvec,BasisNum=11, lambda=10^(-1))
{
  FullResults<-PredictFRegressNormTestUpdatedForRUpdate(PredictorMat =xData,ResponseMat=yData,predVec=predvec,nBasis=BasisNum, lambda=10^(-1))$PredictedResponse
  i_leftout<-matrix(NA,nrow = nrow(xData),ncol = ncol(xData))
  for(i in 1:ncol(i_leftout)){
    i_leftout[,i]<-PredictFRegressNormTestUpdatedForRUpdate(PredictorMat =xData[,-i],ResponseMat=yData[,-i],predVec=predvec,nBasis=BasisNum, lambda=10^(-1))$PredictedResponse #calculate the fitted vallues for predictor vec with each left out
  }
  area_measure<-c()
  absdifference<-matrix(NA,nrow = nrow(i_leftout),ncol = ncol(i_leftout))
  for(i in 1:ncol(i_leftout)){
    absdifference[,i]<-abs(FullResults - i_leftout[,i])
  }
  phi=seq(.01,1,by=.01)
  for(i in 1:10) {
    area_measure<-c(area_measure,trapz(x=1:length(phi),y=quantile(absdifference[,i],phi)))
  }
  out<-area_measure
  return(out)
}

test=area_abs_diffFour(FinalXtL1star,FinalYtL1star,predvec = Oct15CongHt$CongHt)

#adds perturbations if wanted to the response curves in each sample
OUonFloodEM<-function(sigma=1,drift=1,times,endtime,Data){
  dt<-endtime/times; #want this to be around .5 for most situations to keep it stable enough
  objects<-ncol(Data) #number of flood events I have
  normalvec<-matrix(rnorm(objects*times,sd=1),nrow=objects,ncol=times); #this is randn from matlab code
  ouvalue<-matrix(0:0,nrow=objects,ncol=times) #null matrix of 0s
  
  # Defining the cluster sizes and the true clustering structure: No clusters here so just 1
  timevec<-seq(0,endtime,length=times)
  
  # loop to create approximation of Ornstein-Uhlenbeck process	# Using Euler-Maruyama approximation
  ouvalue[1:objects,1]=normalvec[,1] #So we don't start at 0.
  for (timeint in 2:times)
  {
    ouvalue[1:objects,timeint]<- ouvalue[1:objects,timeint-1]-drift*ouvalue[1:objects,timeint-1]*dt + sigma*normalvec[1:objects,timeint]*sqrt(dt) #normalvec[1:clussizes[1],timeint] n random values from normal(0, sigma)
  }
  
  samesizeOU<-matrix(0:0,nrow=nrow(Data),ncol=objects) #makes sure it is the same size as the vector it is being applied to
  for(i in 1:objects){
    samesizeOU[,i]<-approx(x=1:ncol(ouvalue),y = ouvalue[i,], n=nrow(Data))$y
  }
  # Defining and naming the observed curves
  obscurves<-Data+samesizeOU;
  out<-list(1:nrow(Data),obscurves,samesizeOU)
  names(out)<-c("timevec","obscurves","samesizeOU")
  return(out)
}

mar03dphi<-trapz(x=1:length(phi),y=quantile(absdifference[,3],phi))
feb20dphi<-trapz(x=1:length(phi),y=quantile(absdifference[,10],phi))

#Set the bootstrap. no perturbations used. uncomment to add
area_abs_diff=test
alpha = 0.5
theta<-as.vector((1/area_abs_diff)^alpha/sum((1/area_abs_diff)^alpha))
N=ncol(FinalXtL1star)
nboot = 500
boot_area_abs_diff<-c()
for(i in 1:nboot){
  index_samp<-sample(1:N,N,replace=TRUE,prob = theta) #resamples the needed data with the given probability
  # sigma = 0 #runif(1,2.804333,4.2065)
  #  drift = 1 # runif(1,0.5,1)
  # OUyData<- OUonFloodEM(sigma=sigma,drift=drift,times=2000,endtime = 1000,Data=FinalYtL1star[index_samp])
  # new_y<-matrix(NA,nrow=nrow(FinalYtL1star),ncol=ncol(FinalYtL1star))
  # for(i in 1:ncol(FinalYtL1star)){
  #    new_y[,i]<-supsmu(x=OUyData$timevec,y=OUyData$obscurves[,i])$y #generically smooth the OU data to look more like my data
  #  }
  temp_area<-area_abs_diffFour(xData = FinalXtL1star[,index_samp],yData = FinalYtL1star[,index_samp],BasisNum = 11, predvec = Oct15CongHt$CongHt)
  boot_area_abs_diff<-c(boot_area_abs_diff,temp_area)
}
out0<-c(quantile(boot_area_abs_diff,0.9),quantile(boot_area_abs_diff,0.95),quantile(boot_area_abs_diff,0.99))
#This code makes Table 4





##Now looking at Air and Water Temperature
AirTempMat<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/MainAirTemp.csv")
AirTempMat=AirTempMat[,-1]
WaterTempMat<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/MainWaterTemp.csv")
WaterTempMat=WaterTempMat[,-1]
ai_air<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/adakislandair.csv")
ai_air=ai_air[,-1]
ai_water<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/adakislandwater.csv")
ai_water=ai_water[,-1]

kahului_air<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/kahuluiair.csv")
kahului_air=kahului_air[,-1]
kahului_water<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/kahuluiwater.csv")
kahului_water=kahului_water[,-1]

phb_air<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/prudhoebayair.csv")
phb_air=phb_air[,-1]
phb_water<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/prudhoebaywater.csv")
phb_water=phb_water[,-1]

rockport_air<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/rockportair.csv")
rockport_air=rockport_air[,-1]
#rockport_water<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/rockportwater.csv")
#rockport_water=rockport_water[,-1]

sjs_air<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/shipjohnshoalair.csv")
sjs_air=sjs_air[,-1]
sjs_water<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/shipjohnshoalwater.csv")
sjs_water=sjs_water[,-1]

library(fdaconcur)

tempDT<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/tempDT.csv")
tempDT=tempDT[,-1]

##Predicting the water temp curves that were missing. just adjust predvec.
FullTempResults<-PredictFRegressNormTestUpdatedForRUpdate(AirTempMat,WaterTempMat, ai_air,Basis = "BSpline", nBasis = 21, Plot = FALSE)
par(mfrow=c(1,1))
plot(x=as_datetime(tempDT),y=FullTempResults$PredictedResponse,type = "l", ylim = c(0,40), ylab = "Temp C.",col = "blue",main = "SJS: True Temp (red) vs Pred. WaterTemp (blue)")
lines(x=as_datetime(tempDT),y=FullTempResults$Lower,col="green")
lines(x=as_datetime(tempDT),y=FullTempResults$Upper,col="green")
lines(x=as_datetime(sjs_water$V1),y=sjs_water$V2,col = "Red")
lines(x=as_datetime(tempDT),y=i_leftout[,31])

# now we have FullTempResults$PredictedResponse and need to get the same with each left out:
i_leftout<-matrix(NA,nrow = 1000,ncol = 35)
for(i in 1:35){
  i_leftout[,i]<-PredictFRegressNormTestUpdatedForRUpdate(AirTempMat[,-i],WaterTempMat[,-i], ai_air,Basis = "BSpline", nBasis = 21, Plot = FALSE)$PredictedResponse
}

ri<-c()
squaredifference<-matrix(NA,nrow = nrow(i_leftout),ncol = ncol(i_leftout))
for(i in 1:35){
  squaredifference[,i]<-(FullTempResults$PredictedResponse - i_leftout[,i])^2
  ri[i]<-sqrt(trapz(x=1:length(squaredifference[,i]), y = squaredifference[,i])) #need
}
#Table with ri = Delta is in appendix

bSplineyhat<-function(xData,yData,predvec,BasisNum=21, lambda=10^(-1)){
  n<-nrow(xData)
  nPred<-ncol(xData) #Number of predictors (ie flood events)
  
  gaittime <- seq(1:n)
  gaitrange <- c(0,n)
  gaitfine = seq(0,n,.2)
  
  mygaitExp <- array(NA, dim = c(n,nPred,2))
  mygaitExp[1:n, ,] <- seq(1:n)
  
  for(i in 1:nPred){
    mygaitExp[,i, 1] <- xData[,i]
    mygaitExp[,i, 2] <- yData[,i]
  }
  
  gaitbasis <- create.bspline.basis(gaitrange, nbasis=BasisNum,norder=6) #original 20 leaving norder at 6.
  D2fdPar = fdPar(gaitbasis, lambda=lambda)
  gaitSmooth = smooth.basis(gaittime, mygaitExp, D2fdPar)
  betalist  = list(const=D2fdPar, cong=D2fdPar)
  gaitfd = gaitSmooth$fd
  
  names(gaitfd$fdnames) = c("Normalized time", "Event", "Height")
  gaitfd$fdnames[[3]] = c("Cong", "Cedar")
  
  congfd  = gaitfd[,1] #Predictor Stuff
  cedfd = gaitfd[,2] # Response Stuff This seems like a more detailed version of what we get from the other textbook.
  
  xfdlist   = list(const=rep(1,nPred), cong=congfd)
  
  gaitRegress= fRegress(cedfd, xfdlist, betalist)
  betaestlist = gaitRegress$betaestlist
  cedIntercept = predict(betaestlist$const$fd, gaitfine) #B0
  congCoef = predict(betaestlist$cong$fd, gaitfine) #Slope Term B1
  
  Intercept<-approx(x=1:length(cedIntercept), y=cedIntercept, n = n)$y
  Slope<-approx(x=1:length(congCoef), y=congCoef , n = n)$y
  yHat = Intercept + Slope*predvec
  out<-list(Intercept,Slope,yHat)
  names(out)<-c("B0t","B1t","yHatout")
  return(out)
}

deltaPfun<-function(xData,yData,predvec,BasisNum=21, lambda=10^(-1))
{
  FullResults<-bSplineyhat(xData=xData,yData=yData,predvec=predvec,BasisNum=BasisNum, lambda=10^(-1))$yHatout
  i_leftout<-matrix(NA,nrow = nrow(xData),ncol = ncol(xData))
  for(i in 1:ncol(i_leftout)){
    i_leftout[,i]<-bSplineyhat(xData=xData[,-i],yData=yData[,-i],predvec=predvec,BasisNum=BasisNum, lambda=10^(-1))$yHatout #calculate the fitted vallues for predictor vec with each left out
  }
  ri<-c()
  squaredifference<-matrix(NA,nrow = nrow(i_leftout),ncol = ncol(i_leftout))
  for(i in 1:ncol(i_leftout)){
    squaredifference[,i]<-(FullResults - i_leftout[,i])^2
    ri[i]<-sqrt(trapz(x=1:length(squaredifference[,i]), y = squaredifference[,i]))
  }
  out<-ri
  return(out)
}

temp_boot<-function(alpha=0,obsDP,nboot = 10,predvec){
  theta<-as.vector((1/obsDP)^alpha/sum((1/obsDP)^alpha))
  N=ncol(AirTempMat)
  boot_deltaP<-c()
  for(i in 1:nboot){
    index_samp<-sample(1:N,N,replace=TRUE,prob = theta) #resamples the needed data with the given probability
    temp_deltaP<-deltaPfun(xData = AirTempMat[,index_samp],yData = WaterTempMat[,index_samp],BasisNum = 21, predvec = predvec)
    boot_deltaP<-c(boot_deltaP,temp_deltaP)
  }
  out<-c(quantile(boot_deltaP,0.9),quantile(boot_deltaP,0.95),quantile(boot_deltaP,0.99))
  return(out)
}

#Adjust alpha and predvec based on desired station
test<-temp_boot(alpha=0,obsDP=ri,predvec = rockport_air,nboot = 10)
testp1<-temp_boot(alpha=0.1,obsDP=ri,predvec = rockport_air,nboot = 100)
testp3<-temp_boot(alpha=0.3,obsDP=ri,predvec = rockport_air,nboot = 100)
testp5<-temp_boot(alpha=0.5,obsDP=ri,predvec = rockport_air,nboot = 100)
#Results in Table 5


##AIP Air and Water
area_abs_diffBspline<-function(xData,yData,predvec,BasisNum=21, lambda=10^(-1))
{
  FullResults<-PredictFRegressNormTestUpdatedForRUpdate(PredictorMat =xData,ResponseMat=yData,predVec=predvec,Basis = "BSpline",nBasis=BasisNum, lambda=10^(-1))$PredictedResponse
  i_leftout<-matrix(NA,nrow = nrow(xData),ncol = ncol(xData))
  for(i in 1:ncol(i_leftout)){
    i_leftout[,i]<-PredictFRegressNormTestUpdatedForRUpdate(PredictorMat =xData[,-i],ResponseMat=yData[,-i],predVec=predvec,nBasis=BasisNum, Basis = "BSpline",lambda=10^(-1))$PredictedResponse #calculate the fitted vallues for predictor vec with each left out
  }
  area_measure<-c()
  phi=seq(.01,1,by=.01)
  absdifference<-matrix(NA,nrow = nrow(i_leftout),ncol = ncol(i_leftout))
  rawpercentiles<-matrix(NA,nrow = length(phi),ncol = ncol(i_leftout))
  for(i in 1:ncol(i_leftout)){
    absdifference[,i]<-abs(FullResults - i_leftout[,i])
  }
  phi=seq(.01,1,by=.01)
  for(i in 1:ncol(xData)) {
    area_measure<-c(area_measure,trapz(x=1:length(phi),y=quantile(absdifference[,i],phi)))
    rawpercentiles[,i]<-quantile(absdifference[,i],phi)
  }
  out<-list(area_measure,rawpercentiles)
  names(out)=c("areameasure","rawpercentiles")
  return(out)
}

phi=seq(.01,1,by=.01)
out<-area_abs_diffBspline(AirTempMat,WaterTempMat,ai_air)
par(mfrow=c(1,1))

#Figure 13
plot(phi,out$rawpercentiles[,1],type = "l",lwd=2,ylim = c(0,3),ylab = "Absolute Difference (\u00B0C)" ,xlab = "Percentile", main = "Absolute Difference Percentiles")
for(i in 2:35){
  lines(phi,out$rawpercentiles[,i],col=i,lwd=2)
}


OUonFloodEM<-function(sigma=1,drift=1,times,endtime,Data){
  dt<-endtime/times; #want this to be around .5 for most situations to keep it stable enough
  objects<-ncol(Data) #number of flood events I have
  normalvec<-matrix(rnorm(objects*times,sd=1),nrow=objects,ncol=times); #this is randn from matlab code
  ouvalue<-matrix(0:0,nrow=objects,ncol=times) #null matrix of 0s
  
  # Defining the cluster sizes and the true clustering structure: No clusters here so just 1
  timevec<-seq(0,endtime,length=times)
  
  # loop to create approximation of Ornstein-Uhlenbeck process	# Using Euler-Maruyama approximation
  ouvalue[1:objects,1]=normalvec[,1] #So we don't start at 0.
  for (timeint in 2:times)
  {
    ouvalue[1:objects,timeint]<- ouvalue[1:objects,timeint-1]-drift*ouvalue[1:objects,timeint-1]*dt + sigma*normalvec[1:objects,timeint]*sqrt(dt) #normalvec[1:clussizes[1],timeint] n random values from normal(0, sigma)
  }
  
  samesizeOU<-matrix(0:0,nrow=nrow(Data),ncol=objects) #makes sure it is the same size as the vector it is being applied to
  for(i in 1:objects){
    samesizeOU[,i]<-approx(x=1:ncol(ouvalue),y = ouvalue[i,], n=nrow(Data))$y
  }
  # Defining and naming the observed curves
  obscurves<-Data+samesizeOU;
  out<-list(1:nrow(Data),obscurves,samesizeOU)
  names(out)<-c("timevec","obscurves","samesizeOU")
  return(out)
}

##Start here adjusting Alpha and the predvec for AIP boot.
predvec =  kahului_air
aw_areaabsdiff=area_abs_diffBspline(AirTempMat,WaterTempMat,predvec = predvec)$areameasure #gets observed area for that predvec and all observed values
alpha = 0
theta<-as.vector((1/aw_areaabsdiff)^alpha/sum((1/aw_areaabsdiff)^alpha)) #selection probabilities
N=ncol(WaterTempMat)
nboot = 100
boot_area_abs_diff<-c()
for(i in 1:nboot){
  index_samp<-sample(1:N,N,replace=TRUE,prob = theta) #resamples the needed data with the given probability
  # sigma = runif(1,4.4,6.6)
  # drift = runif(1,0.5,1)
  # OUyData<- OUonFloodEM(sigma=sigma,drift=drift,times=1000,endtime = 500,Data=WaterTempMat[index_samp])
  # newY<-matrix(NA,nrow=nrow(WaterTempMat),ncol=ncol(WaterTempMat)) #perturbing the Y
  # for(i in 1:ncol(WaterTempMat)){
  #   newY[,i]<-supsmu(x=OUyData$timevec,y=OUyData$obscurves[,i])$y #generically smooth the OU data to look more like my data
  # }
  temp_area<-area_abs_diffBspline(xData = AirTempMat[,index_samp],yData = WaterTempMat[,index_samp],BasisNum = 21, predvec = predvec)$areameasure
  boot_area_abs_diff<-c(boot_area_abs_diff,temp_area)
}
boot_area_abs_diff1=boot_area_abs_diff #needed this to do in parts bc R crashed with Nboot = 100
boot_area_abs_diff<-c(boot_area_abs_diff1,boot_area_abs_diff)
out0<-c(quantile(boot_area_abs_diff,0.9),quantile(boot_area_abs_diff,0.95),quantile(boot_area_abs_diff,0.99))
c(quantile(boot_area_abs_diff,0.9),quantile(boot_area_abs_diff,0.95),quantile(boot_area_abs_diff,0.99))
#These result in Table 6








##Simulation Plot Code and Map:

###Power and Pvalue plots for AIP and Delta simulations

library(ggplot2)
##plots
power0N10<-c(.39,.4,.13,.04,.02,.06,.15,.36,.32,.37)
se_p0N10<-c(.049,.0492,.0338,.0197,.0141,.0239,.0359,.0482,.0469,.0485)
power01N10<-c(.39,.45,.18,.03,.01,.06,.15,.39,.39,.39)
se_p01N10<-c(.049,.05,.0386,.0171,.01,.0239,.0359,.049,.049,.049)
power03N10<-c(.48,.46,.17,.03,.01,.05,.17,.39,.45,.52)
se_p03N10<-c(.0502,.0501,.0378,.0171,.01,.0219,.0378,.049,.05,.0502)
power05N10<-c(.56,.55,.17,.05,.02,.04,.18,.45,.56,.62)
se_p05N10<-c(.0499,.05,.0378,.0219,.0141,.0197,.0386,.05,.0499,.0488)

alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN10<-c(power0N10,power01N10,power03N10,power05N10)
seN10<-c(se_p0N10,se_p01N10,se_p03N10,se_p05N10)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N10<-data.frame(cbind(lambda,alpha,powerN10,seN10))
p<- ggplot(N10, aes(x=lambda, y=powerN10, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN10-seN10, ymax=powerN10+seN10), width=.1)

one<-p+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


#N20
power0N20<-c(.92,.9,.65,.21,.02,.14,.68,.89,.89,.89)
se_p0N20<-c(.0273,.0302,.0479,.0409,.0141,.0349,0.0469,.0314,.0314,.0314)
power01N20<-c(.94,.9,.67,.17,.02,.12,.66,.87,.91,.91)
se_p01N20<-c(.0239,.0302,.0473,.0378,.0141,.0327,.0476,.0338,.0288,.0288)
power03N20<-c(.94,.92,.65,.21,.02,.12,.63,.91,.92,.93)
se_p03N20<-c(.0239,.0273,.0479,.0409,.0141,.0327,.0485,.0288,.0273,.0256)
power05N20<-c(.96,.92,.67,.19,.02,.13,.63,.9,.93,.95)
se_p05N20<-c(.0197,.0273,.0473,.0394,.0141,.0338,.0485,.0302,.0256,.0219)

alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN20<-c(power0N20,power01N20,power03N20,power05N20)
seN20<-c(se_p0N20,se_p01N20,se_p03N20,se_p05N20)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N20<-data.frame(cbind(lambda,alpha,powerN20,seN20))
p2<- ggplot(N20, aes(x=lambda, y=powerN20, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN20-seN20, ymax=powerN20+seN20), width=.1)

two<-p2+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


#N50
power0N50<-c(.99,.99,.82,.22,.05,.28,.82,.98,.99,1)
se_p0N50<-c(.01,.01,.0386,.0416,.0219,.0451,.0386,.0141,.01,0)
power01N50<-c(.99,.98,.81,.19,.05,.27,.8,.98,.99,1)
se_p01N50<-c(.01,.0141,.0394,.0394,.0219,.0446,.0402,.0141,.01,0)
power03N50<-c(.99,.99,.78,.21,.05,.26,.79,.96,.99,1)
se_p03N50<-c(.01,.01,.0416,.0409,.0219,.0441,.0409,.0197,.01,0)
power05N50<-c(.99,.99,.76,.23,.05,.25,.78,.96,.98,1)
se_p05N50<-c(.01,.01,.0429,.0423,.0219,.0435,.0416,.0197,.0141,0)


alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN50<-c(power0N50,power01N50,power03N50,power05N50)
seN50<-c(se_p0N50,se_p01N50,se_p03N50,se_p05N50)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N50<-data.frame(cbind(lambda,alpha,powerN50,seN50))
p3<- ggplot(N50, aes(x=lambda, y=powerN50, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN50-seN50, ymax=powerN50+seN50), width=.1)

three<-p3+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


#N100
power0N100<-c(1,.98,.81,.31,.08,.31,.82,.99,1,.99)
se_p0N100<-c(0,.0141,.0394,.0465,.0273,.0465,.0386,.01,0,.01)
power01N100<-c(1,.98,.82,.3,.08,.3,.82,.96,1,.99)
se_p01N100<-c(0,.0141,.0386,.0461,.0273,.0461,.0386,.0197,0,.01)
power03N100<-c(1,.98,.81,.29,.08,.27,.82,.96,1,.99)
se_p03N100<-c(0,.0141,.0394,.0456,.0273,.0446,.0386,.0197,0,.01)
power05N100<-c(1,.98,.8,.3,.07,.27,.83,.96,1,.99)
se_p05N100<-c(0,.0141,.0402,.0461,.0256,.0446,.0378,.0197,0,.01)

alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN100<-c(power0N100,power01N100,power03N100,power05N100)
seN100<-c(se_p0N100,se_p01N100,se_p03N100,se_p05N100)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N100<-data.frame(cbind(lambda,alpha,powerN100,seN100))
p4<- ggplot(N100, aes(x=lambda, y=powerN100, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN100-seN100, ymax=powerN100+seN100), width=.1)

four<-p4+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


library(gridExtra)
library(ggpubr)
ggarrange(one, two,three,  four,
          labels = c("N = 10","N = 20", "N = 50", "N = 100"),
          ncol = 2, nrow = 2)

##Power for delta P
pval0N10<-1-c(.925,.921,.877,.705,.549,.684,.868,.916,.925,.93)
se_pval0N10<-c(.0045,.0058,.0109,.0208,.0246,.0225,.0102,.007,.0049,.0041)
pval01N10<-1-c(.927,.921,.879,.697,.544,.683,.867,.917,.929,.931)
se_pval01N10<-c(.0045,.0058,.0107,.0213,.0242,.0225,.0103,.007,.0047,.0041)
pval03N10<-1-c(.93,.925,.876,.688,.54,.675,.864,.918,.933,.938)
se_pval03N10<-c(.0047,.006,.0106,.0213,.0241,.0227,.0105,.0076,.0046,.0038)
pval05N10<-1-c(.939,.928,.873,.684,.532,.669,.862,.921,.941,.947)
se_pval05N10<-c(.0043,.0064,.0108,.0215,.0233,.0227,.0112,.0078,.0047,.0039)

#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN10<-c(pval0N10,pval01N10,pval03N10,pval05N10)
seN10<-c(se_pval0N10,se_pval01N10,se_pval03N10,se_pval05N10)
lambda<-rep(c( .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N10<-data.frame(cbind(lambda,alpha,pvalN10,seN10))
library(ggplot2)
p<- ggplot(N10, aes(x=lambda, y=pvalN10, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN10-seN10, ymax=pvalN10+seN10), width=.1)
#print(p)
# Finished line plot
one<-p+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = .95, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

#N=20
pval0N20<-1-c(.971,.972,.938,.794,.466,.795,.94,.971,.968,.971)
se_pval0N20<-c(.0025,.0027,.0086,.019,.025,.0167,.0083,.0019,.0035,.0022)
pval01N20<-1-c(.973,.973,.935,.793,.464,.79,.937,.972,.969,.973)
se_pval01N20<-c(.0025,.0026,.0086,.0189,.0248,.0167,.0085,.0019,.0034,.0019)
pval03N20<-1-c(.976,.973,.933,.788,.464,.782,.936,.974,.972,.977)
se_pval03N20<-c(.0025,.0029,.0091,.019,.0244,.0166,.0086,.002,.0034,.0021)
pval05N20<-1-c(.979,.976,.933,.783,.464,.776,.935,.977,.975,.981)
se_pval05N20<-c(.0025,.0031,.0097,.019,.024,.0168,.0092,.0019,.0036,.0019)

#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN20<-c(pval0N20,pval01N20,pval03N20,pval05N20)
seN20<-c(se_pval0N20,se_pval01N20,se_pval03N20,se_pval05N20)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N20<-data.frame(cbind(lambda,alpha,pvalN20,seN20))
library(ggplot2)
p2<- ggplot(N20, aes(x=lambda, y=pvalN20, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN20-seN20, ymax=pvalN20+seN20), width=.1)
#print(p)
# Finished line plot
two<-p2+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

#N=50
pval0N50<-1-c(.99,.988,.96,.819,.52,.796,.963,.987,.989,.991)
se_pval0N50<-c(.0006,.0009,.0063,.0186,.0274,.0213,.0055,.0018,.0018,.0004)
pval01N50<-1-c(.991,.988,.958,.818,.52,.793,.962,.987,.989,.992)
se_pval01N50<-c(.0006,.0011,.0067,.0186,.0272,.0212,.0056,.0019,.002,.0003)
pval03N50<-1-c(.992,.989,.957,.815,.52,.789,.959,.987,.991,.993)
se_pval03N50<-c(.0006,.0011,.0068,.0184,.027,.0212,.0061,.0022,.0019,.0003)
pval05N50<-1-c(.994,.99,.955,.812,.52,.785,.957,.988,.992,.995)
se_pval05N50<-c(.0006,.0011,.0072,.0185,.0267,.0211,.0062,.0025,.002,.0002)
#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN50<-c(pval0N50,pval01N50,pval03N50,pval05N50)
seN50<-c(se_pval0N50,se_pval01N50,se_pval03N50,se_pval05N50)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N50<-data.frame(cbind(lambda,alpha,pvalN50,seN50))
library(ggplot2)
p3<- ggplot(N50, aes(x=lambda, y=pvalN50, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN50-seN50, ymax=pvalN50+seN50), width=.1)
#print(p)
# Finished line plot
three<-p3+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

#N = 100
pval0N100<-1-c(.995,.993,.966,.791,.525,.806,.97,.993,.995,.994)
se_pval0N100<-c(.0002,.0014,.007,.0223,.0305,.0228,.0052,.0009,.0003,.0008)
pval01N100<-1-c(.996,.993,.965,.79,.526,.805,.968,.993,.995,.995)
se_pval01N100<-c(.0002,.0014,.0072,.0221,.0305,.0227,.0054,.001,.0004,.001)
pval03N100<-1-c(.996,.992,.963,.787,.527,.802,.967,.993,.996,.995)
se_pval03N100<-c(.0002,.0016,.0076,.0221,.0301,.0228,.0056,.0012,.0004,.0011)
pval05N100<-1-c(.997,.993,.961,.785,.528,.799,.965,.993,.997,.996)
se_pval05N100<-c(.0002,.0016,.008,.0221,.0299,.0229,.006,.0013,.0005,.0011)

#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN100<-c(pval0N100,pval01N100,pval03N100,pval05N100)
seN100<-c(se_pval0N100,se_pval01N100,se_pval03N100,se_pval05N100)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N100<-data.frame(cbind(lambda,alpha,pvalN100,seN100))
library(ggplot2)
p4<- ggplot(N100, aes(x=lambda, y=pvalN100, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN100-seN100, ymax=pvalN100+seN100), width=.1)
#print(p)
# Finished line plot
four<-p4+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

ggarrange(one, two,three,  four,
          labels = c("N = 10","N = 20", "N = 50", "N = 100"),
          ncol = 2, nrow = 2)


##Plotting the new measure simulation results
power0N10<-c(.52,.47,.26,.02,.01,.03,.22,.52,.47,.51)
se_p0N10<-c(.0502,.0502,.0441,.0141,.01,.0171,.0416,.0502,.0502,.0502)
power01N10<-c(.54,.45,.24,.03,.01,0,.18,.56,.52,.54)
se_p01N10<-c(.0501,.05,.0429,.0171,.01,0,.0386,.0499,.0502,.0501)
power03N10<-c(.59,.53,.28,.02,.01,0,.17,.56,.58,.59)
se_p03N10<-c(.0494,.0502,.0451,.0141,.01,0,.0378,.0499,.0496,.0494)
power05N10<-c(.7,.62,.33,.03,.01,.02,.19,.68,.69,.69)
se_p05N10<-c(.0461,.0488,.0473,.0171,.01,.0141,.0394,.0469,.0465,.0465)

alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN10<-c(power0N10,power01N10,power03N10,power05N10)
seN10<-c(se_p0N10,se_p01N10,se_p03N10,se_p05N10)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N10<-data.frame(cbind(lambda,alpha,powerN10,seN10))
p<- ggplot(N10, aes(x=lambda, y=powerN10, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN10-seN10, ymax=powerN10+seN10), width=.1)

one<-p+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


#N20
power0N20<-c(.95,.93,.71,.15,.02,.21,.68,.92,.97,.98)
se_p0N20<-c(.0219,.0256,.0456,.0359,.0141,.0409,.0469,.0273,.0171,.0141)
power01N20<-c(.95,.94,.72,.16,.01,.21,.67,.92,.98,.99)
se_p01N20<-c(.0219,.0239,.0451,.0368,.01,.0409,.0473,.0273,.0141,.01)
power03N20<-c(.96,.95,.71,.14,0,.2,.65,.91,.98,.99)
se_p03N20<-c(.0197,.0219,.0456,.0349,0,.0402,.0479,.0288,.0141,.01)
power05N20<-c(.97,.94,.67,.14,0,.21,.66,.91,.98,.99)
se_p05N20<-c(.0171,.0239,.0473,.0349,0,.0409,.0476,.0288,.0141,.01)

alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN20<-c(power0N20,power01N20,power03N20,power05N20)
seN20<-c(se_p0N20,se_p01N20,se_p03N20,se_p05N20)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N20<-data.frame(cbind(lambda,alpha,powerN20,seN20))
p2<- ggplot(N20, aes(x=lambda, y=powerN20, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN20-seN20, ymax=powerN20+seN20), width=.1)

two<-p2+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


#N50
power0N50<-c(1,.99,.87,.35,.05,.38,.83,.98,1,1)
se_p0N50<-c(0,.01,.0338,.0479,.0219,.0488,.0378,.0141,0,0)
power01N50<-c(1,.99,.87,.33,.04,.35,.83,.98,1,1)
se_p01N50<-c(0,.01,.0338,.0473,.0197,.0479,.0378,.0141,0,0)
power03N50<-c(1,.99,.87,.33,.04,.35,.8,.98,1,1)
se_p03N50<-c(0,.01,.0338,.0473,.0197,.0479,.0402,.0141,0,0)
power05N50<-c(1,.98,.87,.32,.04,.35,.8,.98,1,1)
se_p05N50<-c(0,.0141,.0338,.0469,.0197,.0479,.0394,.0141,0,0)


alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN50<-c(power0N50,power01N50,power03N50,power05N50)
seN50<-c(se_p0N50,se_p01N50,se_p03N50,se_p05N50)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N50<-data.frame(cbind(lambda,alpha,powerN50,seN50))
p3<- ggplot(N50, aes(x=lambda, y=powerN50, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN50-seN50, ymax=powerN50+seN50), width=.1)

three<-p3+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


#N100
power0N100<-c(1,.99,.92,.37,.05,.32,.85,1,1,1)
se_p0N100<-c(0,.01,.0273,.0485,.0219,.0469,.0359,0,0,0)
power01N100<-c(1,.99,.92,.38,.05,.3,.85,1,1,1)
se_p01N100<-c(0,.01,.0273,.0488,.0219,.0461,.0359,0,0,0)
power03N100<-c(1,.99,.91,.38,.05,.3,.84,.99,1,1)
se_p03N100<-c(0,.01,.0288,.0488,.0219,.0461,.0368,0.01,0,0)
power05N100<-c(1,.99,.91,.36,.04,.29,.84,.99,1,1)
se_p05N100<-c(0,.01,.0288,.0482,.0197,.0456,.0368,0.01,0,0)

alpha<-rep(c(0,.1,0.3,0.5),each= 10)
powerN100<-c(power0N100,power01N100,power03N100,power05N100)
seN100<-c(se_p0N100,se_p01N100,se_p03N100,se_p05N100)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N100<-data.frame(cbind(lambda,alpha,powerN100,seN100))
p4<- ggplot(N100, aes(x=lambda, y=powerN100, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=powerN100-seN100, ymax=powerN100+seN100), width=.1)

four<-p4+labs(title=" ", x=expression(lambda), y = "Power")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "Power = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))


library(gridExtra)
library(ggpubr)
ggarrange(one, two,three,  four,
          labels = c("N = 10","N = 20", "N = 50", "N = 100"),
          ncol = 2, nrow = 2)

##pval for new measure
pval0N10<-1-c(.937,.935,.905,.754,.58,.74,.888,.935,.939,.936)
se_pval0N10<-c(.0049,.0038,.0064,.0151,.0246,.0168,.0077,.0049,.003,.0051)
pval01N10<-1-c(.94,.934,.903,.749,.578,.735,.89,.937,.942,.938)
se_pval01N10<-c(.0043,.004,.0066,.015,.0243,.0167,.0074,.0047,.0028,.0046)
pval03N10<-1-c(.945,.938,.903,.741,.575,.728,.889,.94,.949,.944)
se_pval03N10<-c(.0043,.0042,.0068,.0154,.0239,.0169,.0081,.0046,.0026,.0044)
pval05N10<-1-c(.951,.944,.902,.733,.572,.723,.886,.944,.956,.951)
se_pval05N10<-c(.0042,.0041,.0075,.0156,.0234,.0166,.0082,.0047,.0025,.004)

#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN10<-c(pval0N10,pval01N10,pval03N10,pval05N10)
seN10<-c(se_pval0N10,se_pval01N10,se_pval03N10,se_pval05N10)
lambda<-rep(c( .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N10<-data.frame(cbind(lambda,alpha,pvalN10,seN10))
library(ggplot2)
p<- ggplot(N10, aes(x=lambda, y=pvalN10, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN10-seN10, ymax=pvalN10+seN10), width=.1)
#print(p)
# Finished line plot
one<-p+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = .95, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

#N=20
pval0N20<-1-c(.978,.975,.952,.804,.487,.769,.949,.974,.976,.979)
se_pval0N20<-c(.0013,.0021,.0052,.0172,.0264,.022,.0056,.0019,.0021,.001)
pval01N20<-1-c(.979,.975,.951,.801,.497,.769,.949,.975,.977,.98)
se_pval01N20<-c(.0013,.0024,.0054,.0174,.026,.022,.0056,.0017,.002,.0009)
pval03N20<-1-c(.981,.976,.95,.797,.485,.765,.946,.977,.979,.983)
se_pval03N20<-c(.0012,.0029,.0057,.0174,.0255,.0219,.0062,.0019,.0017,.0008)
pval05N20<-1-c(.984,.978,.949,.793,.485,.763,.945,.978,.983,.987)
se_pval05N20<-c(.0012,.0032,.0061,.0175,.025,.0218,.0066,.0021,.0019,.0007)

#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN20<-c(pval0N20,pval01N20,pval03N20,pval05N20)
seN20<-c(se_pval0N20,se_pval01N20,se_pval03N20,se_pval05N20)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N20<-data.frame(cbind(lambda,alpha,pvalN20,seN20))
library(ggplot2)
p2<- ggplot(N20, aes(x=lambda, y=pvalN20, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN20-seN20, ymax=pvalN20+seN20), width=.1)
#print(p)
# Finished line plot
two<-p2+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

#N=50
pval0N50<-1-c(.992,.989,.975,.828,.498,.811,.965,.99,.991,.991)
se_pval0N50<-c(.0003,.0015,.006,.0196,.0285,.0224,.0075,.0008,.0004,.0002)
pval01N50<-1-c(.992,.989,.974,.828,.499,.809,.964,.99,.992,.992)
se_pval01N50<-c(.0003,.0016,.0063,.0195,.0283,.0224,.0075,.0008,.0004,.0002)
pval03N50<-1-c(.993,.99,.973,.827,.502,.808,.961,.991,.993,.994)
se_pval03N50<-c(.0003,.0017,.0064,.0192,.0279,.0221,.008,.0008,.0004,.0002)
pval05N50<-1-c(.995,.991,.973,.825,.505,.806,.961,.992,.995,.996)
se_pval05N50<-c(.0003,.0018,.0069,.0192,.0276,.0221,.0084,.001,.0004,.0002)
#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN50<-c(pval0N50,pval01N50,pval03N50,pval05N50)
seN50<-c(se_pval0N50,se_pval01N50,se_pval03N50,se_pval05N50)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N50<-data.frame(cbind(lambda,alpha,pvalN50,seN50))
library(ggplot2)
p3<- ggplot(N50, aes(x=lambda, y=pvalN50, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN50-seN50, ymax=pvalN50+seN50), width=.1)
#print(p)
# Finished line plot
three<-p3+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

#N = 100
pval0N100<-1-c(.996,.994,.984,.833,.495,.817,.973,.995,.996,.996)
se_pval0N100<-c(.0002,.0007,.0027,.019,.0313,.0204,.0048,.0003,.0001,.0001)
pval01N100<-1-c(.996,.994,.983,.832,.496,.816,.972,.995,.996,.996)
se_pval01N100<-c(.0002,.0008,.0029,.0189,.0311,.0203,.005,.0004,.0001,.0001)
pval03N100<-1-c(.997,.995,.982,.829,.499,.814,.971,.995,.997,.997)
se_pval03N100<-c(.0001,.0009,.003,.0189,.0308,.0201,.0051,.0005,.0002,.0001)
pval05N100<-1-c(.997,.995,.981,.827,.5,.812,.97,.996,.997,.998)
se_pval05N100<-c(.0002,.0011,.0032,.0188,.0304,.0202,.0054,.0007,.0002,.0001)

#ggplot2
alpha<-rep(c(0,.1,0.3,0.5),each= 10)
pvalN100<-c(pval0N100,pval01N100,pval03N100,pval05N100)
seN100<-c(se_pval0N100,se_pval01N100,se_pval03N100,se_pval05N100)
lambda<-rep(c(.25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 1.75,2),4)

N100<-data.frame(cbind(lambda,alpha,pvalN100,seN100))
library(ggplot2)
p4<- ggplot(N100, aes(x=lambda, y=pvalN100, group=alpha, color=as.character(alpha))) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=pvalN100-seN100, ymax=pvalN100+seN100), width=.1)
#print(p)
# Finished line plot
four<-p4+labs(title=" ", x=expression(lambda), y = "Average p-value")+labs(colour="Alpha")+
  theme_classic() + scale_x_continuous(breaks=c( .5, .75, 1, 1.25, 1.5, 1.75,2))+
  geom_hline(yintercept=0.05, linetype='dashed', col = 'black')+
  annotate("text", x = 1.75, y = 0.05, label = "p-value = 0.05", vjust = -0.5)+
  scale_color_manual(values=c('dodgerblue4',"darkorange2","firebrick4","seagreen2"))

ggarrange(one, two,three,  four,
          labels = c("N = 10","N = 20", "N = 50", "N = 100"),
          ncol = 2, nrow = 2)



library(mapproj)
map(database= "world", ylim=c(20,70), xlim=c(-190,-30), col="grey80", fill=TRUE, projection="gilbert", orientation= c(90,0,225))
lon <- c(-68.204 ,-71.050 , -70.246 ,-73.181 ,-76.039 ,-74.418 ,
         -75.548 ,-76.671 ,-77.786 ,-79.924 ,-80.903 ,-81.465 ,
         -80.034 ,-81.807,-82.553 , -82.832 ,-85.880 ,-89.326 ,
         -91.338 ,-93.343 , -118.500 ,-120.754 ,-124.498 ,
         -124.105 ,-123.441 ,-177.361 ,-145.752 ,
         -131.625,-135.327 ,-124.184, -122.039 ,-162.327,-157.79,-97.215,
         -164.067  )  #longitude vector. needs -
lat <- c(44.392 ,42.355 , 43.656 ,41.174 ,38.220 , 39.357 ,
         35.796 ,34.717 ,34.213,32.781 ,32.035 ,30.675 ,
         26.613 ,26.132,27.858 ,27.978 ,30.213 ,30.326 ,
         29.450 , 29.768 ,34.008 ,35.169 ,42.739 ,
         46.904 , 48.125  ,28.215  ,60.558 ,
         55.331 ,59.450 ,41.746 ,38.056,55.062, 21.433, 26.061,
         67.575
)  #latitude vector
coord <- mapproject(lon, lat, proj="gilbert", orientation=c(90, 0, 225))  #convert points to projected lat/long
points(coord, pch=20, cex=1.2, col=c(rep("red",34),"blue"))  #plot converted points


