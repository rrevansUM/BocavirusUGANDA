#####################################################################################
## Created by Josh Petrie on Jan 24, 2016
## Using estimated parameter output from mcmc, this script creates output tables of: 
##      1. mcmc acceptance rate, parameter means and 95% CI
##      2. Serial interval mean and 95% CI
##      3. Cumulative risks of community and household infection with 95% CI for all 
##         combinations of patient characteristics (e.g. <9 years, high risk, unvaccinated)
#####################################################################################


#write acceptance rate
#optimal acceptance rate: ~45-55% in single dimension
accept <-round(colMeans(a$accept,na.rm=T),2)       
# estimated mean of parameters
mu <-round(c(mean(a$para1),mean(a$para2),mean(a$para3), mean(a$para4), 
             mean(a$para5),mean(a$para6), mean(a$para7), mean(a$para8)),6)
# estimated 95% CI of parameters
ci <-round(c(quantile(a$para1,c(0.025,0.975)),quantile(a$para2,c(0.025,0.975)),
             quantile(a$para3,c(0.025,0.975)), quantile(a$para4,c(0.025,0.975)), 
             quantile(a$para5,c(0.025,0.975)),  quantile(a$para6,c(0.025,0.975)), 
             quantile(a$para7,c(0.025,0.975)),quantile(a$para8,c(0.025,0.975))),6)


#put acceptance rate, mean, and cis in a data frame
param.output <- data.frame(
  accept = accept, 
  mean = mu, 
  ci1 = ci[c(1,3,5,7,9,11,13,15)], 
  ci2 = ci[c(2,4,6,8,10,12,14,16)]
  )

#SIfn calculates the serial interval from a given set of parameters
SIfn<-function(params,HHt) {
  SI.FN <- 
    exp(-(HHt/params$para1)^params$para2)-exp(-((HHt+1)/params$para1)^params$para2)
  
  SI <- sum(HHt*SI.FN)/sum(SI.FN)
  SI
}

#SI.apply applies SIfn to parameters generated from each mcmc iteration
SI.apply<-function(params,HHt) {
  SIvec<-vector(length=length(params[[1]]))
  temp<-data.frame(para1=params[[1]],para2=params[[2]])
  for(i in 1:length(temp[,1])){
    SIvec[i]<-SIfn(temp[i,],HHt)
  }
  SIvec
}

#calculate the mean serial interval
muSI <-round(c(mean(SI.apply(a,seq(0,14,1)))),6)
#calculate the 95% CI for the serial interval
ciSI <-round(c(quantile(SI.apply(a,seq(0,14,1)),c(0.025,0.975))),6)
#put mean serial interval and 95% CI in a data frame
SI.output <- data.frame(SI=muSI,ciSIlwr=ciSI[c(1)],ciSIupr=ciSI[c(2)])

#HHRISK calculates the total household risk of infection from a given set of parameters
HHRISK<-function(params,HHt,AGE9to17,AGEgte18,HighRisk,Vax){
  SI.FN <- 
    exp(-(HHt/params$para1)^params$para2)-exp(-((HHt+1)/params$para1)^params$para2)
  
  SI <- sum(HHt*SI.FN)/sum(SI.FN)
  
  HH.HAZ <- 
    (exp(-(HHt/params$para1)^params$para2)-exp(-((HHt+1)/params$para1)^params$para2))*
    params$para3*exp(params$para5*AGE9to17+params$para6*AGEgte18+params$para7*
                       HighRisk+params$para8*Vax)
  eHH.HAZ <- exp(-HH.HAZ)
  
  HHl=loess(1-eHH.HAZ ~ HHt, span=.35)
  HHf<-function(x) predict(HHl, newdata=x)
  HHi=integrate(HHf,0,14)
  
  HHi$value
  
}
#HHRISK.apply applies HHRISK to parameters generated from each mcmc iteration
HHRISK.apply<-function(params,HHt,AGE9to17,AGEgte18,HighRisk,Vax){
  HouseRisk<-vector(length=length(params[[1]]))
  temp<-data.frame(para1=params[[1]],para2=params[[2]],para3=params[[3]],para4=params[[4]],
                   para5=params[[5]],para6=params[[6]],para7=params[[7]],para8=params[[8]])
  for(i in 1:length(temp[,1])){
    HouseRisk[i]<-HHRISK(temp[i,],HHt,AGE9to17,AGEgte18,HighRisk,Vax)
  }
  HouseRisk
}

#calculate mean household risk of infection for each combination of subject characteristics
muHH <-round(
       c(mean(HHRISK.apply(a,seq(0,14,1),0,0,0,0)),
         mean(HHRISK.apply(a,seq(0,14,1),0,0,1,0)),
         mean(HHRISK.apply(a,seq(0,14,1),0,0,0,1)),
         mean(HHRISK.apply(a,seq(0,14,1),0,0,1,1)),
         mean(HHRISK.apply(a,seq(0,14,1),1,0,0,0)),
         mean(HHRISK.apply(a,seq(0,14,1),1,0,1,0)),
         mean(HHRISK.apply(a,seq(0,14,1),1,0,0,1)),
         mean(HHRISK.apply(a,seq(0,14,1),1,0,1,1)),
         mean(HHRISK.apply(a,seq(0,14,1),0,1,0,0)),
         mean(HHRISK.apply(a,seq(0,14,1),0,1,1,0)),
         mean(HHRISK.apply(a,seq(0,14,1),0,1,0,1)),
         mean(HHRISK.apply(a,seq(0,14,1),0,1,1,1))),
       6)
#calculate 95% CI for household risk of infection for each combination of 
#subject characteristics
ciHH <-round(c(quantile(HHRISK.apply(a,seq(0,14,1),0,0,0,0),c(0.025,0.975)),
               quantile(HHRISK.apply(a,seq(0,14,1),0,0,1,0),c(0.025,0.975)),
               quantile(HHRISK.apply(a,seq(0,14,1),0,0,0,1),c(0.025,0.975)), 
               quantile(HHRISK.apply(a,seq(0,14,1),0,0,1,1),c(0.025,0.975)), 
               quantile(HHRISK.apply(a,seq(0,14,1),1,0,0,0),c(0.025,0.975)),  
               quantile(HHRISK.apply(a,seq(0,14,1),1,0,1,0),c(0.025,0.975)), 
               quantile(HHRISK.apply(a,seq(0,14,1),1,0,0,1),c(0.025,0.975)),
               quantile(HHRISK.apply(a,seq(0,14,1),1,0,1,1),c(0.025,0.975)),
               quantile(HHRISK.apply(a,seq(0,14,1),0,1,0,0),c(0.025,0.975)),
               quantile(HHRISK.apply(a,seq(0,14,1),0,1,1,0),c(0.025,0.975)),
               quantile(HHRISK.apply(a,seq(0,14,1),0,1,0,1),c(0.025,0.975)),
               quantile(HHRISK.apply(a,seq(0,14,1),0,1,1,1),c(0.025,0.975))),6)

#CRISK calculates the total community risk of infection from a given set of parameters
CRISK<-function(params,com,Ct,AGE9to17,AGEgte18,HighRisk,Vax){
  C.HAZ<-com$COM*params$para4*(exp(params$para5*AGE9to17+params$para6*AGEgte18+
                                     params$para7*HighRisk+params$para8*Vax))
  eC.HAZ <- exp(-C.HAZ)
  
  Cl=loess(1-eC.HAZ ~ Ct, span=.25)
  Cf<-function(x) predict(Cl, newdata=x)
  Ci=integrate(Cf,0,104)
  
  Ci$value
}
#CRISK.apply applies CRISK to parameters generated from each mcmc iteration
CRISK.apply<-function(params,com,Ct,AGE9to17,AGEgte18,HighRisk,Vax){
  ComRisk<-vector(length=length(params[[1]]))
  temp<-data.frame(para1=params[[1]],para2=params[[2]],para3=params[[3]],para4=params[[4]],
                   para5=params[[5]],para6=params[[6]],para7=params[[7]],para8=params[[8]])
  for(i in 1:length(temp[,1])){
    ComRisk[i]<-CRISK(temp[i,],com,Ct,AGE9to17,AGEgte18,HighRisk,Vax)
  }
  ComRisk
}
#calculate mean community risk of infection for each combination of subject characteristics
muC<-round(
  c(mean(CRISK.apply(a,com,seq(0,104,1),0,0,0,0)),
    mean(CRISK.apply(a,com,seq(0,104,1),0,0,1,0)),
    mean(CRISK.apply(a,com,seq(0,104,1),0,0,0,1)),
    mean(CRISK.apply(a,com,seq(0,104,1),0,0,1,1)),
    mean(CRISK.apply(a,com,seq(0,104,1),1,0,0,0)),
    mean(CRISK.apply(a,com,seq(0,104,1),1,0,1,0)),
    mean(CRISK.apply(a,com,seq(0,104,1),1,0,0,1)),
    mean(CRISK.apply(a,com,seq(0,104,1),1,0,1,1)),
    mean(CRISK.apply(a,com,seq(0,104,1),0,1,0,0)),
    mean(CRISK.apply(a,com,seq(0,104,1),0,1,1,0)),
    mean(CRISK.apply(a,com,seq(0,104,1),0,1,0,1)),
    mean(CRISK.apply(a,com,seq(0,104,1),0,1,1,1))),6)
#calculate 95% CI for community risk of infection for each combination of 
#subject characteristics
ciC<-round(
  c(quantile(CRISK.apply(a,com,seq(0,104,1),0,0,0,0),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),0,0,1,0),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),0,0,0,1),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),0,0,1,1),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),1,0,0,0),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),1,0,1,0),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),1,0,0,1),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),1,0,1,1),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),0,1,0,0),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),0,1,1,0),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),0,1,0,1),c(0.025,0.975)),
    quantile(CRISK.apply(a,com,seq(0,104,1),0,1,1,1),c(0.025,0.975))),6)
#labels for each combination of subject characteristics
labels<-
  c("<9,HighRisk -,Vax -", "<9,HighRisk +,Vax -", "<9,HighRisk -,Vax +", 
    "<9,HighRisk +,Vax +", "9-17,HighRisk -,Vax -", "9-17,HighRisk +,Vax -", 
    "9-17,HighRisk -,Vax +", "9-17,HighRisk +,Vax +", "18+,HighRisk -,Vax -", 
    "18+,HighRisk +,Vax -", "18+,HighRisk -,Vax +", "18+,HighRisk +,Vax +")
#put mean and 95% community and household risks in a data frame
risk.output <- data.frame(labels=labels, HHRisk=muHH, 
                           ciHHlwr=ciHH[c(1,3,5,7,9,11,13,15,17,19,21,23)], 
                           ciHHupr=ciHH[c(2,4,6,8,10,12,14,16,18,20,22,24)],
                           CRisk=muC,
                           ciClwr=ciC[c(1,3,5,7,9,11,13,15,17,19,21,23)], 
                           ciCupr=ciC[c(2,4,6,8,10,12,14,16,18,20,22,24)])
#write output
outdir<-"C:/FILEPATH/"

write.csv(param.output, paste(outdir, Sys.Date(),"mcmc_param.output.csv", sep="" ))
write.csv(SI.output, paste(outdir, Sys.Date(),"mcmc_SI.output.csv", sep="" ))
write.csv(risk.output, paste(outdir, Sys.Date(),"mcmc_risk.output.csv", sep="" ))

# plot of posterior

pdf(paste(outdir,Sys.Date(),"mcmcPostPlot.pdf", sep=""),width=24,height=36)
layout(matrix(1:8,ncol=2,byrow=T))
par(mar=c(2,4,1,1))
plot(1:iter,a$para1,type="l",ylab="", main="para1")
plot(1:iter,a$para2,type="l",ylab="", main="para2")
plot(1:iter,a$para3,type="l",ylab="", main="para3")
plot(1:iter,a$para4,type="l",ylab="", main="para4")
plot(1:iter,a$para5,type="l",ylab="", main="para5")
plot(1:iter,a$para6,type="l",ylab="", main="para6")
plot(1:iter,a$para7,type="l",ylab="", main="para7")
plot(1:iter,a$para8,type="l",ylab="", main="para8")

dev.off()