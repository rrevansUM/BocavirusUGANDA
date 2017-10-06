#####################################################################################
## Created by Sophia Ng on Aug 31, 2015 ## Modified by Josh Petrie on Nov 2, 2015
## This script estimates model parameters from the data formatted in 
## THmodel_premcmc_commented.R using mcmc
##
## Hazards of infection from both the community and household contribute to the 
## likelihood.
##
## The hazard of infection at time t from the community for each cohort member, j, is
## modeled as:
## C.HAZ.SUM.FN = com*para4*exp(para5*AGECAT2+para6*AGECAT3+para7*HR+para8*vax), where
## com(t) is a time varying proxy for the baseline hazard from community informed 
##  by state surveillance
## para4 is a scaling parameter for the baseline community hazard
## para5-para8 are parameters for the effects of various subject characteristics
## 
## The daily hazard of infection for each household contact, j, from each infected 
## individual, i, in the household is modeled as:
## HH.HAZ.SUM.FN = SI.FN(t)*para3*exp(para5*AGECAT2+para6*AGECAT3+para7*HR+para8*vax), where
## SI.FN is a weibull distribution that models the serial interval 
## (time in days between symptom onset of prior and subsequent influenza cases) 
## in those households where influenza has been introduced:
##         exp(-(t/para1)^para2)-exp(-((t)/para1)^para2), where
##         para1 and para2 are shape parameters for the weibull distribution
## para3 is a constant parameter representing the baseline household hazard of infection
## para5-para8 are parameters for the effects of various subject characteristics
##
## The total hazard of infection is the sum of the hazards from the household and from
## the community.
##
## The contribution to the likelihood for those NOT infected during follow up, and
## prior to the date of illness onset for those who were infected is:
##     exp(-1*pre.haz)
##
##     pre.haz is the sum of the hazards of infection from the community (pre.c.haz ) and
##     from the household (pre.hh.haz) prior to infection
##
## The additional contribution to the likelihood on the date of illness onset 
## for those infected during follow up is:
##     1-(exp(-1*on.haz))
##
##     on.haz is the sum of the hazards of infection from the community (ondate.c.haz ) and
##     from the household (ondate.hh.haz) on the date of illness onset for those infected
#####################################################################################

# MESS required for the auc() function used in constraining the shape parameters of the SI.FN
rm(list = ls())
if(!suppressPackageStartupMessages(require(MESS)))
  install.packages("MESS")

# Read and Format the data
setwd("S:/MartinEpi/Projects/BoV-Uganda/Graphs and Reports/HH-transmission-model-master_2/premcmc_data")
ondate.c.exp <- readRDS("ondate.c.exp.rds")
ondate.hh.exp <- readRDS("ondate.hh.exp.rds")
preinf.c.exp <- readRDS("preinf.c.exp.rds")
preinf.hh.exp <- readRDS("preinf.hh.exp.rds")

#### mcmc function

#Josh: para1: weibull shape param 1, para2: weibull shape param 2, para3: hh base hazard, 
#      para4: com hazard scaling, para5: age 9-17, para6: age 18+, para7: High Risk, 
#      para8: vax

# prior function here
# usually we set non-informative (flat) priors, e.g. Uniform(a,b) or Normal(mu,sigma) with 
# relatively large sigma
logprior <- function(current)
{
  # write your own prior function here (CHANGE)
  current$logprior=dunif(current$para1,0,10,log=TRUE)+dunif(current$para2,0,10,log=TRUE)+
    dunif(current$para3, 0, 1, log=TRUE)+dunif(current$para4, 0,1, log=TRUE)+
    dnorm(current$para5, 0,100, log=TRUE)+dnorm(current$para6, 0,100, log=TRUE) +
    dnorm(current$para7, 0,100, log=TRUE) +dnorm(current$para8, 0,100, log=TRUE)  
  current
}

# the likelihood function
loglike <- function(current,PRE.HH, PRE.C, ON.HH, ON.C) 
{
  
    #SI.FN is a weibull function that models the serial interval distribution in households
    #where influenza has been introduced
    SI.FN <-function(current, EXPDAY) 
    { 
      exp(-(EXPDAY/current$para1)^current$para2)-
        exp(-((EXPDAY+1)/current$para1)^current$para2)
    }
    
    #HH.HAZ.SUM.FN models the daily hazard of infection for each household contact, j, 
    #from each infected individual, i, in the household
    HH.HAZ.SUM.FN <-function(INPUT, current) 
    {
      d <- as.data.frame(INPUT)$expday 
      sum (SI.FN(current, as.data.frame(INPUT)$expday) * current$para3 * 
             exp(current$para5*(as.data.frame(INPUT)$AGECAT2) + 
                 current$para6*(as.data.frame(INPUT)$AGECAT3) +
                 current$para7*(as.data.frame(INPUT)$HR) +
                 current$para8*(as.data.frame(INPUT)$vax)), na.rm=TRUE)
    } 
    
    #Apply HH.HAZ.SUM.FN to the data (formatted in C:/FILEPATH/THmodel_pre_mcmc_commented.R)
    #to obtain the hazard of household infection prior to illness onset...
    pre.hh.haz <-sapply(PRE.HH, HH.HAZ.SUM.FN, current=current)
    #...and on the day of illness onset for those infected
    on.hh.haz <-sapply(ON.HH, HH.HAZ.SUM.FN, current=current)
    
    #C.HAZ.SUM.FN models the daily hazard of infection from the community for each 
    #cohort member, j
    C.HAZ.SUM.FN <-function(INPUT,current) {
      com <-as.data.frame(INPUT)$com
      sum(com * current$para4 * (exp(current$para5 * (as.data.frame(INPUT)$AGECAT2) + 
                                     current$para6 * (as.data.frame(INPUT)$AGECAT3) +
                                     current$para7 * (as.data.frame(INPUT)$HR) +
                                     current$para8 * (as.data.frame(INPUT)$vax)
                                    )
                                )
         )
    }
    
    # Apply C.HAZ.SUM.FN to the data (formatted in C:/FILEPATH/THmodel_pre_mcmc_commented.R)
    # to obtain the hazard of community infection prior to illness onset...
    pre.c.haz <- sapply(PRE.C, C.HAZ.SUM.FN, current = current)
    # ...and on the day of illness onset for those infected
    on.c.haz <- sapply(ON.C, C.HAZ.SUM.FN, current = current)
    
    # Sum the household and community hazards to get total daily hazards
    #*******************************************************************
    # NOTE: The hazard vectors need to be the same length to be summed, so 
    # PRE.C and PRE.HH should have 1 element per j
    # ON.C and ON.HH should have 1 element per i
    #*******************************************************************

    pre.haz <- pre.hh.haz+pre.c.haz 
    on.haz <-on.hh.haz+on.c.haz
    
    #The contribution to the likelihood for those NOT infected during follow up, and
    #prior to the date of illness onset for those who were infected is:
    like.not<-exp(-1*pre.haz)

    #The additional contribution to the likelihood on the date of illness onset 
    #for those infected during follow up is:
    like.inf <-1-(exp(-1*on.haz))
    
    #Total log-likelihood
    current$loglike <-sum(log(like.inf+0.000001), log(like.not+0.0000001), na.rm=TRUE)
    current

}

#To assure model convergence, we use goodqk to constrain the Weibull shape parameters 
#such that the corresponding serial interval function would result in >=80% of household 
#secondary infections occurring within 14 days of the onset of illness in the index case
goodqk <-function(current) 
{
  SI.FN <-function(current, EXPDAY) 
  { 
    exp(-(EXPDAY/current$para1)^current$para2)-exp(-((EXPDAY+1)/current$para1)^current$para2)
  }
  temp <- data.frame(t=0:30, pmass=NA)
  temp$pmass <-SI.FN(current, temp$t)
  select <-temp[temp$pmass>=0,]
  allA <-auc(select$t, select$pmass)
  b14 <- auc(select$t[select$t<14], select$pmass[select$t<14])
  (b14/allA)>=0.8 | temp$pmass[temp$t==0] >=0
}
  
# set parameter constrains (when applicable)
bad_config<-function(current)
{
  ok_config=TRUE
  if(current$para1<0 | current$para2<0 | current$para3<0 | current$para4<=0 | current$para3>=1 | goodqk(current)!=TRUE) ok_config=FALSE   # write your constrains here, e.g., para1 should always be non-negative   (CHANGE)
  if(is.na((current$loglike))) ok_config=FALSE 
  
  bad=FALSE;if(!ok_config) bad=TRUE
  bad
}

# Metropolis-Hastings method for updating the parameter values
metropolis <-function(OLD,current,PRE.HH, PRE.C, ON.HH, ON.C)
{
  reject=bad_config(current)
  if(!reject)
  {
    current=logprior(current)
    current=loglike(current,PRE.HH, PRE.C, ON.HH, ON.C)
    lu=log(runif(1))
    if(is.na(current$logprior+sum(current$loglike)-
             OLD$logprior-sum(OLD$loglike))) reject=TRUE
    if(!is.na(current$logprior+sum(current$loglike)-
              OLD$logprior-sum(OLD$loglike))){
      if(lu>(current$logprior+sum(current$loglike)-
             OLD$logprior-sum(OLD$loglike))) reject=TRUE
    }  
  }
  if(reject)current=OLD
  current
}

#main function for MCMC process
mcmc <-function(PRE.HH, PRE.C, ON.HH, ON.C,MCMC_iterations,BURNIN_iterations,THINNING=1)
{
  #give initial values for parameters    
  current=list(para1=2, para2=2, para3=0.2, para4=0.005, para5=0.05, para6=0.05, 
               para7=0.05,  para8=-0.5,logprior=0,loglike=0)
  #write the variable to save the posteriors of the parameters   
  dump=list(para1=c(),para2=c(),para3=c(),para4=c(),para5=c(),para6=c(),
            para7=c(),para8=c(),  
            #Calculate the acceptance rate (45%-55% is good)
            accept=matrix(NA,ncol=8,nrow=MCMC_iterations))
  #set the sigma for updating parameters, adjust until satisfied with acceptance rate
  sigma = list(para1=1,para2=5,para3=.2, para4=0.001, para5=.75, para6=0.5, 
               para7=.75, para8=.5)
  #write logprior for initial values
  current=logprior(current)
  #write loglike for initial values
  current=loglike(current,PRE.HH, PRE.C, ON.HH, ON.C)
  #write posteriors at each iteration
  pb <- txtProgressBar(min = 0, max = MCMC_iterations, style = 3)
  for (iteration in (-BURNIN_iterations+1):MCMC_iterations)
  {
    old=current; current$para1=rnorm(1,current$para1,sigma$para1); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)    
    old=current; current$para2=rnorm(1,current$para2,sigma$para2); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)    
    old=current; current$para3=rnorm(1,current$para3,sigma$para3); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)    
    old=current; current$para4=rnorm(1,current$para4,sigma$para4); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)    
    old=current; current$para5=rnorm(1,current$para5,sigma$para5); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)    
    old=current; current$para6=rnorm(1,current$para6,sigma$para6); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)    
    old=current; current$para7=rnorm(1,current$para7,sigma$para7); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)   
    old=current; current$para8=rnorm(1,current$para8,sigma$para8); current=metropolis(old,current,PRE.HH, PRE.C, ON.HH, ON.C)    
    
    # save the posteriors after discarding burn-in values
    if(iteration>0) 
    {
      rr <- iteration
      dump$para1[rr]=current$para1   
      dump$para2[rr]=current$para2   
      dump$para3[rr]=current$para3   
      dump$para4[rr]=current$para4   
      dump$para5[rr]=current$para5   
      dump$para6[rr]=current$para6   
      dump$para7[rr]=current$para7   
      dump$para8[rr]=current$para8   

      # save the result for accept/reject at each step to help find appropriate sigma
      if(rr>1) {    
        dump$accept[rr, 1] <- 1 * (dump$para1[rr] != dump$para1[rr - 1])    
        dump$accept[rr, 2] <- 1 * (dump$para2[rr] != dump$para2[rr - 1])    
        dump$accept[rr, 3] <- 1 * (dump$para3[rr] != dump$para3[rr - 1])    
        dump$accept[rr, 4] <- 1 * (dump$para4[rr] != dump$para4[rr - 1])    
        dump$accept[rr, 5] <- 1 * (dump$para5[rr] != dump$para5[rr - 1])    
        dump$accept[rr, 6] <- 1 * (dump$para6[rr] != dump$para6[rr - 1])    
        dump$accept[rr, 7] <- 1 * (dump$para7[rr] != dump$para7[rr - 1])    
        dump$accept[rr, 8] <- 1 * (dump$para8[rr] != dump$para8[rr - 1])  
      }
    }
  setTxtProgressBar(pb, iteration)
  }
  dump
}

#start with small number of iterations (about 1000) to tune sigmas, then do full number
iter = 15000; burnin = iter / 3
set.seed(12345)

#run the mcmc
ptm <- proc.time()
a <- mcmc(preinf.hh.exp, preinf.c.exp, ondate.hh.exp, ondate.c.exp, iter, burnin, THINNING = 1)
proc.time() - ptm
#save the generated data
outdir<-"C:/FILEPATH/"
write.csv(a, paste(outdir, Sys.Date(),"mcmc_data.csv", sep="" ))

#generate output
source("C:/Users/Josh/Box Sync/Working Docs/Aim 3 Model/Independent Study/Output.R")