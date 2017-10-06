#########################################################
# title: "TH Model Infection Simulator"
# author: "Sophia Ng, Josh Petrie, Richard R. Evans"
# date: 9/27/2017
#########################################################

suppressPackageStartupMessages(library(tidyverse))

# Description

# Created by Josh Petrie on Feb 7, 2016 
# This script simulates daily community and household-acquired influenza 
# infections using a dataset with household structure and covariate data only,
# proxy data on weekly community influenza incidence and model parameters 
# estimated in THmodel_mcmc_commented.R 
# 
# There is probably a way easier way to do this given how repetitive the code is
# 
# # Read in data sets
# 
# Read in clean dataset, one row containing required info of each individual in 
# the cohort including household ID, individual ID, infection indicator, 
# end of follow up date, and all covariates to be included in regression, 
# but EXCLUDING illness onset dates


setwd("S:/MartinEpi/Projects/BoV-Uganda/Graphs and Reports/HH-transmission-model-master_2/")
# setwd("~/Documents/Martin Lab/HBov-1/Uganda/HH-transmission-model-master/")
formdat1 <- read.csv("THmodel_exampleData.csv")

glimpse(formdat1)

# Read in proxy community hazard data, one row for each date of study follow up 
# including counts of reported influenza cases standardized to the peak 
# (i.e proxy in peak week ==1)

setwd("S:/MartinEpi/Projects/BoV-Uganda/Graphs and Reports/HH-transmission-model-master_2/")
# setwd("~/Documents/Martin Lab/HBov-1/Uganda/HH-transmission-model-master/")
com <- read.csv("THmodel_exampleComData.csv")

# something weird at the end of the file
com <- com[1:105, ]
glimpse(com)

# theta: vector of parameter estimates from THmodel_mcmc_commented
# 
# they are as follows: q (weibull shape), k (weibull scale), HH baseline hazard,
# community hazard scale, age 9-17 (coeff), age 18+ (coeff), 
# high risk condition (coeff), vaccination (coeff)


# theta <- c(q, k, hh base haz, com haz scale, 
#            Age 9-17, Age 18+, HighRisk, Vaccination)
# q and k = weibull shape params

theta <- c(2.869142, # q - weibull shape
           3.296173, # k - weibull scale
           0.239762, # HH baseline hazard
           0.001353, # community hazard scale
           -1.26998, # age 9-17 estimate
           -1.03735, # age 18+ estimate
           -0.05432, # high risk condition estimate
           -0.31826) # vaccination estimate

# number of simulation iterations

iter <- 1000

# Simulation

# Simulation function to estimate total numbers of community and 
# household-acquired infections

Sim <- function(formdat1, com, theta, iter){
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  # empty matrix to store the daily number of infections in each iter
  Infections <- matrix(data = NA, nrow = iter, ncol = length(com[,1]))
  
  # set empty vectors/data frame
  vectors <- c("PRI",     # total number of 1st infection
               "PRIexp",  # total number at risk of 1st infection
               "PRIrisk", # calculated risk of 1st infection
               "SEC",     # total number of 2nd infections
               "SECexp",  # total number at risk of 2nd infx (exposed to PRI)
               "SECrisk", # calculated risk of 2nd infection
               "TER",     # total number of 3rd infections
               "TERexp",  # total number at risk of 3rd infx (exposed to SEC)
               "TERrisk", # calculated risk of 3rd infection
               "QUA",     # total number of quartenary infections
               "QUAexp",  # total number at risk of 4th infx (exposed to TER)
               "QUArisk", # calculated risk of 4th infection
               "HHrisk")  # calculated total risk of HH-acquired infection
  
  InfTypes <- data.frame(matrix(nrow = iter, ncol = length(vectors)))
  colnames(InfTypes) <- vectors
  # Can similarly set up empty datasets to capture number of infections stratified by 
  # covariates such as vaccination here:
  # Vax1 <- rep(NA,iter)
  # Vax0 <- rep(NA,iter)
  # InfVax <- data.frame(Vax1, Vax0)
  
  #loop over each iteration
  for (z in 1:iter) {
    
    #################################################
    # write in d(j) of community exposures starting #
    # from start of season to end                   #
    #################################################
    
    # create an empty list to write in daily community 
    # infection risk for each subject
    c.risk <- vector("list", nrow(formdat1))
    # loop through each individual j
    for (j in 1:nrow(formdat1)) { 
      # create empty vectors to keep track of current j's ID and hhID
      curhhID <- curID <- temp <- NA 
      # set the day that follow up ends
      end <- formdat1$endfudate[j] 
      # write in hhID and ID of the current j
      curhhID <- formdat1$hhID[j]
      curID <- formdat1$ID[j] 
      
      c.risk[[j]] <- data.frame(
        date = seq(0, end, by = 1),
        ID = formdat1$ID[j],
        hhID = formdat1$hhID[j], 
        AGECAT2 = formdat1$AGECAT2[j], 
        AGECAT3 = formdat1$AGECAT3[j], 
        HR = formdat1$HR[j],
        vaxdate = formdat1$vaxdate[j],
        vax = NA,
        cRisk = NA
      )
      # loop through each follow-up day to calculate daily 
      # community infection risk
      for (i in 1:length(c.risk[[j]][,1])) {
        #calculate time-varying vaccination status
        c.risk[[j]]$vax[i] <- ifelse(c.risk[[j]]$vaxdate[i] > c.risk[[j]]$date[i],
                                     0, 1)
        #calculate daily risk of infection from the community
        inner <- exp(theta[5] * c.risk[[j]]$AGECAT2[i] + 
                       theta[6] * c.risk[[j]]$AGECAT3[i] +
                       theta[7] * c.risk[[j]]$HR[i] +
                       theta[8] * c.risk[[j]]$vax[i])
        c.risk[[j]]$cRisk[i] <- 1 - exp(-(com$COM[i] * theta[4] * inner))
      } 
    } # one day per row, from date of index onset to end
    
    #####################################################
    # Simulate onset dates for infection from community #
    #####################################################
    
    #create empty variable to store onset dates of primary cases
    formdat1$ondate <- NA
    for ( j in 1:nrow(formdat1)) {
      c.risk[[j]]$ondate <- NA
      for (i in 1:length(c.risk[[j]][, 1])) {
        #sample value from uniform distribution
        cProb <- runif(1, 0, 1)
        # if the sampled value is less than the community risk 
        # of infection on day 1 of follow up then record the onset date
        if(i == 1) {
          if(cProb < c.risk[[j]]$cRisk[i]) { 
            c.risk[[j]]$ondate[i] <- i - 1
            formdat1$ondate[j] <- c.risk[[j]]$ondate[i]
          }
        } else {
          # simulate infections for remaining days of follow up
          # if not already infected
          if(is.na(c.risk[[j]]$ondate[i - 1])) {
            if(cProb < c.risk[[j]]$cRisk[i]) {
              c.risk[[j]]$ondate[i] <- i - 1
              formdat1$ondate[j] <- c.risk[[j]]$ondate[i]
            } # end inner if
          } # end outer if
        } # end if-else
      } # end for over i in c.risk[[j]]
    } # end for over j
    
    #store data of those that were infected
    pos <- formdat1[which(!is.na(formdat1$ondate)), ]
    if(nrow(pos) > 0) {
      #identify first individual infected in the household
      pos1 <- aggregate(ondate ~ hhID, min, data = pos)
      #record date of first household infection in dataset
      formdat1$firstpos <- NA
      for(i in 1:length(formdat1[,1])) {
        for(j in 1:length(pos1[,1])) {
          if(pos1$hhID[j] == formdat1$hhID[i]) {
            formdat1$firstpos[i] <- pos1$ondate[j]
          }
        }
      }
      
      rm(i); rm(j); rm(pos); rm(pos1)
      #list of households with at least 1 infection so far
      ihhIDs <-unique(formdat1$hhID[which(!is.na(formdat1$firstpos))])
      #create list of all infected individuals and their data
      infected <- vector("list", length(ihhIDs))
      
      for (i in 1:length(ihhIDs)) {
        pos <- formdat1[!is.na(formdat1$ondate),]  
        infected[[i]] <- data.frame(
          hhID = ihhIDs[i],
          ID = pos$ID[pos$hhID %in% ihhIDs[i]],
          ondate = pos$ondate[pos$hhID %in% ihhIDs[i]], 
          endfudate = pos$endfudate[pos$hhID %in% ihhIDs[i]]
        ) 
      }
      
      ###################################################
      # generate daily risk of HH infx of j from each i #
      # who were infected before j                      #
      ###################################################
      
      hh.risk <- vector("list", nrow(formdat1))
      # loop through each j
      for ( j in 1:nrow(formdat1)) { 
        # create empty vectors to keep track of current j's ID and hhID
        curhhID <- curID <- temp <- NA 
        # set the day that follow up ends
        end <- formdat1$endfudate[j] 
        # write in hhID and ID of the current j
        curhhID <- formdat1$hhID[j]
        curID <- formdat1$ID[j] 
        
        # if onset of i was after follow up ended, or there was no i, or 
        # j was already infected then hh exposure <- NA
        if (formdat1$firstpos[j] > end  | 
            is.na(formdat1$firstpos[j]) |
            !is.na(formdat1$ondate[j])) {
          hh.risk[[j]] <- NA
          next 
        } else {
          # get data of all i in the same household
          temp <- infected[[which(ihhIDs %in% formdat1$hhID[j])]] 
          # exclude i if his onset date was later than the end of follow-up
          temp <- temp[temp$ondate <= end, ] 
          # do not include j himself as i
          temp <- temp[!(temp$ID %in% curID), ] 
          # if no eligible i, hh exposure <- NA
          if (nrow(temp) == 0) {
            hh.risk[[j]] <- NA
            next 
          } else {
            # create an empty list, one element contains one eligible i for j
            exp <- vector("list", nrow(temp)) 
            # Loop through each i
            for (i in 1:nrow(temp)) { 
              # get season dates, from i's onset to i's onset +14 days
              mindays <- min(end, (temp$ondate[i] + 14))
              tempdays <- seq(temp$ondate[i], mindays, by = 1) 
              # one row per each season date
              index <- rep(row.names(temp[i, ]), length(tempdays))
              exp[[i]] <- temp[i,][index, ]
              # write season dates
              exp[[i]]$date <- tempdays 
              # write exposure dates, d(i,j), 
              # as season date minus i's onset date 
              exp[[i]]$expday <- as.numeric(exp[[i]]$date - exp[[i]]$ondate) 
              # write covariate (here, 9 >= AGE < 18 of j)
              exp[[i]]$AGECAT2 <- formdat1$AGECAT2[j] 
              # write covariate (here, AGE =< 18 of j)
              exp[[i]]$AGECAT3 <- formdat1$AGECAT3[j] 
              # write covariate (here, HIGH RISK STATUS of j)
              exp[[i]]$HR <- formdat1$HR[j] 
              # write time varying covariate (here, VACCINATION STATUS of j)
              exp[[i]]$vaxdate <- formdat1$vaxdate[j]
              for ( k in 1:nrow(exp[[i]]))
                exp[[i]]$vax[k] <- ifelse(exp[[i]]$date[k] < exp[[i]]$vaxdate[k], 0, 1)
              
              ###################################################
              # calculate daily js risk of infection from       #
              # household exposure to i                         #
              # Fixed an error in this part of code calculating #
              # hhRisk 12-12-2016                               #
              # Generated from Weibull distribution             #
              ###################################################
              
              linpred <- exp(theta[5] * exp[[i]]$AGECAT2 +
                               theta[6] * exp[[i]]$AGECAT3 +
                               theta[7] * exp[[i]]$HR +
                               theta[8] * exp[[i]]$vax) * theta[3]
              dayof <- (exp[[i]]$expday / theta[1]) ^ theta[2]
              dayb4 <- ((exp[[i]]$expday + 1) / theta[1]) ^ theta[2]
              exp[[i]]$hhRisk <- 1 - exp(-(exp(-dayof) - exp(-dayb4)) * linpred)
            }  
            # collaspe 1st to ith exp into the jth element in hh.risk
            hh.risk[[j]] <- do.call(rbind, exp) # make data frame from list
          }
        }
      }
      
      ###################################################
      #  Simulate onset dates for secondary infections  # 
      ###################################################
      
      #create empty variable to store onset dates of secondary cases
      formdat1$ondate2 <- NA
      
      for ( j in 1:nrow(formdat1)){
        if (is.na(hh.risk[[j]])) next  # warnings introduced here: condition has length > 1
        hh.risk[[j]]$ondate <- NA
        for (i in 1:length(hh.risk[[j]][,1])) {
          #sample value from uniform distribution
          hhProb <- runif(1,0,1)
          #if the sampled value is less than the household risk of infection on day 1 of
          #follow up then record the onset date
          if(i == 1) {
            if(hhProb < hh.risk[[j]]$hhRisk[i]) {
              hh.risk[[j]]$ondate[i] <- hh.risk[[j]]$date[i]
              formdat1$ondate2[j] <- hh.risk[[j]]$ondate[i]
            }
          } else {
            #Simulate infections for remaining days of follow up if not already infected
            if(is.na(hh.risk[[j]]$ondate[i - 1])){
              if(hhProb < hh.risk[[j]]$hhRisk[i]) {
                hh.risk[[j]]$ondate[i] <- hh.risk[[j]]$date[i]
                formdat1$ondate2[j] <- hh.risk[[j]]$ondate[i]
              }
            }
          } # end if-else
        } # end for over i in hh.risk[[j]]
      } # end for over j individuals
      
      # store data of those that were infected
      HHpos <- formdat1[which(!is.na(formdat1$ondate2)), ] 
      if(nrow(HHpos) > 0) {
        # identify first secondary infection in the household
        HHpos1 <- aggregate(ondate2 ~ hhID, min, data = HHpos)
        # record date of first secondary infection in dataset
        formdat1$secondpos <- NA
        for(i in 1:length(formdat1[, 1])) {
          for(j in 1:length(HHpos1[, 1])) {
            if(HHpos1$hhID[j] == formdat1$hhID[i]) {
              formdat1$secondpos[i] <- HHpos1$ondate2[j]
            }
          }
        }
        
        rm(i); rm(j); rm(HHpos); rm(HHpos1)
        
        #########################################################
        #                 Secondary Infections                  #
        #########################################################
        
        # list of households with at least 2 infections so far
        ihhIDs2 <- unique(formdat1$hhID[which(!is.na(formdat1$secondpos))])
        infected2 <- vector("list", length(ihhIDs2))
        # create list of all infected individuals in households  
        # with at least 2 infections and their data             
        for (i in 1:length(ihhIDs2)) {
          pos <- formdat1[!is.na(formdat1$ondate2), ]  
          infected2[[i]] <- data.frame(
            hhID = ihhIDs2[i],
            ID = pos$ID[pos$hhID %in% ihhIDs2[i]],
            ondate2 = pos$ondate2[pos$hhID %in% ihhIDs2[i]],
            endfudate = pos$endfudate[pos$hhID %in% ihhIDs2[i]]
          ) 
        }
        # empty list to write in the daily risk of household infection of
        # j from each i who were infected before j; i are secondary infections
        hh.risk2 <- vector("list", nrow(formdat1))
        # loop through each j
        for (j in 1:nrow(formdat1)) { 
          # create empty vectors to keep track of current j's ID and hhID
          curhhID <- curID <- temp <- NA 
          # set the day that follow up ends
          end <- formdat1$endfudate[j] 
          # write in hhID and ID of the current j
          curhhID <- formdat1$hhID[j]
          curID <- formdat1$ID[j] 
          
          # if onset of i was after follow up ended, or there was no i, or 
          # j was already infected then hh exposure <- NA
          if (formdat1$secondpos[j] > end  |
              is.na(formdat1$secondpos[j]) |
              !is.na(formdat1$ondate[j])   |
              !is.na(formdat1$ondate2[j])) {
            hh.risk2[[j]] <- NA
            next 
          } else {
            # get data of all i in the same household
            temp <- infected2[[which(ihhIDs2%in%formdat1$hhID[j])]] 
            # exclude i if his onset date was later than the end of follow-up
            temp <- temp[temp$ondate2 <= end, ] 
            # do not include j himself as i
            temp <- temp[!(temp$ID %in% curID), ] 
            # if no eligible i, hh exposure <- NA
            if (nrow(temp) == 0) {
              hh.risk2[[j]] <- NA
              next 
            } else {
              # create an empty list, one element contains one eligible i for j
              exp <- vector("list", nrow(temp)) 
              # Loop through each i
              for ( i in 1:nrow(temp)) { 
                # get season dates, from i's onset to i's onset +14 days
                mindays <- min(end, (temp$ondate2[i] + 14))
                tempdays <- seq(temp$ondate2[i], mindays, by = 1) 
                # one row per each season date
                index <- rep(row.names(temp[i, ]), length(tempdays))
                exp[[i]] <- temp[i,][index, ] 
                # write season dates
                exp[[i]]$date <- tempdays 
                # write exposure dates, d(i,j), as season date minus i's onset date
                exp[[i]]$expday <- as.numeric(exp[[i]]$date - exp[[i]]$ondate2) 
                # write covariate (here, 9 >= AGE < 18 of j)
                exp[[i]]$AGECAT2 <- formdat1$AGECAT2[j] 
                # write covariate (here, AGE =< 18 of j)
                exp[[i]]$AGECAT3 <- formdat1$AGECAT3[j] 
                # write covariate (here, HIGH RISK STATUS of j)
                exp[[i]]$HR <- formdat1$HR[j]  
                # write time varying covariate (here, VACCINATION STATUS of j)
                exp[[i]]$vaxdate <- formdat1$vaxdate[j]
                for ( k in 1:nrow(exp[[i]]))
                  exp[[i]]$vax[k] <- ifelse(exp[[i]]$date[k] < exp[[i]]$vaxdate[k], 0, 1)
                
                ###################################################
                # calculate daily js risk of infection from       #
                # household exposure to i                         #
                # Fixed an error in this part of code calculating #
                # hhRisk 12-12-2016                               #
                # Generated from Weibull distribution             #
                ###################################################
                
                linpred <- exp(theta[5] * exp[[i]]$AGECAT2 +
                                 theta[6] * exp[[i]]$AGECAT3 + 
                                 theta[7] * exp[[i]]$HR + 
                                 theta[8] * exp[[i]]$vax) * theta[3]
                dayof <- (exp[[i]]$expday / theta[1]) ^ theta[2]
                dayb4 <- ((exp[[i]]$expday + 1) / theta[1]) ^ theta[2]
                exp[[i]]$hhRisk2 <- 1 - exp( -(exp(-dayof) - exp(-dayb4)) * linpred)
              }  
              # collaspe 1st to ith exp into the jth element in hh.risk2
              # make data frame of secondary hh infx
              hh.risk2[[j]] <- do.call(rbind, exp)
            }
          } 
        }
        
        #########################################################
        #                 Tertiary Infections                  #
        #########################################################
        
        #create empty variable to store onset dates of tertiary cases
        formdat1$ondate3 <- NA
        #Simulate onset dates for tertiary infections
        for (j in 1:nrow(formdat1)) {
          if (is.na(hh.risk2[[j]])) next
          hh.risk2[[j]]$ondate <- NA
          for (i in 1:length(hh.risk2[[j]][,1])) {
            #sample value from uniform distribution
            hhProb <- runif(1, 0, 1)
            #if the sampled value is less than the household risk of infection on day 1 of
            #follow up then record the onset date
            if(i == 1) {
              if(hhProb<hh.risk2[[j]]$hhRisk2[i]){
                hh.risk2[[j]]$ondate[i] <- hh.risk2[[j]]$date[i]
                formdat1$ondate3[j] <- hh.risk2[[j]]$ondate[i]
              }
            } else {
              #Simulate infections for remaining days of follow up if not already infected
              if(is.na(hh.risk2[[j]]$ondate[i - 1])){
                if(hhProb < hh.risk2[[j]]$hhRisk2[i]){
                  hh.risk2[[j]]$ondate[i] <- hh.risk2[[j]]$date[i]
                  formdat1$ondate3[j] <- hh.risk2[[j]]$ondate[i]
                }
              }
            } # end if-else
          } # end loop over i in j
        } # end loop over j individuals
        
        #store data of those that were infected
        HH2pos <- formdat1[which(!is.na(formdat1$ondate3)), ]
        if(nrow(HH2pos) > 0) {
          #identify first tertiary infection in the household
          HH2pos1 <- aggregate(ondate3 ~ hhID, min, data = HH2pos)
          
          #record date of first tertiary infection in dataset
          formdat1$thirdpos = NA
          for(i in 1:length(formdat1[, 1])) {
            for(j in 1:length(HH2pos1[, 1])) {
              if(HH2pos1$hhID[j] == formdat1$hhID[i]) {
                formdat1$thirdpos[i] <- HH2pos1$ondate3[j]
              }
            }
          }
          
          rm(i); rm(j); rm(HH2pos); rm(HH2pos1)
          #list of households with at least 3 infections so far
          ihhIDs3 <-unique(formdat1$hhID[which(!is.na(formdat1$thirdpos))])
          #create list of all infected individuals in households with at least 3 infections
          #so far and their data
          infected3 <- vector("list", length(ihhIDs3))
          
          for ( i in 1:length(ihhIDs3)) {
            pos <- formdat1[!is.na(formdat1$ondate3), ]  
            infected3[[i]] <- data.frame(
              hhID = ihhIDs3[i], 
              ID = pos$ID[pos$hhID %in% ihhIDs3[i]],
              ondate3 = pos$ondate3[pos$hhID %in% ihhIDs3[i]],
              endfudate = pos$endfudate[pos$hhID %in% ihhIDs3[i]]
            )
          }
          # create an empty list to write in the daily risk of household infection of
          #j from each i who were infected before j; i are tertiary infections
          hh.risk3 <- vector("list", nrow(formdat1))
          
          # loop through each j
          for (j in 1:nrow(formdat1)) { 
            # create empty vectors to keep track of current j's ID and hhID
            curhhID <- curID <- temp <- NA 
            # set the day that follow up ends
            end <- formdat1$endfudate[j] 
            # write in hhID and ID of the current j
            curhhID <- formdat1$hhID[j]
            curID <- formdat1$ID[j] 
            
            # if onset of i was after follow up ended, or there was no i, or 
            # j was already infected then hh exposure <- NA
            if (formdat1$thirdpos[j] > end  |
                is.na(formdat1$thirdpos[j]) |
                !is.na(formdat1$ondate[j])  |
                !is.na(formdat1$ondate2[j]) |
                !is.na(formdat1$ondate3[j])) {
              hh.risk3[[j]]<-NA
              next 
            } else {
              # get data of all i in the same household
              index <- which(ihhIDs3 %in% formdat1$hhID[j])
              temp <- infected3[[index]] 
              # exclude i if his onset date was later than either j's end
              temp <- temp[temp$ondate3 <= end, ] 
              # do not include j himself as i
              temp <- temp[!(temp$ID %in% curID), ] 
              # if no eligible i, hh exposure <- NA
              if (nrow(temp) == 0) {
                hh.risk3[[j]] <- NA
                next 
              } else {
                # create an empty list, one element contains one eligible i for j
                exp <- vector("list", nrow(temp)) 
                # Loop through each i
                for ( i in 1:nrow(temp)) { 
                  # get season dates, from i's onset to i's onset +14 days
                  mindays <- min(end, (temp$ondate2[i] + 14))
                  tempdays <- seq(temp$ondate3[i], mindays, by = 1) 
                  # one row per each season date
                  index <- rep(row.names(temp[i,]), length(tempdays))
                  exp[[i]] <- temp[i,][index, ] 
                  # write season dates
                  exp[[i]]$date <- tempdays 
                  # write exposure dates, d(i,j), as season date minus i's onset date
                  exp[[i]]$expday <- as.numeric(exp[[i]]$date - exp[[i]]$ondate3) 
                  # write covariate (here, 9 >= AGE < 18 of j)
                  exp[[i]]$AGECAT2 <-formdat1$AGECAT2[j] 
                  # write covariate (here, AGE =< 18 of j)
                  exp[[i]]$AGECAT3 <- formdat1$AGECAT3[j] 
                  # write covariate (here, HIGH RISK STATUS of j)
                  exp[[i]]$HR <- formdat1$HR[j]  
                  # write time varying covariate (here, VACCINATION STATUS of j)
                  exp[[i]]$vaxdate <- formdat1$vaxdate[j]
                  for ( k in 1:nrow(exp[[i]]))
                    exp[[i]]$vax[k] <- ifelse(exp[[i]]$date[k] < exp[[i]]$vaxdate[k], 0, 1)
                  
                  ###################################################
                  # calculate daily js risk of infection from       #
                  # household exposure to i                         #
                  # Fixed an error in this part of code calculating #
                  # hhRisk 12-12-2016                               #
                  # Generated from Weibull distribution             #
                  ###################################################
                  
                  linpred <- exp(theta[5] * exp[[i]]$AGECAT2 +
                                   theta[6] * exp[[i]]$AGECAT3 + 
                                   theta[7] * exp[[i]]$HR + 
                                   theta[8] * exp[[i]]$vax) * theta[3]
                  dayof <- (exp[[i]]$expday / theta[1]) ^ theta[2]
                  dayb4 <- ((exp[[i]]$expday + 1) / theta[1]) ^ theta[2]
                  exp[[i]]$hhRisk3 <- 1 - exp( -(exp(-dayof) - exp(-dayb4)) * linpred)
                } # end loop over 
                # collaspe 1st to ith exp into the jth element in hh.risk3
                hh.risk3[[j]] <- do.call(rbind, exp)
              } # end if-else-else 
            } # end if-else
          } # end loop over jth individual in formdat1 data
          
          #################################################
          # Simulate onset dates for quatenary infections #
          #################################################
          
          #create empty variable to store onset dates of quartenary cases
          formdat1$ondate4 <- NA
          #Simulate onset dates for quatenary infections
          for (j in 1:nrow(formdat1)) {
            if (is.na(hh.risk3[[j]])) next
            hh.risk3[[j]]$ondate <- NA
            for (i in 1:length(hh.risk3[[j]][,1])) {
              #sample value from uniform distribution
              hhProb <- runif(1, 0, 1)
              #if the sampled value is less than the household risk of infection on day 1 of
              #follow up then record the onset date
              if(i == 1) {
                if(hhProb < hh.risk3[[j]]$hhRisk3[i]) {
                  hh.risk3[[j]]$ondate[i] <- hh.risk3[[j]]$date[i]
                  formdat1$ondate4[j] <- hh.risk3[[j]]$ondate[i]
                }
              } else {
                #Simulate infections for remaining days of follow up if not already infected
                if(is.na(hh.risk3[[j]]$ondate[i - 1])) {
                  if(hhProb < hh.risk3[[j]]$hhRisk3[i]){
                    hh.risk3[[j]]$ondate[i] <- hh.risk3[[j]]$date[i]
                    formdat1$ondate4[j] <- hh.risk2[[j]]$ondate[i]
                  }
                }
              } # end if-else
            } # end loop over i in hh.risk3   
          } # end loop over jth individual 
        } # end loop over secondary infection individuals
      } # end loop over primary infection individuals
    } # end loop over positive individuals
    
    # total number of primary infections
    InfTypes$PRI[z] <- sum(!is.na(formdat1$ondate))
    # total number at risk of primary infection 
    InfTypes$PRIexp[z] <- #Total number of individuals in cohort
      # calculate risk of primary infection
      InfTypes$PRIrisk[z] <- InfTypes$PRI[z] / InfTypes$PRIexp[z]
    # total number of secondary infections
    InfTypes$SEC[z] <- sum(!is.na(formdat1$ondate2))
    # total number at risk of secondary infection (exposed to primary case)
    InfTypes$SECexp[z] <- sum(!is.na(formdat1$firstpos) & is.na(formdat1$ondate))
    # calculate risk of secondary infection
    InfTypes$SECrisk[z] <- InfTypes$SEC[z] / InfTypes$SECexp[z]
    # total number of tertiary infections
    InfTypes$TER[z] <- sum(!is.na(formdat1$ondate3))
    # total number at risk of tertiary infection (exposed to secondary case)
    InfTypes$TERexp[z] <- sum(!is.na(formdat1$secondpos) &
                                is.na(formdat1$ondate)     &
                                is.na(formdat1$ondate2))
    # calculate risk of tertiary infection
    InfTypes$TERrisk[z] <- InfTypes$TER[z] / InfTypes$TERexp[z]
    # total number of quartenary infections
    InfTypes$QUA[z] <- sum(!is.na(formdat1$ondate4))
    # total number at risk of quartenary infection (exposed to tertiary case)
    InfTypes$QUAexp[z] <- sum(!is.na(formdat1$thirdpos) & 
                                is.na(formdat1$ondate)    &  
                                is.na(formdat1$ondate2)   &
                                is.na(formdat1$ondate3[j]))
    # calculate risk of quartenary infection
    InfTypes$QUArisk[z] <- InfTypes$QUA[z] / InfTypes$QUAexp[z]
    # calculate total risk of household infection
    all_infx <- InfTypes$SEC[z] + InfTypes$TER[z] + InfTypes$QUA[z]
    InfTypes$HHrisk[z] <- all_infx / InfTypes$SECexp[z]
    
    # overall onset dates to single variable
    formdat1$infdate <- formdat1$ondate
    index <- is.na(formdat1$infdate) & !is.na(formdat1$ondate2)
    formdat1$infdate[index] <- formdat1$ondate2[index]
    index <- is.na(formdat1$infdate) & !is.na(formdat1$ondate3)
    formdat1$infdate[index] <- formdat1$ondate3[index]
    index <- is.na(formdat1$infdate) & !is.na(formdat1$ondate4)
    formdat1$infdate[index] <- formdat1$ondate4[index]
    
    # get total number of infections on days that infections occurred
    infdate <- as.data.frame(table(formdat1$infdate))
    infdate
    
    # get total number of infections for all season days of follow-up
    for (d in 1:(formdat1$endfudate[1] + 1)) {
      Infections[z, d] <- 0
      for (r in 1:nrow(infdate)) {
        if(infdate$Var1[r] == (d - 1)) {
          Infections[z, d] <- infdate$Freq[r]
        }
      } 
    }
    setTxtProgressBar(pb, z)
  } # end loop over z in iter
  
  # output simulated infection data
  OUT <- list(Infections, InfTypes)
  OUT
}

# Run Simulation

system.time(PredInf <- Sim(formdat1, com, theta, iter))

