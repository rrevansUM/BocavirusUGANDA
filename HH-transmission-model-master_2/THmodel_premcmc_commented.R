#####################################################################################
## Created by Sophia Ng on Aug 31, 2015 ## Modified by Josh Petrie on Nov 2, 2015
## This script formats the data to be used in the likelihood function of the 
## individual-based transmission hazard model for analysis of influenza infection
## risks from the household and from the community in a household COHORT study. 
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
##
## The data are formatted into 4 lists with 1 element (dataframe) per individual, j 
## The first 2, preinf.hh.exp and preinf.c.exp, will be applied to exp(-1*pre.haz)
##     preinf.hh.exp -- contains data for each day of household exposure, d(i,j),
##        from each infected household contact, i, prior to j's illness onset date or 
##        follow up end.
##     preinf.c.exp -- contains data for each day of community exposure, d(j), prior to 
##        j's illness onset date or follow up end.
##        
## The final 2, ondate.hh.exp and ondate.c.exp, will be applied to 1-(exp(-1*on.haz))
##     ondate.hh.exp -- contains data for household exposure, d(i,j), from each infected 
##        household contact, i, on j's illness onset date.
##     ondate.c.exp -- contains data for community exposure, d(j), on 
##        j's illness onset date.
##
## Each dataframe within the 4 lists also contains co-variate info, eg. Age, etc.
#####################################################################################

# Read in clean dataset, one row containing required info of each individual in the cohort
# including household ID, individual ID, illness onset date, end of follow up date, and all
# covariates to be included in regression
formdat1 <-read.csv("C:/FILEPATH/FILE.csv")

# Read in proxy community hazard data, one row for each date of study follow up including
# counts of reported influenza cases standardized to the peak (i.e proxy in peak week ==1)
com<-read.csv("C:/FILPATH/FILE.csv")#home
#------------------------------------------------------------------------------------------

# create a list of all hhIDs for later indexing
hhIDs <-unique(formdat1$hhID)

# create a list of hhIDs where hh introduction of influenza occurred for later indexing
ihhIDs <-unique(formdat1$hhID[which(!is.na(formdat1$introduced))])

# Create a list of all the i (infected individuals) in each household
# each element contains 1 household, each row in each element contains 1 infected individual

infected <- vector("list", length(ihhIDs))

for ( i in 1:length(ihhIDs)) {
  pos <-formdat1[!is.na(formdat1$ondate),]  #  keep only i's w/ illness onset dates
# write each i's household ID, individual ID, own onset date, their end of followup date
  infected[[i]] <- data.frame(hhID=ihhIDs[i], ID=pos$ID[pos$hhID %in% ihhIDs[i]], 
                              ondate=pos$ondate[pos$hhID %in% ihhIDs[i]], 
                              endfudate=pos$endfudate[pos$hhID %in% ihhIDs[i]]) 
} 
#------------------------------------------------------------------------------------------

# Here we format the data for the contribution to the hazard function of time 
# at risk of household (d(i,j)) infection prior to infection or end of follow up 

# create an empty list to write in d(i,j) of hh exposure of j from each i infected prior to j
preinf.hh.exp <- vector("list", nrow(formdat1))

# loop through each j
for ( j in 1:nrow(formdat1)) {
  # create empty vectors to keep track of current j's ID and hhID
  curhhID <-curID<- temp <-NA 
  # empty vector for the last day of j's exposure  
  end <- as.Date(NA) 
  # write hhID and ID of the current j
  curhhID <-formdat1$hhID[j] ; curID <-formdat1$ID[j] 
  # end is either the studyend date (if j wasn't infected) or 
  # the day before ondate of j (if j was infected)
  if (is.na(formdat1$ondate[j])) {end <-formdat1$endfudate[j]
  } else {end <-(formdat1$ondate[j]-1)} 
  
  # write data for each hh exposure date, d(i,j), starting from onset date of i 
  # through +14 days after i's onset date or end, whichever comes first
  
  # if j has onset before first i or no hh introduction, hh exposure <- NA
  if (formdat1$firstpos[j]>end|is.na(formdat1$firstpos[j])) {preinf.hh.exp[[j]]<-NA; next 
  } else {
    # get info of all i in the same household
    temp <- infected[[which(ihhIDs%in%formdat1$hhID[j])]] 
    # exclude i if their onset date was later than j's end
    temp <- temp[temp$ondate<=end,] 
    # do not include current j as i
    temp <- temp[!(temp$ID %in% curID),] 
    
    # if no eligible i, hh exposure <-NA
    if (nrow(temp)%in%0) {preinf.hh.exp[[j]] <-NA ; next 
    } else {
        # create an empty list, one element contains one i for j
        exp <-vector("list", nrow(temp)) 
        # Loop through each i
        for ( i in 1:nrow(temp)) {
          # get season dates from i's onset to i's onset +14 days
          tempdays <-seq(temp$ondate[i], min(end,(temp$ondate[i]+14)), by=1)
          # one row per each season date
          exp[[i]] <- temp[i,][rep(row.names(temp[i,]),length(tempdays)),]
          # write season dates
          exp[[i]]$date <-tempdays 
          # write exposure dates, d(i,j), as season date minus i's onset date 
          exp[[i]]$expday <-as.numeric(exp[[i]]$date-exp[[i]]$ondate) 
          # write covariate (here, 9>=AGE<18 of j)
          exp[[i]]$AGECAT2 <-formdat1$AGECAT2[j] 
          # write covariate (here, AGE=<18 of j)
          exp[[i]]$AGECAT3 <-formdat1$AGECAT3[j] 
          # write covariate (here, HIGH RISK STATUS of j)
          exp[[i]]$HR <- formdat1$HR[j] 
          # write time varying covariate (here, VACCINATION STATUS of j)
          exp[[i]]$vaxdate <- formdat1$vaxdate[j]
          if (exp[[i]]$date<exp[[i]]$vaxdate){
            exp[[i]]$vax<-0
          } else {exp[[i]]$vax<-1}
          }  
        # collaspe 1st to ith exp into the jth element in preinf.hh.exp
        preinf.hh.exp[[j]] <- do.call(rbind,exp)} 
    
  }
}

#------------------------------------------------------------------------------------------

# Here we format the data for the contribution to the hazard function of time 
# at risk of community (d(j)) infection prior to infection or end of follow up

# create an empty list to write in d(j) of community exposure of each j 
preinf.c.exp <- vector("list", nrow(formdat1)) 

# loop through each j
for ( j in 1:nrow(formdat1)) { 
  # create empty vectors to keep track of current j's ID and hhID
  curhhID <-curID<- temp <-NA 
  # empty vector for the last day of j's exposure  
  end <- as.Date(NA) 
  # write hhID and ID of the current j
  curhhID <-formdat1$hhID[j] ; curID <-formdat1$ID[j] 
  # end is either the studyend date (if j wasn't infected) or 
  # the day before j's illness onset date (if j was infected)
  if (is.na(formdat1$ondate[j])) {end <-formdat1$endfudate[j]
  } else {end <-(formdat1$ondate[j]-1)} 
    
    # write data for each com exposure date, d(j), starting from the start of the season
    # through the end of the season or end of follow up, whichever comes first.
    # There is 1 element in the list per j, and 1 row in each element per exposure day, d(j)
    # Data include d(j), and all covariates to be included in regression model,
    # here age, high risk status, and time-varying vaccination status.
    preinf.c.exp[[j]] <-data.frame(date=seq(0, end, by=1 ), AGECAT2=formdat1$AGECAT2[j], 
                                   AGECAT3=formdat1$AGECAT3[j],HR=formdat1$HR[j],
                                   vaxdate=formdat1$vaxdate[j]) 
    for (i in 1:length(preinf.c.exp[[j]][,1])){
      preinf.c.exp[[j]]$com[i]<-com$COM[i]
      if (preinf.c.exp[[j]]$vaxdate[i]>preinf.c.exp[[j]]$date[i]) {
        preinf.c.exp[[j]]$vax[i]=0
      } else {preinf.c.exp[[j]]$vax[i]=1}
    }
    
}

#------------------------------------------------------------------------------------------

# Here we format the data for the contribution to the hazard function of time 
# at risk of household (d(i,j)) infection on the day of illness onset of infected js

# create empty list to write in d(i,j) of hh exposure of j from each i infected prior to j
ondate.hh.exp <- vector("list", nrow(pos))

# loop through each j
for ( j in 1:nrow(pos)) {  
  # create empty vectors to keep track of current j's ID and hhID
  curhhID <-curID<- temp  <-NA
  # empty vector for j's date of illness onset 
  curondate <- as.Date(NA)
  # write hhID and ID of the current j
  curhhID <-pos$hhID[j] ; curID <-pos$ID[j] 
  # write j's date of illness onset
  curondate <-pos$ondate[j] 
  
  # write data for each hh exposure, d(i,j) on j's date of illness onset
  
  # if j is the first individual infected in the household, hh exposure <- NA
  if (pos$index[j]==1 | pos$index2[j]==1) {ondate.hh.exp[[j]]<-NA; next
  } else { 
  # get info of all i in the same household  
  temp <- infected[[which(ihhIDs%in%pos$hhID[j])]]
  # exclude i if their onset date was later than j's illness onset date
  temp <- temp[temp$ondate<=curondate,] 
  # do not include current j as i
  temp <- temp[!(temp$ID%in%curID),]
  
  # create an empty list, one element contains one i for j
  exp <-vector("list", nrow(temp))
  
  # Loop through each i
  for ( i in 1:nrow(temp)) {
    # one row of data for j's illness onset date
    exp[[i]] <- temp[i,]
    # write season date of j's illness onset
    exp[[i]]$date <- curondate
    # write exposure date, d(i,j), as j's illness onset date minus i's illness onset date 
    exp[[i]]$expday <-as.numeric(exp[[i]]$date-exp[[i]]$ondate)
    # write covariate (here, 9>=AGE<18 of j)
    exp[[i]]$AGECAT2 <-pos$AGECAT2[j] 
    # write covariate (here, AGE=<18 of j)
    exp[[i]]$AGECAT3 <-pos$AGECAT3[j] 
    # write covariate (here, HIGH RISK STATUS of j)
    exp[[i]]$HR <- pos$HR[j]
    # write time varying covariate (here, VACCINATION STATUS of j)
    exp[[i]]$vaxdate <- pos$vaxdate[j]
    if (exp[[i]]$date<exp[[i]]$vaxdate){
      exp[[i]]$vax<-0
    } else {exp[[i]]$vax<-1}

  # collaspe 1st to ith exp into the jth element in ondate.hh.exp
  ondate.hh.exp[[j]] <- do.call(rbind,exp)}
  }
}

#------------------------------------------------------------------------------------------

# Here we format the data for the contribution to the hazard function of time 
# at risk of community (d(j)) infection on the date of illness onset of infected j's

# create an empty list to write in d(j) of community exposure of each j 
ondate.c.exp <- vector("list", nrow(pos)) 

# loop through each j
for ( j in 1:nrow(pos)) {  
  # create empty vectors to keep track of current j's ID and hhID
  curhhID <-curID<- temp  <-NA
  # empty vector for j's date of illness onset 
  curondate <- as.Date(NA)
  # write hhID and ID of the current j
  curhhID <-pos$hhID[j] ; curID <-pos$ID[j] 
  # write j's date of illness onset
  curondate <-pos$ondate[j]
  
  # Write data for com exposure, d(j), on the day of illness onset of j.
  # There is 1 element in the list per j, and 1 row in each element for exposure on the
  # date of illness onset, d(j).
  # Data include d(j), and all covariates to be included in regression model,
  # here age, high risk status, and time-varying vaccination status.
  ondate.c.exp[[j]] <-data.frame(date=curondate, AGECAT2=pos$AGECAT2[j], 
                                 AGECAT3=pos$AGECAT3[j], HR=pos$HR[j], 
                                 vaxdate=pos$vaxdate[j])
  for (i in 1:length(ondate.c.exp[[j]][,1])){
    ondate.c.exp[[j]]$com[i]<-com$COM[which(ondate.c.exp[[j]]$date[i]==com$DAY)]
    if (ondate.c.exp[[j]]$vaxdate[i]>ondate.c.exp[[j]]$date[i]) {
      ondate.c.exp[[j]]$vax[i]=0
    } else {ondate.c.exp[[j]]$vax[i]=1}
  }
}
#------------------------------------------------------------------------------------------
