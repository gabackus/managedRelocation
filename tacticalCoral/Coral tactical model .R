rm(list=ls())
library(ggplot2)
library(reshape2)
require(RCurl)
## Basicl model structure

##set some intial paramteres
  ##length of time series (years)
  simulation.length <- 100
  ##the size of the reef (number of niches wide)
  reef.size <- 20

##source in the generic coral data set (from eco let paper)
  source(getURL("https://raw.githubusercontent.com/chrit88/Changes-bowie/master/bowie%204.txt"))


##simulate an initial structure for the reef (both age and spatial). Current set to randomly select from a uniform distribution for structure
##and from a gamma distrubition for age rounded to whole numbers 
##(we should estimate these from data by fitting a gamm distribution to the age data from a reef, if possible)
initial.spp.structure <- sample(coral.types$unique.ID, size=reef.size, replace=T)
initial.age.structure <- round(rgamma(reef.size, shape=1, scale=40))

##plot out the starting structure of the reef
ggplot(data=data.frame(x=1:length(initial.spp.structure), initial.spp.structure), aes(x=x))+geom_bar(aes(fill=initial.spp.structure), width=1)+theme_bw()
##and age structure
hist(initial.age.structure)

##create blank matrices for simulation results
results.str<-matrix(NA, nrow=simulation.length, ncol=reef.size)
results.age<-matrix(NA, nrow=simulation.length, ncol=reef.size)

##the simulation:

### then at each time step:

##spp.structure.t will be the object that is manipulated in each simulation, then saved out at the end
  spp.structure.t<-initial.spp.structure
  age.structure.t<-initial.age.structure
  ##for the moment lets add in some empty niches so we can code:
  spp.structure.t[c(2, 12, 16)]<-NA
  age.structure.t[c(2, 12, 16)]<-NA

################################################################################################################################    
## first the population ages:
age.structure.t<-age.structure.t+1
################################################################################################################################  
##then we have adult mortality (needs to be temp dependant):

  
################################################################################################################################
##then those that are left reproduce:
##needs to be temp dependant
##at the moment this is just a random amount drawn from a unifrom distribution
larvae.t<-data.frame(ID=coral.types$unique.ID, number=round(runif(length(coral.types$unique.ID),0, 100000)))

################################################################################################################################
##then mortality in the larval stages (probably just uniform mortality (e.g. 0.01% survive)). This will introduce allee effects
##needs to be temperature dependant 
##currently just fixed to 0.1%
mort<-0.001

##calculate the remining larvae availble to settle
larvae.t$surviving<-round(larvae.t$number*mort)

################################################################################################################################
##then settlement (not temperature dependant):

##number of vacant niches:
free.niches.t<-length(which(is.na(spp.structure.t)))

##create a vector of larvae to settle, and randomly sample from this, add this into the reef (all vectorized):
spp.structure.t[which(is.na(spp.structure.t))]<-as.character(sample(rep(larvae.t$ID, larvae.t$surviving), size=free.niches.t, replace=T))

##update the age structure to have those nices that are filled (which will be all, is that ok?) as age 0:
age.structure.t[which(is.na(age.structure.t))]<-0
################################################################################################################################
##save out results as two data frames:
  ##1) age structure at each time step
  ##2) spp structure at each time step
##saved as rows in a maxtrix (created above), coded to be used in a loop with year=i

results.str[i,]<-spp.structure.t
results.age[i,]<-age.structure.t
################################################################################################################################
#}##end loop