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

##source in the generic coral data set (from eco let paper). 
  ##needs to link to the github repository, but I cant do this at the moment...
  ##source(getURL("https://raw.githubusercontent.com/chrit88/Changes-bowie/master/bowie%204.txt"))
  source('~/Desktop/Work/managedRelocation/tacticalCoral/Generic coral data.R', encoding = 'UTF-8')
  ##look at the dat:
  coral.types

##simulate an initial structure for the reef (both age and spatial). Current set to randomly select (by weighting)
##and from a gamma distrubition for age rounded to whole numbers 
##(we should estimate these from data by fitting a gamm distribution to the age data from a reef, if possible)
initial.spp.structure <- sample(coral.types$unique.ID, prob=coral.types$norm.prob, size=reef.size, replace=T)
initial.age.structure <- round(rgamma(reef.size, shape=1, scale=40))

##lognormal distribution for drawing sizes from...or just interpolate from the ages? might make more sense...
##this is simpler for the moment so lets stick with this:
m<-coral.types$size[match(initial.spp.structure, coral.types$unique.ID)]
s<-coral.types$size.sd[match(initial.spp.structure, coral.types$unique.ID)]
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
initial.size.structure<-rlnorm(coral.types$size[match(initial.spp.structure, coral.types$unique.ID)], meanlog = location, sdlog = shape)

##make a vector of the generation times of each of the "species" in the community, based on data in coral.types:
##drawn from normal distribution with mean and sd set via the coral.types table
##rounded to whole numbers for simulations
initial.gen.time<-round(rnorm(coral.types$gen.time[match(initial.spp.structure, coral.types$unique.ID)], coral.types$gen.time[match(initial.spp.structure, coral.types$unique.ID)], coral.types$gen.time.sd[match(initial.spp.structure, coral.types$unique.ID)]))

################################################################################################################################
##the inital vectors:
initial.spp.structure
initial.age.structure
initial.size.structure
initial.gen.time
    
##plot out the starting structure of the reef
ggplot(data=data.frame(x=1:length(initial.spp.structure), initial.spp.structure), aes(x=x))+geom_bar(aes(fill=initial.spp.structure), width=1)+theme_bw()
##and age structure
#hist(initial.age.structure)

##create blank matrices for simulation results
results.str<-matrix(NA, nrow=simulation.length, ncol=reef.size)
results.age<-matrix(NA, nrow=simulation.length, ncol=reef.size)
results.size<-matrix(NA, nrow=simulation.length, ncol=reef.size)
results.gen<-matrix(NA, nrow=simulation.length, ncol=reef.size)


################################################################################################################################
##the simulation:

### then at each time step:

##spp.structure.t will be the object that is manipulated in each simulation, then saved out at the end
##initial for the first round of added to the results and read back in from there on, hence 2:simulation.length
  results.str[1,]<-initial.spp.structure
  results.age[1,]<-initial.age.structure
  results.size[1,]<-initial.size.structure
  results.gen[1,]<-initial.gen.time
  
for(i in 2:simulation.length){
#i=2
################################################################################################################################    
##the structures from the previous year, these will then get updated in the loop:
spp.structure.t<-results.str[i-1,]
age.structure.t<-results.age[i-1,]
size.structure.t<-results.size[i-1,]
gen.time.t<-results.gen[i-1,]
################################################################################################################################  
## first the population ages:
age.structure.t<-age.structure.t+1
################################################################################################################################  
##then we have adult mortality (needs to be temp dependant):
##currently just trandom  
  # # Some thermal performance data
  # temp = c(17, 21, 24, 28, 31, 33, 
  #          17, 21, 24, 28, 31, 33)
  # surv = c(0.3, 0.4, 0.68, 0.82, 0.78, 0.3, 
  #            0.4, 0.45, 0.58, 0.75, 0.83, 0.6)
  # library(mgcv)
  # surv.temp<-gam(surv~temp, family=binomial)  
  # new.temps<-seq(10, 40, 0.5)
  # surv.temp<-data.frame(temps=new.temps, surv=predict(surv.temp, newdata=data.frame(temp=new.temps), type="response"))
  # 
  ##for the moment lets add in some empty niches so we can code:
  ##the lcoations in space of thos that have died
  mort.locations<-sample(1:reef.size, 3)
  spp.structure.t[mort.locations]<-NA
  age.structure.t[mort.locations]<-NA
  size.structure.t[mort.locations]<-NA
  gen.time.t[mort.locations]<-NA
  
  ##mort.locations needs to be kep for the rest of the code!
################################################################################################################################
##then those that are left reproduce:
##needs to be temp dependant? or is this just acccounted for in mortality in the adult and polyp stages?
  
##at the moment this is just an amount drawn from a log normal distribution, roughly based on coral.types
fec.t<-coral.types$fecundity[match(spp.structure.t, coral.types$unique.ID)]
fec.sd.t<-coral.types$fecundity.sd[match(spp.structure.t, coral.types$unique.ID)]
##only keep the ones that are old enough to reproduce  
fec.t[which(age.structure.t<gen.time.t)]<-NA
fec.sd.t[which(age.structure.t<gen.time.t)]<-NA
fec.location <- log(fec.t^2 / sqrt(fec.sd.t^2 + fec.t^2))
fec.shape <- sqrt(log(1 + (fec.sd.t^2 / fec.t^2)))
larvae.t<-na.omit(data.frame(ID=spp.structure.t, number=rlnorm(fec.t, meanlog = fec.location, sdlog = fec.shape)))
################################################################################################################################
##then mortality in the larval stages (probably just uniform mortality (e.g. 0.01% survive)). This will introduce allee effects
##needs to be temperature dependant 
##currently just fixed to 50%
mort<-0.5

##calculate the remining larvae availble to settle
larvae.t$surviving<-round(larvae.t$number*mort)

################################################################################################################################
##then settlement (not temperature dependant):
##number of vacant niches:
free.niches.t<-length(which(is.na(spp.structure.t)))

##create a vector of larvae to settle, and randomly sample from this, add this into the reef (all vectorized):
spp.structure.t[mort.locations]<-as.character(sample(rep(larvae.t$ID, larvae.t$surviving), size=free.niches.t, replace=F))

##update the age structure to have those nices that are filled (which will be all, is that ok?) as age 0:
age.structure.t[mort.locations]<-0

################################################################################################################################
##then insert some size structure for the next time point, reall this should be a growth function based on the per year growth 
##rates in coral.types...
##add in the zeroes (just settled corals)
size.structure.t[mort.locations]<-0
##then add on growh for each type based on the table 

growth.t<-coral.types$growth[match(spp.structure.t, coral.types$unique.ID)]
growth.sd.t<-coral.types$growth.sd[match(spp.structure.t, coral.types$unique.ID)]
grow.location<-log(growth.t^2 / sqrt(growth.sd.t^2 + growth.t^2))
grow.shape <- sqrt(log(1 + (growth.sd.t^2 / growth.t^2)))
size.structure.t<-size.structure.t+rlnorm(growth.t, meanlog = grow.location, sdlog = grow.shape)
################################################################################################################################
##generate a generation time for the new individuals and add it into the results
gen.time.t[mort.locations]<-round(rnorm(coral.types$gen.time[match(spp.structure.t[mort.locations], coral.types$unique.ID)], coral.types$gen.time[match(spp.structure.t[mort.locations], coral.types$unique.ID)], coral.types$gen.time.sd[match(spp.structure.t[mort.locations], coral.types$unique.ID)]))
################################################################################################################################
##save out results as two data frames:
  ##1) age structure at each time step
  ##2) spp structure at each time step
##saved as rows in a maxtrix (created above), coded to be used in a loop with year=i

results.str[i,]<-spp.structure.t
results.age[i,]<-age.structure.t
results.size[i,]<-size.structure.t
results.gen[i,]<-gen.time.t
################################################################################################################################
}##end loop