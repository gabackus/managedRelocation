rm(list=ls())
library(ggplot2)
library(reshape2)
require(RCurl)
#set.seed(32)
## Basicl model structure

##set some intial paramteres
  ##length of time series (years)
  simulation.length <- 200
  ##the size of the reef (number of niches wide)
  reef.size <- 15
  ##gb: needed for new mortality calculations
  locations <- 1:reef.size


##source in the generic coral data set (from eco let paper). 
  ##needs to link to the github repository, but I cant do this at the moment...
  #source(getURL("https://raw.githubusercontent.com/gabackus/managedRelocation/master/Coral-tactical/Generic%20coral%20data.R"), encoding = 'UTF-8')
  #source('~/Desktop/Work/managedRelocation/Coral-tactical/Generic coral data.R', encoding = 'UTF-8')
  ##gb: This might work a little bit better
  cor.dat <- getURL("https://raw.githubusercontent.com/gabackus/managedRelocation/master/Coral-tactical/Generic%20coral%20data.R", ssl.verifypeer = FALSE)
  eval(parse(text = cor.dat))

  ##look at the dat:
  coral.types
  coral.types$numeric.code<-as.numeric(row.names(coral.types))

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
  #ggplot(data=data.frame(x=1:length(initial.spp.structure), initial.spp.structure), aes(x=x))+geom_bar(aes(fill=initial.spp.structure), width=1)+theme_bw()
##and age structure
  #hist(initial.age.structure)

##create blank matrices for simulation results
results.str<-matrix(NA, nrow=simulation.length, ncol=reef.size)
results.age<-matrix(NA, nrow=simulation.length, ncol=reef.size)
results.size<-matrix(NA, nrow=simulation.length, ncol=reef.size)
results.gen<-matrix(NA, nrow=simulation.length, ncol=reef.size)

## gb: Each individual adult will have a some chance of dying each time step
adult.mortality<-0.1
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
  #for(i in 2:47){
  #i=48
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
    #mort.locations<-sample(1:reef.size, rpois(1, 5))
    ## gb: To make mortality scale better with the reef size, we can give each individual some chance of dying
    ## gb: For now, we can take a uniform random number between 0 and 1 for each space.
    ## gb: Every time that random number is less than adult.mortality, the individual dies.
    mort.locations<-locations[runif(reef.size)<adult.mortality]
    spp.structure.t[mort.locations]<-NA
    age.structure.t[mort.locations]<-NA
    size.structure.t[mort.locations]<-NA
    gen.time.t[mort.locations]<-NA
    
    ##mort.locations needs to be kep for the rest of the code!
  ################################################################################################################################
  ##then those that are left reproduce:
  ##needs to be temp dependant? or is this just acccounted for in mortality in the adult and polyp stages?
    
  ##at the moment this is just an amount drawn from a log normal distribution, roughly based on coral.types
  ##gb: Coral size doesn't do much with the current model. For simplicity, I think we should multiply the fecundity by the coral's size. 
    fec.t<-coral.types$fecundity[match(spp.structure.t, coral.types$unique.ID)]*size.structure.t
  ##if there are individuals that are old enough to reproduce:
  if(length(which(!is.na(fec.t)))>0){
    fec.sd.t<-coral.types$fecundity.sd[match(spp.structure.t, coral.types$unique.ID)]
    ##only keep the ones that are old enough to reproduce  
    fec.t[which(age.structure.t<gen.time.t)]<-NA
    ##if there are indivudals old enough to reproduce:
    if(length(which(!is.na(fec.t)))>0){
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
      free.niches.t<-length(mort.locations)
      
    ##if there are enough survivng polyps to repopulation all the free niches:
      if(sum(larvae.t$surviving)>=free.niches.t){
      ##create a vector of larvae to settle, and randomly sample from this, add this into the reef (all vectorized):
      spp.structure.t[mort.locations]<-as.character(sample(rep(larvae.t$ID, larvae.t$surviving), size=free.niches.t, replace=F))
    }else{
      ##if there arent enough surviing offspring to settle, then doa bit of a hack around.
      ##calculate the number of larvae free to settle
      free.larvae<-rep(larvae.t$ID, larvae.t$surviving)
      if(length(free.larvae)>0){
        ##create a vector of NA's to fill
        na.free.niches<-rep(NA, length=free.niches.t)
        ##fill as many as there are free polyps, leaving some NA's left over
        na.free.niches[1:length(rep(larvae.t$ID, larvae.t$surviving))]<-as.character(rep(larvae.t$ID, larvae.t$surviving))
        ##randomize the order of these so that the settlement in space is randome
        spp.structure.t[mort.locations]<-sample(na.free.niches, length(na.free.niches), replace=F)
      }
    }
    
    ##update the age structure to have those nices that are filled (which will be all, is that ok?) as age 0:
    age.structure.t[mort.locations]<-0
    ################################################################################################################################
    ##then insert some size structure for the next time point, reall this should be a growth function based on the per year growth 
    ##rates in coral.types...
    ##add in the zeroes (just settled corals)
    size.structure.t[mort.locations]<-0
  ##close the if statement about reproduction
    }
  }
  
  ##then add on growh for each type based on the table 
  growth.t<-coral.types$growth[match(spp.structure.t, coral.types$unique.ID)]
  growth.sd.t<-coral.types$growth.sd[match(spp.structure.t, coral.types$unique.ID)]
  grow.location<-log(growth.t^2 / sqrt(growth.sd.t^2 + growth.t^2))
  grow.shape <- sqrt(log(1 + (growth.sd.t^2 / growth.t^2)))
  size.structure.t<-size.structure.t+rlnorm(growth.t, meanlog = grow.location, sdlog = grow.shape)
  ################################################################################################################################
  ##generate a generation time for the new individuals and add it into the results
  gen.time.t[mort.locations]<-round(rnorm(length(coral.types$gen.time[match(spp.structure.t[mort.locations], coral.types$unique.ID)]), coral.types$gen.time[match(spp.structure.t[mort.locations], coral.types$unique.ID)], coral.types$gen.time.sd[match(spp.structure.t[mort.locations], coral.types$unique.ID)]))

  ################################################################################################################################
  ##save out results as two data frames:
    ##1) age structure at each time step
    ##2) spp structure at each time step
  ##saved as rows in a maxtrix (created above), coded to be used in a loop with year=i
  size.structure.t[which(is.nan(size.structure.t))]<-NA
  age.structure.t[which(is.nan(age.structure.t))]<-NA
  gen.time.t[which(is.nan(gen.time.t))]<-NA
  
  results.str[i,]<-spp.structure.t
  results.age[i,]<-age.structure.t
  results.size[i,]<-size.structure.t
  results.gen[i,]<-gen.time.t
  ################################################################################################################################
  if(length(which(!is.na(spp.structure.t)))==0){break("all dead")}
  ################################################################################################################################
}##end loop 

##gb: I added some stuff to help me visualize the model a little bit better
##This could probably look a little cleaner. I'm not completely fluent in R yet.

##This matrix is to convert the results to population size
##We might want to convert this to proportions (by simply dividing by the reef size)
##or think about other appropriate diversity indices (richness, inverse Simpson's, etc.)
##There is a column for each unique coral type and a row for each time step
results.pop.size<-matrix(0,simulation.length,length(coral.types$unique.ID))

##This converts results.str to population sizes
for (j in 1:length(coral.types$unique.ID)){
  results.pop.size[,j]<-rowSums(results.str==coral.types$unique.ID[j])
}

##This makes it a dataframe
pop.size.df<-as.data.frame(results.pop.size)
##Names the columns
colnames(pop.size.df)<-coral.types$unique.ID
##Adds time as a column
pop.size.df$Time=1:simulation.length

##I ended up plotting it as proportions
ggplot((melt(pop.size.df,id='Time')),aes(x=Time,y=value/reef.size,fill=variable))+
  geom_area(colour="black", size=2, alpha=.4) +
  scale_fill_manual(values=rainbow(12))+
  ylab('Proportion of reef')+theme_bw()

##the other option is as an image:
##convert to numeric
num.str<-matrix(coral.types$numeric.code[match(results.str, coral.types$unique.ID)],nrow=simulation.length, ncol=reef.size) 
image(num.str)
