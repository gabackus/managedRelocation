##build a data set of the different types of corals (structure and stratergy for growth)
##I have removed "Folio" because it isnt in the Darling et al Eco let paper which can allow us to pramtererise this table with probabilities
##also changed "massive" to "domed" in line with their useage
## reading the description of the speacies in the paper it seems like brooding vs spawing is a big deal, do we definatly want to not include this?
coral.types<-data.frame(expand.grid("Structure"=c("Domed", "Branching", "Plating"), "Stratergy"=c("Competitive", "Weedy","Stress tolerant", "Generalist")))
coral.types$unique.ID<-paste(coral.types$Structure, coral.types$Stratergy, sep="+")  
coral.types$probability.occurance<-c(0, 95.7, 39.1, 56.5, 39.1, 21.7, 100, 0, 1.7, 100, 35.7, 100)
##not quite sure what these proabbilities refer to at the moment, so I have normalised by stratergy them and used them to draw the random structure:
coral.types$norm.prob<-NA
coral.types$norm.prob[which(coral.types$Stratergy=="Competitive")]<-coral.types$probability.occurance[which(coral.types$Stratergy=="Competitive")]/sum(coral.types$probability.occurance[which(coral.types$Stratergy=="Competitive")])
coral.types$norm.prob[which(coral.types$Stratergy=="Weedy")]<-coral.types$probability.occurance[which(coral.types$Stratergy=="Weedy")]/sum(coral.types$probability.occurance[which(coral.types$Stratergy=="Weedy")])
coral.types$norm.prob[which(coral.types$Stratergy=="Stress tolerant")]<-coral.types$probability.occurance[which(coral.types$Stratergy=="Stress tolerant")]/sum(coral.types$probability.occurance[which(coral.types$Stratergy=="Stress tolerant")])
coral.types$norm.prob[which(coral.types$Stratergy=="Generalist")]<-coral.types$probability.occurance[which(coral.types$Stratergy=="Generalist")]/sum(coral.types$probability.occurance[which(coral.types$Stratergy=="Generalist")])

##add in the generation time from the eco let paper, and its SD. We can use this later to draw random generation times for each individual generated
coral.types$gen.time<-NA
coral.types$gen.time.sd<-NA
coral.types$gen.time[which(coral.types$Stratergy=="Competitive")]<-9.72
coral.types$gen.time.sd[which(coral.types$Stratergy=="Competitive")]<-1.11
coral.types$gen.time[which(coral.types$Stratergy=="Weedy")]<-9.41
coral.types$gen.time.sd[which(coral.types$Stratergy=="Weedy")]<-1.56
coral.types$gen.time[which(coral.types$Stratergy=="Generalist")]<-10
coral.types$gen.time.sd[which(coral.types$Stratergy=="Generalist")]<-0
coral.types$gen.time[which(coral.types$Stratergy=="Stress tolerant")]<-10.17
coral.types$gen.time.sd[which(coral.types$Stratergy=="Stress tolerant")]<-1.52
##if we add in the mean size, and sd size for each, as well as fecundity and sd fecundity then we can assign
##random "sizes" and thus fecundities for each coral after the 1st year of "growth"
##note that size needs to be bounded at the minimum (i.e. we cant have a negative size, or indeed no size). Minimum threshold value?
#thus:
coral.types$size<-NA
coral.types$size.sd<-NA
coral.types$size[which(coral.types$Stratergy=="Competitive")]<-259.32
coral.types$size.sd[which(coral.types$Stratergy=="Competitive")]<-314.68
coral.types$size[which(coral.types$Stratergy=="Weedy")]<-105.54
coral.types$size.sd[which(coral.types$Stratergy=="Weedy")]<-138.54
coral.types$size[which(coral.types$Stratergy=="Generalist")]<-248.77
coral.types$size.sd[which(coral.types$Stratergy=="Generalist")]<-243.39
coral.types$size[which(coral.types$Stratergy=="Stress tolerant")]<-137.81
coral.types$size.sd[which(coral.types$Stratergy=="Stress tolerant")]<-285.26
##Fecundity This is measured in eggs per polyp. I dont know what the relationship between size and number of polyps is though? 
##We need that to link size to reproduction. We could also consider specifically modelling growth, given that we have mean growth rates per year, 
##although we would need to explicity link this to temperature which could be a pain...
coral.types$fecundity<-NA
coral.types$fecundity.sd<-NA
coral.types$fecundity[which(coral.types$Stratergy=="Competitive")]<-20.04
coral.types$fecundity.sd[which(coral.types$Stratergy=="Competitive")]<-32.39
coral.types$fecundity[which(coral.types$Stratergy=="Weedy")]<-19.70
coral.types$fecundity.sd[which(coral.types$Stratergy=="Weedy")]<-14.26
coral.types$fecundity[which(coral.types$Stratergy=="Generalist")]<-18.46
coral.types$fecundity.sd[which(coral.types$Stratergy=="Generalist")]<-9.26
coral.types$fecundity[which(coral.types$Stratergy=="Stress tolerant")]<-372.45
coral.types$fecundity.sd[which(coral.types$Stratergy=="Stress tolerant")]<-694.15
##also add in growth rates, incase we want to make it growth dependant too...
coral.types$growth<-NA
coral.types$growth.sd<-NA
coral.types$growth[which(coral.types$Stratergy=="Competitive")]<-49.83
coral.types$growth.sd[which(coral.types$Stratergy=="Competitive")]<-41.28
coral.types$growth[which(coral.types$Stratergy=="Weedy")]<-11.35
coral.types$growth.sd[which(coral.types$Stratergy=="Weedy")]<-10.04
coral.types$growth[which(coral.types$Stratergy=="Generalist")]<-19.18
coral.types$growth.sd[which(coral.types$Stratergy=="Generalist")]<-11.58
coral.types$growth[which(coral.types$Stratergy=="Stress tolerant")]<-7.98
coral.types$growth.sd[which(coral.types$Stratergy=="Stress tolerant")]<-6.69
##then make it into cm (as the sizes are)
coral.types$growth<-coral.types$growth/10
coral.types$growth.sd<-coral.types$growth.sd/10