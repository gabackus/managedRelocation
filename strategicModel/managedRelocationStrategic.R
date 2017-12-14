library(pracma)
library(rmutil)
library(grid)
library(gridExtra)
library(gtable)
library(ggplot2)
library(reshape)

amSetup <- function(S=50,XL=2^9,tempInc=0.04,zVals='random',gamMean=0.005,gamSD=0){
  # To setup the basic model of architecture of the assisted migration model, run this function

  # First, specify parameter values for all of the S species.
  # Optimal temperature for each species. These can be randomly picked from a uiform distribution or evenly spread over the range of acceptable values.
  if (zVals=='random'){
    z <- runif(S,10,30)
  } else if (zVals=='even'){
    z <- seq(10,30,S+2) # Make a sequence from the lowest temperature to the highest temperature, with extra space at either end.
    z <- z[2:(S+1)]     # Remove the padding on either side so there are S optimal temperatures.
  }
  
  # Mean dispersal distance for each species. These are randomly (lognormal) distributed around gamMean with standard deviation gamSD.
  # See Urban et al. 2012 for similar methodology and reasoning.
  gam <- rlnorm(S,log(gamMean),log(1+gamSD))

  # sig is the width of thermal tolerance. Roughly, the species has positive growth for temperatures between z-sig and z+sig, though this is skewed left.
  sig <- array(4,S)
  # lam is the skewness in thermal tolerance. For simplicity, we stick with -3.4 (as in Urban et al. 2012) to skew the thermal tolerance to have a mild decrease moving toward colder temperatures and a sharp decrease moving toward warmer temperatures.
  lam <- array(-3.4,S)
  # delt is the width of the competitive thermal tolerance. Roughly, analogous to sig but for determining competition coefficients
  delt <- array(4,S)
  # a is the relative strength of competitive ability in each species.
  # This differs from delt because a species might be a strong competitor (high a) over a small range (low delt) or a weak competitor (low a) over a high range (high delt).
  a <- array(1,S)

  
  # Next, we specify some of the model architecture
  # Spatial domain: from -1 to 1 over XL nodes (I like to use powers of 2 for XL)
  # The distance between nodes is the length of the domain divided by the number of nodes
  dx <- 2/XL
  # Building the spatial domain nodes (we end up cutting off -1)
  x <- seq(-1,1,length=XL+1)
  x <- x[-1]
  
  # When dispersing over a space into bins, it is helpful to have a vector of where to cutoff for these bins.
  xB <- seq(x[2]-dx/2,x[XL-2]+dx/2,length=XL-2)
  # Indexing the dispersal of the offspring if easier with these parameters
  kL <- c(length(xB)-1,length(xB)/2)
  
  # Temperature: linear increase from 10 (at x=-1) to 30 (at x=1)
  temp <- seq(10,30,length=XL)
  # Temperature increase rate is tempInc (constant for now)
  
  # Habitat quality: constant over space (for now)
  Q <- array(50,XL)
  
  # Matrix versions of certain values for easier matrix arithmetic
  tempM <- t(kronecker(matrix(1,1,S),temp))
  zM <- kronecker(matrix(1,1,XL),z)
  sigM <- kronecker(matrix(1,1,XL),sig)
  lamM <- kronecker(matrix(1,1,XL),lam)
  deltM <- kronecker(matrix(1,1,XL),delt)
  aM <- kronecker(matrix(1,1,XL),a)
  xM <- t(kronecker(matrix(1,1,S),x))
  # Reproductive rates are adjusted so minimum is 0 and maximum is 2

  r <- exp(-(tempM-zM)^2/sigM^2)*(1+erf(lamM*(tempM-zM)/sigM))
  rMin <- min(r)
  rMax <- max(r)
  # Adjusted reproductive rate is calculated as C[1]+C[2]*(unadjusted rate)
  const <- c(-2*rMin/(rMax-rMin),2/(rMax-rMin))
  
  # Export it all as a list
  return(list(P=list(S=S,           # Biological parameters
                     z=z,
                     gam=gam,
                     sig=sig,
                     lam=lam,
                     delt=delt,
                     a=a,
                     const=const),
              X=list(XL=XL,         # Abiotic parameters and model architecture
                     dx=dx,
                     x=x,
                     xB=xB,
                     kL=kL,
                     temp=temp,
                     tempInc=tempInc,
                     Q=Q),
              M=list(tempM=tempM,   # Matricized values for computation
                     zM=zM,
                     sigM=sigM,
                     lamM=lamM,
                     deltM=deltM,
                     aM=aM,
                     xM=xM)))

}

amInitialize <- function(P,X,M,iTime=10,lcOrder=list(lc=1,mr=0,cl=0)){
  # To find the initial distribution of the population sizes, we run the model with constant temperature for iTime years.
  # Some species might go extinct in the process. Those are eliminated from the model for subsequent simulations.
  
  # Initially set all species to be at half the maximum carrying capacity. This is fairly arbitrary, as the model will initialize long enough that the initial seed should not matter.
  n <- matrix(max(X$Q)/2,P$S,X$XL)

  for (i in 1:iTime){
    n<-amTimeStep(n,P,X,M,0,lcOrder)
  }

  # Calculate the total population size to find which species are extinct
  N <- rowSums(n)
  # Identify which species are still extant after initializing time steps
  notDead <- which(N>0)
  
  # Update n and the P and M lists to incorporate this change
  n <- n[notDead,]
  
  P=list(S=length(notDead),
         z=P$z[notDead],
         gam=P$gam[notDead],
         sig=P$sig[notDead],
         lam=P$lam[notDead],
         delt=P$delt[notDead],
         a=P$a[notDead],
         const=P$const)
  M=list(tempM=t(kronecker(matrix(1,1,P$S),X$temp)),
         zM=kronecker(matrix(1,1,X$XL),P$z),
         sigM=kronecker(matrix(1,1,X$XL),P$sig),
         lamM=kronecker(matrix(1,1,X$XL),P$lam),
         deltM=kronecker(matrix(1,1,X$XL),P$delt),
         aM=kronecker(matrix(1,1,X$XL),P$a),
         xM=t(kronecker(matrix(1,1,P$S),X$x)))
  
  return(list(n=n,
              P=P,
              M=M))
}

amSimulate <- function(n,P,X,M,V=0,lcOrder=list(lc=1,mr=0,cl=0),ccYears=100,plotYears=0){
  # This simulates an initialized community, n, over ccYears
  # The function returns the final population distribution as well as the population size and center over time
  # If you want to see how the population distribution changes over time, set the plotYears to a positive integer
    
  # First, we set up a matrix to save the total population size over time
  N <- matrix(0,P$S,ccYears+1)
  # Record the initial population size
  N[,1] <- rowSums(n)
  
  # We also save the weighted center of the population
  xC <- matrix(0,P$S,ccYears+1)
  # Record the intiial population center
  xC[,1] <- rowSums(n*M$xM)/rowSums(n)
  
  # Run the model over ccYears
  for (i in 1:ccYears){
    
    # Temperature changes before each time step
    X$temp=X$temp+X$tempInc
    M$tempM=M$tempM+X$tempInc
    
    # Run the time step after the change in temperature
    n<-amTimeStep(n,P,X,M,V,lcOrder)
    
    # Record the population size and center
    N[,i+1]<-rowSums(n)
    xC[,i+1] <- rowSums(n*M$xM)/rowSums(n)
    
    # If it's a plotting year, show the population distribution
    if (plotYears!=0){
      if (i%%plotYears==0){
        plotCom(n,P,X)
      }
    }
  }
  return(list(n=n,N=N,xC=xC))
}

compCoef <- function(P,X,M){
  # For every pair of species and point in space, this function calculates the competition coefficient at the current time step with an output of an SxSxXL matrix
  # These are similar to Lotka-Volterra competition coefficients
  if (sum(P$a)==0){                       # If all a values are set to 0, there is no interspecific competition
    alph<-array(0,dim=c(P$S,P$S,X$XL))
    for (i in 1:iTime){
      alph[i,i,]=1                        # All coefficients are 0, except for intraspecific competition, which is 1
    }
  } else{                                 # Otherwise, the alph values are a ratio of a function that is similar to the reproductive function
    alphi<-M$aM*exp(-(M$tempM-M$zM)^2/M$delt^2)*(1+erf(M$lamM*(M$tempM-M$zM)/M$deltM))+1 # Individual competition functions
    alph<-array(0,dim=c(P$S,P$S,X$XL))
    for (i in 1:P$S){
      for (j in 1:P$S){
        alph[i,j,]=alphi[j,]/alphi[i,]            # The ratio of these individual competition functions is the competition coefficients
      }
    }
  }
  return(alph)
}

reproduce <- function(n,P,X,M){
  # The number of offspring born for each species in each location is a Poisson random variable with mean r*n
  
  # The base reproductive rate is a skewed function, adjust such that min(r)=0 and max(r)=2
  # Each species will have a different reproductive rate depending on the temperature at that space.
  r <- P$const[1]+P$const[2]*exp(-(M$tempM-M$zM)^2/M$sigM^2)*(1+erf(M$lamM*(M$tempM-M$zM)/M$sigM))
  
  # For each species and location, the random number of offspring are born
  for (i in 1:P$S){
    for (j in 1:X$XL){
      if (n[i,j]>0){
        n[i,j]=rpois(1,r[i,j]*n[i,j])
      } else{
        n[i,j]=0
      }
    }
  }
  return(n)
}

dispersal <- function(n,P,X,M){
  # Each individual spreads throughout the spatial landscape with a random Laplacian dispersal kernel determined by the species' mean dispersal distance, gam[i].

  # Preallocate the output array
  nD=array(0,c(P$S,X$XL))

  # For each species in each location, disperse!
  for (i in 1:P$S){
    for (j in 1:X$XL){
      if (n[i,j]>0){
        # Each individual lands a certain distance away from its origin, based on a Laplace random variable
        spread=rlaplace(n[i,j],0,P$gam[i])
        # Because the Laplace distribution is continuous, we discretize the output by placing each individual into a bin, based on its distance away from the origin
        spreadC=histc(spread,X$xB)
        # Indexing gets a a little complicated, as some individuals will disperse beyond the boundary of the spatial domain
        # Unlike the base Urban et al. 2012 model, this model does not have reflecting boundary conditions
        # Instead, those that spread beyond the spatial boundary are lost from the system
        # The indexing helps with shifting dispersal shadow and varying window sizes
        jSi = max(1,X$kL[2]+1-j); jSf = min(X$kL[1],X$XL+X$kL[2]-j);
        jBi = max(1,j-X$kL[2]+1); jBf = min(X$XL,X$kL[2]-1+j);
        
        # Add the individuals to the location that the arrived at
        nD[i,jBi:jBf]=nD[i,jBi:jBf]+spreadC$cnt[jSi:jSf]
      }
    }
  }
  return(nD)
}

densDepen <- function(n,P,X,M){
  # The density dependence in this model is roughly a Beverton-Holt model that includes both interspecific and intraspecific competition
  # Each individual has a random chance of survival based on a variety of conditions
  
  # Preallocate the output array
  nDD=array(0,c(P$S,X$XL))

  # Competition coefficients depend on interactions between each species and the temperature at the location at the time
  # These can be thought of as temperature-varying Lotka-Volterra competition coefficients
  alph=compCoef(P,X,M)
  
  # For each species in each location, compete!
  for (i in 1:P$S){
    for (j in 1:X$XL){
      if (n[i,j]>0){
      # Each individual has a mortality chance based on the density dependent competition
        nDD[i,j]=n[i,j]-rbinom(1,n[i,j],1-1/(1+1/X$Q[j]*t(alph[i,,j]%*%n[,j])))
      }
    }
  }
  return(nDD)
}

assistMig <- function(n,P,X,V){}
cullDomin <- function(n,P,X,V){}

amTimeStep <- function(n,P,X,M,V=0,lcOrder=list(lc=1,mr=0,cl=0)){
  # Go through each part of the model in a specified order.
  # Each time step could be one "year" or one "generation", but ultimately it runs through each part of the life cycle in an order determined by lcOrder.
  # Reproduction, dispersal, and density dependence are all required steps in lcOrder
  # The various management techniques can optionally be added to the model between any two of the required steps
  if (lcOrder$lc==1){
    # When lcOrder is 1
    # reproduction -> dispersal -> density dependence
    n <- reproduce(n,P,X,M)
    n <- dispersal(n,P,X,M)
    n <- densDepen(n,P,X,M)
  } else if (lcOrder$lc==2){
    # when lcOrder is 2
    # reproduction -> density dependence -> dispersal
    n <- reproduce(n,P,X,M)
    n <- densDepen(n,P,X,M)
    n <- dispersal(n,P,X,M)
  }
}

plotCom <- function(n,P,X){
  pC <- as.data.frame(t(n))
  pC$x = X$temp
  pCM <- melt(pC,id=c("x"))
  print(ggplot(pCM,aes(x=x,y=value,colour=variable))+
          geom_line(size=2)+
          ylab('Population size')+
          theme_classic()+
          theme(legend.position="none")+
          scale_color_manual(name="variable",values = rainbow(P$S)))
  
}

plotComTime <- function(N,xC,P,X){
  pCT <- as.data.frame(t(log(N+1)))
  pCT$Time = 0:100
  pCTM <- melt(pCT,id=c("Time"))
  p1<-ggplot(pCTM,aes(x=Time,y=value,colour=variable))+
    geom_line(size=2)+
    ylab('Log opulation size')+
    xlab('Time')+
    theme_classic()+
    theme(legend.position="none")+
    scale_color_manual(name="variable",values = rainbow(P$S))
  g1<-ggplotGrob(p1)
  
  pCx <- as.data.frame(t(xC))
  pCx$Time = 0:100
  pCx <- melt(pCx,id=c("Time"))
  p2<-ggplot(pCx,aes(x=Time,y=value,colour=variable))+
      geom_line(size=2)+
      ylab('Population center')+
      xlab('Time')+
      theme_classic()+
      theme(legend.position="none")+
      scale_color_manual(name="variable",values = rainbow(P$S))  
  g2<-ggplotGrob(p2)
  
  grid.arrange(g1,g2,ncol=1)
}

###########################################################################################
# Here's an example simulation

# First, the model parameters and model architecture are created with 50 species.
setup <- amSetup(50)
# The list of parameters (P), the list of space related terms (X), and the matrix version of some important terms (M) are stored in Ps, Xs, and Ms (respectively)
Ps <- setup$P; Xs <- setup$X; Ms <- setup$M

# Initial model conditions are obtained by running the model with stable temperature over a certain number of years.
# The population distributions after the initialization are used as the initial distribution for the simulation.
# Species that go extinct during this initialization step are removed from the model.
# Normally, we'd want to run the initialization step until all species that are likely to go extinct in the absence of temperature change do go extinct.
# To speed up this example, we can stick with an initialization time of 50.

init <- amInitialize(Ps,Xs,Ms,iTime=50)
# We separate the outputs and save them for the simulation.
ni <- init$n; Pi <- init$P; Mi <- init$M; Xi <- Xs

# This simulates the model under 100 years of climate change.
# Additionally, it will plot the population size of every species in the community every 4 years.
# Make plotYears a smaller integer to see plots more frequently or set it to 0 to see no plots at all.
simu <- amSimulate(ni,Pi,Xi,Mi,0,ccYears=100,plotYears=4)
# Separate and outputs to population size over species and space (n) and total population size of each species over time (N).
n <- simu$n; N <- simu$N; xC <- simu$xC

# Plot the community change over time
plotComTime(N,xC,Pi,Xi)

# If we want to compare the same community reacting to a different temperature change scenario, let's change the temperature increase rate
Xi$tempInc <- 0.02
simu2 <- amSimulate(ni,Pi,Xi,Mi,0,ccYears=100,plotYears=4)
# Separate and outputs to population size over species and space (n) and total population size of each species over time (N).
n2 <- simu2$n; N2 <- simu2$N; xC2 <- simu2$xC

# Plot the community change over time
plotComTime(N2,xC2,Pi,Xi)


# Number of species set
Ps$S
# Number of species after initialization
Pi$S
# Number of species after temperature change scenario 1
sum(rowSums(n)>0)
# Number of species after temperature change scenario 2
sum(rowSums(n2)>0)

