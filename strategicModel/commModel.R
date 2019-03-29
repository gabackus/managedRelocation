rm(list=ls())
require(pracma)
require(parallel)

commSetup <- function(S=32, L=512, W=8,
                      zo=NULL, gam=NULL, sig=NULL, A=NULL,
                      gamMean=2.5, gamSD=2.5,
                      sigMean=5, sigSD=5,
                      lam=-2.7, B=10,
                      compType="lottery",
                      XW=seq(129,384),
                      temp2d=NULL,
                      tempLow=9.78, tempHigh=30.22,
                      tempRev=F,
                      tempGH=1, tempGSD=0, tempLSD=1,
                      years=1000,
                      y=1,
                      tempY=NULL,
                      tau=0.04,
                      tempYAC=0.767, tempYSD=0.1639,
                      Q=NULL,
                      QMean=8,
                      QGH=1, QGSD=0, QLSD=0){
  
  # This sets up the basic structure of the model.
  # It creates a list of biological parameters in P (S randomized species and their parameters)
  # and it creates a list of environmental parameters in X
  
  # S:        Total number of species created with hetSetup.
  # L:        Total number of patches in the metacommunity.
  # W:        Number of microhabitats in each patch.
  
  # zo:       A vector of pre-defined optimal temperature values. Only works if length(zo) is S.
  # gam:      A vector of pre-defined mean dispersal distances. Only works if length(gam) is S. If no vector is specified, gamMean and gamSD are used to randomly generate the vector gam.
  # sig:      A vector of pre-defined thermal tolerance breadths. Only works if length(sig) is S. If no vector is specified, sigMean and sigSD are used to randomly generate the vector sig.
  # A:        Matrix of the relative competition coefficients between species.
  # gamMean:  Mean dispersal distance for randomized species. Default is based on Urban et al. 2012.
  # gamSD:    Standard deviation of dispersal distance for randomized species. Default is based on Urban et al. 2012.
  # sigMean:  Mean thermal tolerance breadth for randomized species. Default is based on Urban et al. 2012.
  # sigSD:    Standard deviation of thermal tolerance breadth for randomized species. Default is based on Urban et al. 2012.
  # lam:      Skewness in thermal tolerance. Default is based on Urban et al. 2012. (to have a mild decrease moving toward colder temperatures and a sharp decrease moving toward warmer temperatures).
  # B:        Area of integrated birth rate over all T for each species.
  # compType: The type of competition in a string. Must be either "lottery" or "temp".
  
  # XW:       Window of analysis (to remove edge effects)
  
  # temp2d:   A matrix of pre-defined temperatures over x. Only works if nrow(temp2d) is L and ncol(temp2d) is W.
  # tempLow:  Lowest mean temperature on linear temperature gradient. temp1d(L)=tempLow
  # tempHigh: Highest mean temperature on linear temperature gradient. temp1d(1)=tempHigh
  # tempRev:  If tempRev=T, then temp1d(1)=tempLow and temp1d(L)=tempHigh.
  # tempGH:   Hurst exponent for global temperature heterogeneity. H=0.5 is Brownian motion; 0.5<H<=1 is long-term positive autocorrelation and 0<=H<0.5 is long-term negative autocorrelation.
  # tempGSD:  Controls the magnitude of the global temperature heterogeneity.
  # tempLSD:  Standard deviation in temperature between microhabitats in each patch.
  # years:    Maximum number of years for initialization + climate change.
  # y:        Current year.
  # tempY:    A vector of pre-defined temperatures over time. Only works if length(tempY) is years.
  # tau:      Average temperature change per year.
  # tempYAC:  Temperature autocorrelation over time. Default is global temperature AC from 1880-1979.
  # tempYSD:  Temperature standard deviation over time. Default is global temperature SD from 1880-1979.
  
  # Q:        A matrix of pre-defined habitat quality over x. Only works if nrow(Q) is L and ncol(Q) is W.
  # QMean:    Average habitat quality for any microhabitat
  # QGH:      Hurst exponent for global habitat quality heterogeneity. H=0.5 is Brownian motion; 0.5<H<=1 is long-term positive autocorrelation and 0<=H<0.5 is long-term negative autocorrelation.
  # QGSD:     Controls the magnitude of the global habitat quality heterogeneity.
  # QLSD:     Standard deviation in habitat quality between microhabitats in each patch.
  
  
  
  
  ##########################################################
  # First, specify biological parameter values for all of the S species.
  
  # Optimal temperature for each species. These can be randomly picked from a uniform distribution or pre-defined.
  if(is.null(zo)){
    zo <- runif(S,9.9,30.1)
  } else {
    if(length(zo)!=S){
      stop("zo does not match the number of species!")
    }
  }
  
  # Dispersal distance for each species. These can be randomly picked from a lognormal distribution or pre-defined.
  if(is.null(gam)){
    # To use the lognormal distribution, we need to convert mean and SD values
    gamMu <- log(gamMean/sqrt(1+gamSD^2/gamMean^2))
    gamSig <- sqrt(log(1+gamSD^2/gamMean^2))
    gam <- rlnorm(S,gamMu,gamSig)
  } else {
    if(length(gam)!=S){
      stop("gam does not match the number of species!")
    }
  }
  
  # Thermal tolerance breadth for each species. These can be randomly picked from a lognormal distribution or pre-defined.
  if(is.null(sig)){
    # To use the lognormal distribution, we need to convert mean and SD values
    sigMu <- log(sigMean/sqrt(1+sigSD^2/sigMean^2))
    sigSig <- sqrt(log(1+sigSD^2/sigMean^2))
    sig <- rlnorm(S,sigMu,sigSig)
  } else {
    if(length(sig)!=S){
      stop("sig does not match the number of species!")
    }
  }
  
  # Competition coefficients between each pair of species species. By default, all coefficients are 1 (lottery competition), but A can be pre-defined.
  if(is.null(A)){
    A <- matrix(1,S,S)
  }else{
    if(!(nrow(A)==S & ncol(A)==S)){
      stop("A does not match the number of species!")
    }
  }
  
  ##########################################################
  # Using the randomized species parameters, we derive other variables needed for computation.
  
  # zo helps define where the species' thermal optimum is, but mathematically this is not completely correct.
  # If we let zo be the true optimum, z is the value we plug into the reproduction function so that argmax(b)==zo 
  # We apply the function zAdjust to all values of zo to calculate z
  z <- mapply(zAdjust, sig, zo, lam, 2^13)
  
  # To speed up computation time, we define full dispersal kernels now.
  # The dispersal kernel uses q, a transformation of gam
  q <- sapply(1:S, function(i) 1+1/gam[i]-sqrt(1+1/gam[i]^2))
  l <- (-L):(L)
  k <- t(sapply(1:S, function(i) doubGeom(l,q[i])))
  k[k<10^(-15)] <- 0.0
  # K is a list of separate L by 2L dispersal kernel matrices
  # K[[s]] is the dispersal kernel matrix of species s
  # K[[s]][i,] is a vector of probabilities for a propagule in patch i to spread to patch j-L/2 (this is extra long to account for a propagule spreading beyond the limits of the ecosystem)
  K <- rep(list(matrix(0,L,2*L)),S)
  for(i in 1:S){
    Ki<-matrix(0,2*L,4*L)
    for(j in 1:(2*L)){
      Ki[j,j:(j+2*L)]<-k[i,]
      Ki[j,3/2*L]<-sum(Ki[j,1:(3/2*L)])
      Ki[j,5/2*L+1]<-sum(Ki[j,(5/2*L+1):(4*L)])
    }
    K[[i]]<-Ki[(L/2+1):(3*L/2),(3/2*L):(5/2*L)]
  }
  # Tolerance and reproductive strength have a tradeoff
  # Birth rate is adjusted so each species has roughly equal birth rate when integrated over all T
  # ro is a constant that adjusts to this reproductive output
  ro <- sapply(sig,function(x) rAdjust(x,B,lam,1e-06,2^13))
  
  ##########################################################
  # Put the biological parameters together into a single list, P
  P=list(S=S,
         z=z,
         gam=gam,
         sig=sig,
         lam=lam,
         A=A,
         ro=ro,
         zo=zo,
         K=K,
         compType=compType)
  
  ##########################################################
  # Next, we define environmental parameters.
  # Discrete spatial domain: from 1 to L (integers)
  x <- seq(1,L)
  # temp2d is the current temperature over all x and microhabitats
  if(is.null(temp2d)){
    
    temp1dr <- seq(tempHigh,tempLow,length=L)
    if(tempRev){
      temp1dr<-rev(temp1dr)
    }
    tempG <- tempVarH(L,tempGH)
    temp1d <- temp1dr+tempG*tempGSD
    if(W>1){
      temp2dr <- matrix(temp1d,L,W)
      tempLH <- matrix(rnorm(L*W,0,tempLSD),L,W)
      tempLH <- t(matrix(sapply(1:L, function(x) sort(tempLH[x,]-mean(tempLH[x,]))),W,L))
      temp2d <- temp2dr+tempLH
    } else{
      temp2d<-temp1d
    }

  } else {
    if(!(nrow(temp2d)==L & ncol(tempsd)==W)){
      stop("temp2d does not match environment size!")
    }
    temp1d<-rowMeans(temp2d)
  } 
  # tempY is a vector of the temperature over time
  if(is.null(tempY)){
    tempY<-0:years
    for (i in 1:years){
      # A new epsi is calculated for each time step
      tempY[i+1] <- tempYAC*tempY[i]+rnorm(1,0,tempYSD)*sqrt(1-tempYAC^2)
    }
  } else{
    if(length(tempY)!=years+1){
      stop("tempY does not match years!")
    }
  }
    
  # Habitat quality could differ in space or with species, but we will keep it constant for now
  if(is.null(Q)){
    Qr1 <- matrix(QMean,L,W)
    QG <- tempVarH(L,QGH)
    Qr2 <- Qr1+QG*QGSD
    if(W>1){
      QLH <- matrix(rnorm(L*W,0,QLSD),L,W)
      QLH <- t(matrix(sapply(1:L, function(x) sort(QLH[x,]-mean(QLH[x,]))),W,L))
    }
    else{
      QLH<-0
    }
    Q <- Qr2+QLH
    
  } else {
    if(!(nrow(Q)==L & ncol(A)==W)){
      stop("Q does not match environment size!")
    }
  }  
  
  ##########################################################
  # Put the abiotic parameters together into a single list, X
  
  X=list(L=L,
         x=x,
         XW=XW,
         temp1d=temp1d,
         tau=tau,
         Q=Q,
         tempY=tempY,
         tempYAC=tempYAC,
         tempYSD=tempYSD,
         temp2d=temp2d,
         W=W,
         y=y)
  
  ##########################################################
  # Export it all as a list
  
  return(list(P=P,X=X))
}

tempVarH <-  function(L,H,cZero=T){
  # This function adds some heterogeneity to the temperature gradient with fractional Brownian noise
  # See Keitt (2000)
  # Spectral representation of neutral landscapes
  # Landscape Ecology 15
  
  # H:      Hurst exponent (should be between 0 and 1). It relates to the autocorrelation.
  ###        When H is near 1, this function has positive long-term positive autocorrelation and will look relatively smooth.
  ###        When H=0.5, this function is Brownian motion.
  ###        When H is near 0, then autocorrelation is negative and positive values will more often be followed by negative values (and vice versa).
  # L:     Length of the temperature gradient.
  # cZero: If T, it will center the whole output so the mean is 0.

  
  # random phases uniformly distributed on [0,2pi]
  phif <- runif(L)*2*pi
  
  # adjusted exponent for amplitudes
  betaH <- 1+2*H
  
  # uniformly distributed random numbers
  xf <- rnorm(L)
  # to form the amplitudes
  af <- 1/seq(1,L)^(betaH/2)*xf
  
  # complex coeffcients
  cf <- af*exp(1i*phif)
  
  #  real part of the inverse fourier transform
  tH <- Re(ifft(cf))
  
  # center it around zero?
  if(cZero){
    tH <- tH-mean(tH)
  }
  
  # multiply the output to increase the magnitude of the heterogeneity
  # add that to the the temperature gradient
  return(tH)
}


doubGeom<-function(x,q){
  # Probability mass function for "double geometric" distribution
  # x: distance from origin to landing spot
  # q: probability of remaining in a given patch (and not continuing to move); see supplemental
  
  return((q/(2-q)*(1-q)^abs(x)))
}

rAdjust<-function(sig,B,lam=-2.7,eps=1e-06,len=2^13){
  # This function creates a constant to adjust the reproduction rate so that the area under the curve is roughly equal for all species
  # sig: Thermal tolerance width of a species
  # B:   Desired total integrated area of positive growth
  # lam: Skewness in thermal tolerance
  # eps: Precision of estimate
  # len: Length of temperature vector. Higher values are more precise
  
  # Set up an extended version of a linear tempereature gradient
  temp <- seq(-100,100,length=len)
  
  # The actual optimal temperature is not important here, so we use the center of the temperature gradient
  z <- 20
  r <- exp(-(temp-z)^2/sig^2)*(1+erf(lam*(temp-z)/sig))-1
  
  bL <- -125; bH <- 125
  
  # Binary search for a value of ro such that exp(ro*r) integrates to B over all temperature values where exp(ro*r) is positive
  
  for(i in 1:500){
    bM <- (bL+bH)/2
    R <- exp(bM*r)
    G <- trapz(temp,(R-1)*(R>1))
    differ<- G-B
    if(abs(differ)<eps){
      break
    } else{
      if(differ>0){
        bH<-bM
      } else{
        bL<-bM
      }
    }
  }
  return(bM)
}

zAdjust<-function(sig,zo,lam=-2.7,L=2^13){
  # The reproduction function in Urban et al. 2012 is useful for creating the shape of the reproduction rate over temperature
  # However, the z_i "optimal temperature" doesn't end up where we might expect it to be
  # This function adjusts so that argmax_{temp1d}(R_i)=z_i
  # sig: The thermal tolerance width of a species
  # z:   Optimal temperature of species
  # lam: Skewness in thermal tolerance
  # len: Length of temperature vector. Higher values are more precise
  
  # Set up an extended version of a linear tempereature gradient
  temp <- seq(-100,100,length=L)
  
  # We need to calculate the difference between the expected optimal temperature and the actual optimal temperature
  # To do so, we begin with a baseline at zc=20
  zc <- 20
  
  # Calculate the baseline reproductive rate
  r<-exp(-(temp-zc)^2/sig^2)*(1+erf(lam*(temp-zc)/sig))-1
  
  # index for which temperature has the maximum reproductive output with the baseline
  iZ<-which.max(r)
  # index for baseline optimal temperature
  oZ<-which.min(abs(temp-zc))
  # index for desired optimal temperature
  tZ<-which.min(abs(temp-zo))
  
  # adjusted z to make optimal temperature in the right place
  z<-temp[tZ+oZ-iZ]
  
  return(z)
}

commSimulate <- function(n,P,X,years=100,extInit=F,extThresh=100){
  # This simulates a community, n, over yars
  # n:         Initial population sizes. SxLxW array of population nonnegative integers.
  # P:         List of biotic variables
  # X:         List of abiotic variables
  # years:     How many time steps to run the model.
  # extInit:   If T, the simulation stops running after extThresh time steps without any extinctions. When attempting to initialize a stable community, consider setting extInit to T.
  # extThresh: If extInit==T, then the simulation stops once extInit time steps have passed without any extinctions.
  y <- X$y
  
  # First, we set up a matrix to save the total population size over time
  N <- matrix(0,P$S,years+1)
  # Record the total initial population size of each species across the whole ecosystem
  N[,1] <- apply(n,1,sum)

  # Temperature changes over time, so we need to adjust this over the course of the model
  temp0 <- X$temp1d-X$tempY[y]
  temp2d0 <- X$temp2d-X$tempY[y]
  
  # For output, we want to keep track of the average temperature over time, and we do that with temps
  temps <- seq(0,years)
  temps[1] <- mean(temp0)
  
  # Keep track of tempY and tau outside of X
  tempY <- X$tempY
  tau <- X$tau

  # Used to track the last time an extinction occurred
  lastExtinction <- 0

  # Run the model for a number of time steps equal to 'years'
  for (i in 1:years){
    # Temperature changes before each time step
    X$temp1d=temp0+tau*i+tempY[i+y]
    X$temp2d=temp2d0+tau*i+tempY[i+y]

    # save the mean temperature
    temps[i+1]=mean(X$temp1d)
    
    # Run the time step function for population adjustment after the change in temperature
    if(sum(N[,i])>0){
      n <- timeStep(n,P,X)
    }
    # Record the population size
    N[,i+1]<-apply(n,1,sum)
    
    
    if(any(N[,i]>0 & N[,i+1]==0)){
      lastExtinction <- 0
    } else{
      lastExtinction <- lastExtinction+1
    }
    if(extInit==T & lastExtinction>extThresh){
      break()
    }
  }
  X$y<-X$y+i
  i<-i+1
  return(list(n=n,N=N[,1:i],temps=temps[1:i],X=X))
}

timeStep <- function(n,P,X){
  # Cycle through each step of the model.
  # Each time step could be one "year" or one "generation", but ultimately it runs through each part of the life cycle in an order determined by lcOrder.
  # The various management techniques can optionally be added to the model between any two of the required steps
  # reproduction -> dispersal -> density dependence
  n1 <- reproduce(n,X$L,X$W,P$S,P$z,P$sig,P$ro,P$lam,X$temp2d)
  
  dS<-which(rowSums(n1)>0)
  if(!isempty(dS)){
    dispn<-n1[dS,,]
    dispn2 <- disperse(dispn,X$L,X$W,P$S,P$K[dS])
    n2 <- n1
    n2[dS,,]<-dispn2
  } else {
    n2<-n1
  }
  n <- compete(n2,X$L,X$W,P$S,X$Q,P$A,P$compType,P$z,P$sig,P$ro,P$lam,X$temp2d)
  return(n)
}

bi <- function(z,sig,ro,lam,temp){
  # reproductive rate function
  op<-ro*(exp(-((temp-z)/sig)^2)*(1+erf(lam*(temp-z)/sig))-1)
  return(op)
}
  
reproduce <- function(n,L,W,S,z,sig,ro,lam,temp2d){
  # The number of offspring born for each species in each location is a Poisson random variable with mean r*n
  
  # The base reproductive rate is a skewed function, adjust such that min(r)=0 and max(r)=2
  # Each species will have a different reproductive rate depending on the temperature at that space.

  r <- sapply(1:S, function(i) bi(z[i],sig[i],ro[i],lam,temp2d))
  R <- exp(r)
  R <- aperm(array(R,c(L,W,S)),c(3,1,2))
  
  # Mean number of offspring
  rn <-c(R*n)
  
  # The number of offspring is a Poisson random variable with lambda=r*n
  nr<-array(sapply(rn, function(x) rpois(1,x)),c(S,L,W))

  return(nr)
}

disperse <- function(n,L,W,S,K){
  # Each individual spreads throughout the spatial landscape with a random double geometric dispersal kernel determined by the species' mean dispersal distance, gam[i].
  # For each species in each location, disperse!

  Si<-nrow(n)
  if(is.null(Si)){
    n<-array(n,c(S,L,W))
    Si<-S
  }
  
  n1<-apply(n,c(1,2),sum)
  n2 <- t(sapply(1:Si,function(j) disperseDoubGeom(n1[j,],L,K[[j]])))
  n3 <- c(sapply(1:Si, function(i) t(sapply(1:L,function(j) rebin(sample(1:W,n2[i,j],replace=T),W)))))
  n4 <- aperm(array(n3,c(L,W,Si)),c(3,1,2))

  return(n4)
}

disperseDoubGeom <- function(n,L,K){
  y <- which(n>0)
  n1<-sapply(y,function(x) rmultinom(1,n[x],K[x,]))
  if(length(y)>1){
    n2<-rowSums(n1)
  } else{
    n2<-n1
  }
  n3<-n2[2:(L+1)]
  return(n3)
}

compete <- function(n,L,W,S,Q,A,compType='lottery',z=NULL,sig=NULL,ro=NULL,lam=NULL,temp2d=NULL){
  # The density dependence in this model is roughly a Beverton-Holt model that includes both interspecific and intraspecific competition
  # Each individual has a random chance of survival based on a variety of conditions
  
  # Competition coefficients depend on interactions between each species and the temperature at the location at the time
  # These can be thought of as temperature-varying Lotka-Volterra competition coefficients
  # Probability of survival depends on competition coefficients, number of individuals of each different species at that location, and the quality of the habitat at that location
  
  if(compType=="temp"){
    r <- sapply(1:S, function(i) bi(z[i],sig[i],ro[i],lam,temp2d))
    R<- aperm(array(exp(r),c(L,W,S)),c(3,1,2))
  } else if(compType=="temp"){
    R<-1
  }
  Qrep<-array(rep(Q,each=S),c(S,L,W))
  QR<-1/(Qrep*R)
  nR<-R*n
  anR<-sapply(1:S, function(s) colSums(A[s,]*nR))
  anR<-aperm(array(anR,c(L,W,S)),c(3,1,2))
  
  p<-1/(1+QR*anR)

  nc <- (sapply(1:S,function(s) mapply(rbinom,1,c(n[s,,]),p[s,,])))
  nc2 <- array(t(nc),c(S,L,W))

  return(nc2)
}


unbin <- function(v){
  # Convert vector of population sizes over x into a vector of the location for each individual
  L<-sum(v)
  ub<-matrix(0,L)
  j<-0
  for(i in which(v>0)){
    ub[(j+1):(j+v[i])]<-i
    j<-j+v[i]
  }
  return(c(ub))
}

rebin <- function(ub,L){
  # Converts a vector of individual locations into a vector of population sizes over x
  v<-1:L
  for(i in 1:L){
    v[i]<-length(which(ub==i))
  }
  return(v)
}


diversityIndex <- function(n,type="alpha",index="invSimp"){
  # Calculates the diversity index of a matrix of population sizes
  if(type=='alpha'){
    p <- t(n)/colSums(n)
    p2 <- rowSums(p^2)
    p2[is.nan(p2)] <- 1
    if(index=="invSimp"){
      D <- mean(1/p2)
    } else if(index=="giniSimp"){
      D <- 1- mean(p2)
    }
  } else if (type=='gamma'){
    p <- rowSums(n)/sum(n)
    D <- 1/sum(p^2)
    if(index=="invSimp"){
      D <- 1/sum(p^2)
    } else if(index=="giniSimp"){
      D <- 1- sum(p^2)
    }  
  }
  return(D)
}

commTrim <- function(n,P,X){
  # Remove extinct species from n and P
  nFlat <- t(sapply(1:P$S, function(s) rowSums(n[s,,,drop=FALSE])))
  extant <- which(rowSums(nFlat)>0)
  P$S <- length(extant)
  P$z <- P$z[extant,drop=FALSE]
  P$gam <- P$gam[extant,drop=FALSE]
  P$sig <- P$sig[extant,drop=FALSE]
  P$A <- P$A[extant,extant,drop=FALSE]
  P$ro <- P$ro[extant,drop=FALSE]
  P$zo <- P$zo[extant,drop=FALSE]
  P$K <- P$K[extant]
  n <- n[extant,,,drop=FALSE]
  ct <- list(n=n,P=P)
  return(ct)
}


###### Useful functions for visualization
topSpecies <- function(n,L,W){
  tSpec <- matrix(NA,L,W)
  for(i in 1:L){
    for(j in 1:W){
      nij <- n[,i,j]
      if(sum(nij)>0){
        nij<-nij+runif(32,0,0.001)
        tSpec[i,j]<-which.max(nij)
      }
    }
  }
  return(tSpec)
}

vComTop<- function(tSpec,S){
  image(tSpec,col=rainbow(S))
}

vComSide<-function(n,S){
  matplot(sapply(1:S, function(s) rowSums(n[s,,])),type="l",lwd=2,lty=1,ylim=c(0,90),col=rainbow(S))
}


#### Example simulation
id<-1
set.seed(id)

S <- 1
L <- 512
W <- 1
tau <- 0.04
tempYSD <- 0.2
tempLSD <- 1
QMean=4

iYears <- 200
ccYears <- 100

cSetup<-commSetup(S=S,L=L,W=W,compType="temp",tau=0,years=iYears+ccYears,tempYSD=tempYSD,tempLSD=tempLSD,QMean=QMean)
P0<-cSetup$P
X0<-cSetup$X

n0 <- array(4,c(P0$S,X0$L,X0$W))


cSim1<-commSimulate(n0,P0,X0,years=iYears)

n0f <- cSim1$n

ct1 <- commTrim(n0f,P0,X0)
n1 <- ct1$n
P1 <- ct1$P

X1 <- X0
X1$tau <- tau

cSim2<-commSimulate(n1,P1,X1,years=ccYears)

n1f <- cSim2$n

ct2 <- commTrim(n1f,P1,X1)
n2 <- ct2$n
P2 <- ct2$P
