rm(list=ls())
require(pracma)
require(parallel)

commSetup <- function(S=32, L=512, W=8,
                      zo=NULL, gam=NULL, sig=NULL, A=NULL,
                      gamMean=2.5, gamSD=2.5,
                      sigMean=5, sigSD=5,
                      lam=-2.7, B=10,
                      XW=seq(129,384),
                      temp2d=NULL,
                      tempLow=9.78, tempHigh=30.22,
                      tempRev=F,
                      tempGH=1, tempGSD=0, tempLSD=1,
                      years=1000,
                      y=0,
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
  z <- mapply(zAdjust, sig, zo, lam, 2^11)
  
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
  ro <- sapply(sig,function(x) rAdjust(x,B,lam,0.00001,2^11))
  
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
         K=K)
  
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
    
    temp2dr <- matrix(temp1d,L,W)
    tempLH <- matrix(rnorm(L*W,0,tempLSD),L,W)
    tempLH <- t(matrix(sapply(1:L, function(x) sort(tempLH[x,]-mean(tempLH[x,]))),W,L))
    temp2d <- temp2dr+tempLH
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
    QLH <- matrix(rnorm(L*W,0,QLSD),L,W)
    QLH <- t(matrix(sapply(1:L, function(x) sort(QLH[x,]-mean(QLH[x,]))),W,L))
    
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


doubGeom<-function(k,p){
  # Probability mass function for "double geometric" distribution
  return((p/(2-p)*(1-p)^abs(k)))
}

rAdjust<-function(sig,B,lam,eps,L){
  # This function creates a constant to adjust the reproduction rate so that the area under the curve is roughly equal for all species
  
  # Set up an extended version of a linear tempereature gradient
  temp <- seq(-100,100,length=L)
  
  # The actual optimal temperature is not important here, so we use the center of the temperature gradient
  z <- 20
  r <- exp(-(temp-z)^2/sig^2)*(1+erf(lam*(temp-z)/sig))-1
  
  bL <- -125; bH <- 125
  
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

zAdjust<-function(sig,z,lam,L){
  # The reproduction function in Urban et al. 2012 is useful for creating the shape of the reproduction rate over temperature
  # However, the z_i "optimal temperature" doesn't end up where we might expect it to be
  # This function adjusts so that argmax_{temp1d}(R_i)=z_i
  
  # Set up an extended version of a linear tempereature gradient
  temp1d <- seq(-100,100,length=L)
  
  # We need to calculate the difference between the expected optimal temperature and the actual optimal temperature
  # To do so, we begin with a baseline at zo=20
  zo <- 20
  
  # Calculate the baseline reproductive rate
  r<-exp(-(temp1d-zo)^2/sig^2)*(1+erf(lam*(temp1d-zo)/sig))-1
  
  # index for which temperature has the maximum reproductive output with the baseline
  iZ<-which.max(r)
  # index for baseline optimal temperature
  oZ<-which.min(abs(temp1d-zo))
  # index for desired optimal temperature
  tZ<-which.min(abs(temp1d-z))
  
  # adjusted zi to make optimal temperature in the right place
  zf<-temp1d[tZ+oZ-iZ]
  
  return(zf)
}

commSimulate <- function(n,P,X,years=100){
  # This simulates a community, n, over yars
  
  # First, we set up a matrix to save the total population size over time
  N <- matrix(0,P$S,years+1)
  # Record the initial population size
  N[,1] <- apply(n,1,sum)

  # Temperature changes over time, so we need to adjust this over the course of the model
  temp0 <- X$temp1d
  temp2d0 <- X$temp2d
  
  # For output, we want to keep track of the average temperature over time, and we do that with temps
  temps <- seq(0,years)
  temps[1] <- mean(temp0)
  
  # Keep track of tempY and tau outside of X
  tempY <- X$tempY
  tau <- X$tau
  y <- X$y
  
  lastExtinction <- 0

  # Run the model for years
  
  for (i in 1:years){
    # Temperature changes before each time step
    X$temp1d=temp0+tau*i+tempY[i+y+1]
    X$temp2d=temp2d0+tau*i+tempY[i+y+1]

    # save the mean temperature
    temps[i+1]=mean(X$temp1d)
    
    # Run the time step function for population adjustment after the change in temperature
    n <- timeStep(n,P,X)
    
    # Record the population size
    N[,i+1]<-apply(n,1,sum)
    
    
    if(any(N[,i]>0 & N[,i+1]==0)){
      lastExtinction <- 0
    } else{
      lastExtinction <- lastExtinction+1
    }
    #if(lastExtinction>100){
    #  break()
    #}
  }
  X$y<-i+1
  return(list(n=n,N=N[,1:i],temps=temps[1:i],X=X))
}

tempStoch <- function(years,X){
  epsi<-0:years
  for (i in 1:years){
    # A new epsi is calculated for each time step
    epsi[i+1] <- X$tempYAC*epsi[i]+rnorm(1,0,X$tempYSD)*sqrt(1-X$tempYAC^2)
  }
  return(epsi)
}

timeStep <- function(n,P,X){
  # Cycle through each step of the model.
  # Each time step could be one "year" or one "generation", but ultimately it runs through each part of the life cycle in an order determined by lcOrder.
  # The various management techniques can optionally be added to the model between any two of the required steps
  # reproduction -> dispersal -> density dependence
  
  n1 <- reproduce(n,P,X)

  n2 <- disperse(n1,P,X)
  
  n <- compete(n2,P,X)
  
  return(n)
}

bi <- function(z,sig,ro,lam,temp1d){
  # reproductive rate function
  op<-ro*(exp(-((temp1d-z)/sig)^2)*(1+erf(lam*(temp1d-z)/sig))-1)
  return(op)
}
  
reproduce <- function(n,P,X){
  # The number of offspring born for each species in each location is a Poisson random variable with mean r*n
  
  # The base reproductive rate is a skewed function, adjust such that min(r)=0 and max(r)=2
  # Each species will have a different reproductive rate depending on the temperature at that space.
  r <- sapply(1:P$S, function(i) bi(P$z[i],P$sig[i],P$ro[i],P$lam,X$temp2d))
  R <- exp(r)
  R <- aperm(array(R,c(X$L,X$W,P$S)),c(3,1,2))
  
  # Mean number of offspring
  rn <-c(R*n)
  
  # The number of offspring is a Poisson random variable with lambda=r*n
  nr<-array(sapply(rn, function(x) rpois(1,x)),c(P$S,X$L,X$W))

  return(nr)
}

disperse <- function(n,P,X){
  # Each individual spreads throughout the spatial landscape with a random double geometric dispersal kernel determined by the species' mean dispersal distance, gam[i].
  # For each species in each location, disperse!
  dSpecies<-which(rowSums(n)>0)
  nd<-n
  
  if(length(dSpecies)){
    ndd<-apply(nd,c(1,2),sum)
    
    Si<-length(dSpecies)
    ndd <- t(sapply(dSpecies,function(j) disperseDoubGeom(ndd[j,],X$L,P$K[[j]])))
    nddd <- c(sapply(1:Si, function(i) t(sapply(1:X$L,function(j) rebin(sample(1:8,ndd[i,j],replace=T),8)))))
    ndd <- aperm(array(nddd,c(X$L,8,Si)),c(3,1,2))
    nd[dSpecies,,]<-ndd
  }
  return(nd)
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

compete <- function(n,P,X){
  # The density dependence in this model is roughly a Beverton-Holt model that includes both interspecific and intraspecific competition
  # Each individual has a random chance of survival based on a variety of conditions
  
  # Competition coefficients depend on interactions between each species and the temperature at the location at the time
  # These can be thought of as temperature-varying Lotka-Volterra competition coefficients
  # Probability of survival depends on competition coefficients, number of individuals of each different species at that location, and the quality of the habitat at that location
  
  r <- sapply(1:P$S, function(i) bi(P$z[i],P$sig[i],P$ro[i],P$lam,X$temp2d))
  R <- exp(r)

  an <- sapply(0:(X$W-1), function(j) sapply(1:X$L, function(i) sum(R[i+j*X$L,]*n[,i,j+1])))

  p <- sapply(1:P$S, function(s) 1/(1+c(an)/(c(R[,s])*X$Q)))

  
  # For each species in each location, compete!
  # Number of individuals remaining after competition is a biomial random variable with
  # n is the number of individuals of species i in location x
  # p is probability of survival
  nc <- (sapply(1:P$S,function(s) mapply(rbinom,1,c(n[s,,]),p[,s])))
  nc2 <- array(t(nc),c(P$S,X$L,8))

  return(nc2)
}


unbin <- function(v){
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
  v<-1:L
  for(i in 1:L){
    v[i]<-length(which(ub==i))
  }
  return(v)
}


invSimp <- function(n,type='alpha'){
  if(type=='alpha'){
    p <- t(n)/colSums(n)
    p2 <- rowSums(p^2)
    p2[is.nan(p2)] <- 1
    D <- mean(1/p2) 
  } else if (type=='gamma'){
    p <- rowSums(n)/sum(n)
    D <- 1/sum(p^2)
  }
  return(D)
}

commData <- function(XW,sim,P){
  nFlat<-t(sapply(1:P$S, function(s) rowSums(sim$n[s,,])))
  extInt<-which( rowSums(nFlat[,XW])>0)
  nFlatInt<-nFlat[extInt,XW]
  specRich<-c(mean(colSums(nFlatInt>0)),length(extInt))
  simp<-c(invSimp(nFlatInt,'alpha'),invSimp(nFlatInt,'gamma'))
  disp<-c(mean(P$gam[extInt]),sd(P$gam[extInt]))
  toler<-c(mean(P$sig[extInt]),sd(P$sig[extInt]))
  popSize<-c(mean(rowSums(nFlatInt)),sd(rowSums(nFlatInt)))
  fRange<-sapply(1:length(extInt), function(s) max(which(nFlatInt[s,]>0))-min(which(nFlatInt[s,]>0)))
  qRange<-sapply(1:length(extInt), function(s) quantile(unbin(nFlatInt[s,]),0.975)-quantile(unbin(nFlatInt[s,]),0.025))
  fRanges<-c(mean(fRange),sd(fRange))
  qRanges<-c(mean(qRange),sd(qRange))
  outPut<-c(specRich,simp,disp,toler,popSize,fRanges,qRanges)
  return(outPut)
}

outputData <- function(id,action,tempYSD,lHSD,ccYears,P2,X2,sim1,sim2,nF2,nF3){
  opd<-data.frame(bID=0)
  
  opd$ID<-id
  opd$action<-action
  
  opd$ES<-tempYSD
  opd$het<-lHSD
  opd$S<-P2$S
  
  opd$gam <- P2$gam[i]
  opd$sig <- P2$sig[i]
  opd$zo <- P2$zo[i]
  
  opd$N0 <- sim2$N[i,1]
  opd$Nf <- sim2$N[i,ccYears+1]
  zE <- P2$zo[P2$zo>P2$zo[i]]
  zP <- P2$zo[P2$zo<P2$zo[i]]
  opd$zDifE <- zE[which.min(abs(P2$zo[i]-zE))]-P2$zo[i]
  opd$zDifP <- P2$zo[i]-zP[which.min(abs(P2$zo[i]-zP))]
  opd$sdN <- sd(sim2$N[i,])
  opd$minN <- min(sim2$N[i,])
  opd$maxN <- max(sim2$N[i,])
  opd$meanN <- mean(sim2$N[i,])
  opd$ttExtinct <- max(which(sim$N[i,]>0))

  opd$extinct <- opd$Nf==0
  if(opd$ttExtinct==ccYears+1){opd$ttExtinct <- NaN}
  
  XW <- seq(1*X$L/4+1,3*X$L/4)
  
  nFw1 <- nF2[,XW]; nFw2 <- nF3[,XW]
  Nw1 <- rowSums(nFw1);  Nw2 <- rowSums(nFw2); 
  wS1 <- which(Nw1>0);   wS2 <- which(Nw2>0)

  opd$Ntot1 <- sum(Nw1); opd$Ntot2 <- sum(Nw2)
  opd$gR1 <- sum(Nw1>0); opd$gR2 <- sum(Nw2>0)
  opd$aR1 <- mean(colSums(nFw1>0)); opd$aR2 <- mean(colSums(nFw2>0))
  opd$gS1 <- invSimp(nFw1,'gamma'); opd$gS2 <- invSimp(nFw2,'gamma')
  opd$aS1 <- invSimp(nFw1,'alpha'); opd$aS2 <- invSimp(nFw2,'alpha')
    
  opd$gamMean1 <- mean(P$gam[wS1]); opd$gamMean2 <- mean(P$gam[wS2])
  opd$gamSD1 <- sd(P$gam[wS1]); opd$gamSD2 <- sd(P$gam[wS2])
  opd$sigMean1 <- mean(P$sig[wS1]); opd$sigMean2 <- mean(P$sig[wS2])
  opd$sigSD1 <- sd(P$sig[wS1]); opd$sigSD2 <- sd(P$sig[wS2])
  
  opd$popSizeMean1 <- mean(Nw1); opd$popSizeMean2 <- mean(Nw2)
  opd$popSizeSD1 <- sd(Nw2); opd$popSizeSD2 <- sd(NW2)
  
  fRange1<-sapply(1:S, function(s) max(which(nFw1[s,]>0))-min(which(nFw1[s,]>0)))
  qRange1<-sapply(1:S, function(s) quantile(unbin(nFw1[s,]),0.975)-quantile(unbin(nFw1[s,]),0.025))
  fRange2<-sapply(1:S, function(s) max(which(nFw2[s,]>0))-min(which(nFw2[s,]>0)))
  qRange2<-sapply(1:S, function(s) quantile(unbin(nFw2[s,]),0.975)-quantile(unbin(nFw2[s,]),0.025))
  
  opd$fRangeMean1 <- mean(fRange1); opd$fRangeMean2 <- mean(fRange2);
  opd$fRangeSD1 <- sd(fRange1); opd$fRangeSD2 <- sd(fRange2);
  opd$qRangeMean1 <- mean(qRange1); opd$qRangeMean2 <- mean(qRange2);
  opd$qRangeSD1 <- sd(qRange1); opd$qRangeSD2 <- sd(qRange2);
  
  opdf <- data.matrix(opd)
  colnames(opdf) <- colnames(opd)
  return(opd)
}

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

vComTop<- function(tSpec){
  image(tSpec,col=rainbow(32))
}

vComSide<-function(n,S){
  matplot(sapply(1:S, function(s) rowSums(n[s,,])),type="l",lwd=2,lty=1,ylim=c(0,90))
}

commTrim <- function(n,P){
  nFlat <- t(sapply(1:P$S, function(s) rowSums(n[s,,])))
  extant <- which(rowSums(nFlat)>0)
  P$S <- length(extant)
  P$z <- P$z[extant]
  P$gam <- P$gam[extant]
  P$sig <- P$sig[extant]
  P$A <- P$A[extant,extant]
  P$ro <- P$ro[extant]
  P$zo <- P$zo[extant]
  P$K <- P$K[extant]
  return(P)
}

#### Example simulation
  id<-1
  set.seed(id)
  
  S <- 64
  L <- 512
  W <- 8
  tau <- 0.04
  tempYSD <- 0.2
  tempLSD <- 1
  
  iYears <- 200
  ccYears <- 100
  
  cSetup<-commSetup(S=S,L=L,W=W,tau=0,years=iYears+ccYears,tempYSD=tempYSD,tempLSD=tempLSD)
  P0<-cSetup$P
  X0<-cSetup$X
  
  n0 <- array(4,c(P0$S,X0$L,X0$W))
  
  
  cSim1<-commSimulate(n0,P0,X0,years=iYears)
  
  n1 <- cSim1$n
  
  P1 <- commTrim(n1,P0)
  X1 <- X0
  X1$tau <- tau

  
  
  runif1<-runif(2)
  tempYSD<-runif1[1]
  lHSD<-runif1[2]*2
  testSetup <- hetSetup(S,tempYSD=tempYSD,lHSD=lHSD)
  X1<-testSetup$X; P1<-testSetup$P
  X1$tau<-0
  
  tempt<-tempStoch(iYears+ccYears+2,X1)
  
  
  # Initialize the initial species distributions
  # Popuations start everywhere
  n1 <- array(4,c(P1$S,X1$L,X1$W))
  
  sim1<-hetSimulate(n1,P1,X1,ccYears=iYears,ccTemps=tempt[1:(iYears+1)])
  
  XW<-seq(129,384)
  
  cData1<-commData(XW,sim1,P1)
  
  nFlat1<-t(sapply(1:P1$S, function(s) rowSums(sim1$n[s,,])))
  extant1<-which(rowSums(nFlat1)>0)
  
  P2<-P1
  P2$S<-length(extant1)
  P2$z<-P1$z[extant1]
  P2$gam<-P1$gam[extant1]
  P2$sig<-P1$sig[extant1]
  P2$A<-P1$A[extant1,extant1]
  P2$ro<-P1$ro[extant1]
  P2$zo<-P2$zo[extant1]
  P2$Rmax<-P2$Rmax[extant1]
  P2$k<-P2$k[extant1,]
  P2$K<-P2$K[extant1]
  
  X2<-X1
  X2$tau<-dTemp # Yearly temperature increase is now 0.04 degrees per year (on average)
  X2$temp1d<-X2$temp1d+sim1$temps[sim1$tFinal]-mean(X2$temp1d)
  
  n2<-sim1$n[extant1,,]
  
  sim2<-hetSimulate(n2,P2,X2,ccYears=ccYears+1,ccTemps=tempt[(iYears+1):(iYears+ccYears+2)])
  
  cData2<-commData(XW,sim2,P2)
  
  pt2<-proc.time()
