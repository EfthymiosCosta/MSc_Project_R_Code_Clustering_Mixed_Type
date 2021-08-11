require(MixSim)
require(clustrd)

# Simulation study: Effect of kappa in CDR for varying number of categorical levels

# Define grid of kappa values
kappas <- seq(0, 1, length.out=100)

# Empty matrices to store simulation results
kappasmatlevels <- matrix(NA, nrow=length(kappas), ncol=5)
kappasmatlevels2 <- matrix(NA, nrow=length(kappas), ncol=20)

niter <- 20

# Simulation for 1 artificial data set
set.seed(1234)
# Construction of artificial data set with:
  # Average overlap: 0.05
  # Maximum Overlap: 0.0525
  # Number of clusters: 3
  # Number of observations: 480
  # Number of variables: 10
mixsim1 <- MixSim(BarOmega = 0.05, MaxOmega = 0.0525, PiLow=1.0, K = 3, p = 10, resN = 100000)
mixdt1 <- simdataset(n = 480, Pi = mixsim1$Pi, Mu = mixsim1$Mu, S = mixsim1$S)
for (i in 1:5){
  # Discretise first 2 attributes using i+2 categorical levels
  for (k in 1:2){
    mixdt1$X[,k] <- discretise(mixdt1$X[,k], nlevels=i+2)
  }
  mixdt1df <- as.data.frame(mixdt1$X)
  for (k in 1:2){
    mixdt1df[,k] <- as.factor(mixdt1df[,k])
  }
  index <- 1
  for (kappa in kappas){
    # Initialise vector of ARIs
    arivec <- c()
    for (j in 1:niter){
      # CDR with specified kappa
      outmixcdr = cluspcamix(mixdt1df, nclus = 3, ndim = 2, alpha = kappa)
      # Update vector of ARIs
      arivec <- c(arivec, adjustedRandIndex(outmixcdr$cluster, mixdt1$id))
    }
    # Save mean ARI for given kappa
    kappasmatlevels[index,i] <- mean(arivec)
    save(kappasmatlevels, file='kappasmatlevels.RData')
    cat('Iteration',index,'of',i,'complete \n')
    index <- index+1
  }
}

# Simulation for 4 different artificial data sets
for (l in 1:4){
  set.seed(1234+l)
  # Construction of artificial data set
  mixsim1 <- MixSim(BarOmega = 0.05, MaxOmega = 0.0525, PiLow=1.0, K = 3, p = 10, resN = 100000)
  mixdt1 <- simdataset(n = 480, Pi = mixsim1$Pi, Mu = mixsim1$Mu, S = mixsim1$S)
  # Discretise first 2 attributes using i+2 categorical levels
  for (i in 1:5){
    for (k in 1:2){
      mixdt1$X[,k] <- discretise(mixdt1$X[,k], nlevels=i+2)
    }
    mixdt1df <- as.data.frame(mixdt1$X)
    for (k in 1:2){
      mixdt1df[,k] <- as.factor(mixdt1df[,k])
    }
    index <- 1
    for (kappa in kappas){
      arivec <- c()
      for (j in 1:niter){
        # CDR with specified kappa
        outmixcdr = cluspcamix(mixdt1df, nclus = 3, ndim = 2, alpha = kappa)
        # Update vector of ARIs
        arivec <- c(arivec, adjustedRandIndex(outmixcdr$cluster, mixdt1$id))
      }
      # Save mean ARI for given kappa
      kappasmatlevels2[index,(l-1)*5+i] <- mean(arivec)
      save(kappasmatlevels2, file='kappasmatlevels2.RData')
      cat('Iteration',index,'of',i,'for data set',l+1,'complete \n')
      index <- index+1
    }
  }
}