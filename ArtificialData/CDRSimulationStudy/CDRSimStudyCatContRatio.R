require(MixSim)
require(clustrd)

# Simulation study: Effect of kappa in CDR for varying categorical:continuous attributes ratio

# Define grid of kappa values
kappas <- seq(0, 1, length.out=100)

# Empty matrices to store simulation results
kappasmat <- matrix(NA, nrow=length(kappas), ncol=9)
kappasmat2 <- matrix(NA, nrow=length(kappas), ncol=36)

niter <- 20

# Simulation for 1 artificial data set
set.seed(1234)
# Construction of artificial data set with:
  # Average overlap: 0.05
  # Maximum Overlap: 0.0525
  # Number of clusters: 3
  # Number of observations: 480
  # Number of variables: 10
mixsim1 <- MixSim(BarOmega = 0.05, MaxOmega = 0.0525, PiLow=1.0, K = 3, p = 10, resN = 10000)
mixdt1 <- simdataset(n = 480, Pi = mixsim1$Pi, Mu = mixsim1$Mu, S = mixsim1$S)
for (i in 1:9){
  # Discretise first i attributes
  for (k in 1:i){
    mixdt1$X[,k] <- discretise(mixdt1$X[,k], nlevels=3)
  }
  mixdt1df <- as.data.frame(mixdt1$X)
  for (k in 1:i){
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
    kappasmat[index,i] <- mean(arivec)
    cat('Iteration',index,'of',i,'complete \n')
    index <- index+1
    save(kappasmat, file='kappasmat.RData')
  }
}

# Simulation for 4 different artificial data sets
for (l in 1:4){
  set.seed(1234+l)
  # Construction of artificial data set
  mixsim1 <- MixSim(BarOmega = 0.05, MaxOmega = 0.0525, PiLow=1.0, K = 3, p = 10, resN = 10000)
  mixdt1 <- simdataset(n = 480, Pi = mixsim1$Pi, Mu = mixsim1$Mu, S = mixsim1$S)
  # Discretise first i attributes
  for (i in 1:9){
    for (k in 1:i){
      mixdt1$X[,k] <- discretise(mixdt1$X[,k], nlevels=3)
    }
    mixdt1df <- as.data.frame(mixdt1$X)
    for (k in 1:i){
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
      kappasmat2[index,(l-1)*9+i] <- mean(arivec)
      cat('Iteration',index,'of',i,'for data set',l+1,'complete \n')
      index <- index+1
      save(kappasmat2, file='kappasmat2.RData')
    }
  }
}