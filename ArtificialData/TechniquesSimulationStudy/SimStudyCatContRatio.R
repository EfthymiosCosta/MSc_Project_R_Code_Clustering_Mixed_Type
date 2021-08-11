require(MixSim)
require(cluster)
require(clustrd)
require(clustMixType)
require(mclust)
require(kmed)
require(FactoMineR)
require(fpc)
require(mclust)
require(kamila)
require(homals)

# Simulation study on effect of categorical:continuous attributes ratio

niter <- 100

# Empty matrices to store simulation results
mixsimresmat <- matrix(NA, nrow=100, ncol=72)
mixsimresmat2 <- matrix(NA, nrow=100, ncol=288)

# Empty matrices to store artificial data sets constructed
fulldtmat <- matrix(NA, nrow=nrow(mixdt1$X), ncol=90)
fulldtmat2 <- matrix(NA, nrow=480, ncol=360)

# Simulation for 1 artificial data set
set.seed(1234)
# Construction of artificial data set with:
  # Average overlap: 0.05
  # Maximum Overlap: 0.0525
  # Number of clusters: 3
  # Number of observations: 480
  # Number of variables: 10
mixsim1 <- MixSim(BarOmega = 0.05, MaxOmega = 0.0525, PiLow=1.0, K = 3, p = 10, resN=10000)
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
  # Store final data set
  fulldtmat[, ((i-1)*10+1):(i*10)] <- as.matrix(mixdt1df)
  save(fulldtmat, file='fulldtmat.RData')
  # Calculate Gower's dissimilarities and Ahmad distances
  gower_dist <- daisy(mixdt1df, metric = "gower")
  if (i>1){
    mix_sim <- distmix(mixdt1df, method = "ahmad", idcat = c(1:i), idnum = c((i+1):10))
  }
  # Specify continuous & categorical attributes
  conDf <- data.frame(scale(mixdt1df[,(i+1):10]))
  catDf <- dummyCodeFactorDf(data.frame(mixdt1df[,1:i]))
  # Number of dimensions to be retained using GCV
  howmany <- estim_ncp(data.matrix(mixdt1df))$ncp
  if (howmany<=1){
    howmany <- 2
  }
  outpcamix <- FAMD(mixdt1df, ncp = howmany, graph=FALSE)
  nlpca <- homals(mixdt1df, active=TRUE, verbose=1, rank=1,
                  level=c(rep('nominal', i), rep('numerical', (10-i))), ndim = howmany)
  for (j in 1:niter){
    # PAM with Gower's
    pam_fit <- pam(gower_dist, diss = TRUE, k = 3, do.swap = FALSE, cluster.only = TRUE,
                   medoids=sample(1:nrow(mixdt1df), 3))
    mixsimresmat[j,((i-1)*8+1)] <- adjustedRandIndex(pam_fit, mixdt1$id)
    cat('PAM done for iteration',j,'of',i,'\n')
    # K-prototypes
    outk <- kproto(mixdt1df, 3)
    mixsimresmat[j,((i-1)*8+2)] <- adjustedRandIndex(outk$cluster, mixdt1$id)
    cat('K-prototypes done for iteration',j,'of',i,'\n')
    # Mixed K-means
    if (i>1){
      kmedres <- fastkmed(mix_sim, 3, iterate = 100, init = sample(1:nrow(mixdt1df), 3))
      mixsimresmat[j,((i-1)*8+3)] <- adjustedRandIndex(kmedres$cluster, mixdt1$id)
      cat('Mixed K-means done for iteration',j,'of',i,'\n')
    }
    # Modha-Spangler
    if (i<9){
      msRes <- gmsClust(conDf, catDf, nclust = 3, searchDensity = 5)
      mixsimresmat[j,((i-1)*8+4)] <- adjustedRandIndex(msRes$results$cluster, mixdt1$id)
      cat('Modha-Spangler done for iteration',j,'of',i,'\n')
    }
    # FAMD + K-means
    outkm <- kmeans(outpcamix$ind$coord, 3, iter.max=100000, algorithm = "MacQueen")
    mixsimresmat[j,((i-1)*8+5)] <- adjustedRandIndex(outkm$cluster, mixdt1$id)
    cat('FAMD done for iteration',j,'of',i,'\n')
    # NLPCA + K-means
    nlpcakmeans <- kmeans(nlpca$objscores, 3, iter.max=100000, algorithm = "MacQueen")
    mixsimresmat[j,((i-1)*8+6)] <- adjustedRandIndex(nlpcakmeans$cluster, mixdt1$id)
    cat('NLPCA done for iteration',j,'of',i,'\n')
    # Mixed RKM
    outmix = cluspcamix(mixdt1df, nclus = 3, ndim = 2)
    mixsimresmat[j,((i-1)*8+7)] <- adjustedRandIndex(outmix$cluster, mixdt1$id)
    cat('Mixed RKM done for iteration',j,'of',i,'\n')
    # Mixed FKM
    outfmix = cluspcamix(mixdt1df, nclus = 3, ndim = 2, method = 'mixedFKM')
    mixsimresmat[j,((i-1)*8+8)] <- adjustedRandIndex(outfmix$cluster, mixdt1$id)
    cat('Mixed FKM done for iteration',j,'of',i,'\n')
    cat('Iteration',j,'of',i,'complete \n')
  }
  save(mixsimresmat, file='mixsimresmat.RData')
}

# Simulation for 4 different artificial data sets
for (l in 1:4){
  set.seed(1234+l)
  # Construct artificial data set
  mixsim1 <- MixSim(BarOmega = 0.05, MaxOmega = 0.0525, PiLow=1.0, K = 3, p = 10, resN=10000)
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
    # Store final data set
    fulldtmat2[, (90*(l-1)+(i-1)*10+1):(90*(l-1)+(i*10))] <- as.matrix(mixdt1df)
    save(fulldtmat2, file='fulldtmat2.RData')
    # Calculate Gower's dissimilarities
    gower_dist <- daisy(mixdt1df, metric = "gower")
    if (i>1){
      mix_sim <- distmix(mixdt1df, method = "ahmad", idcat = c(1:i), idnum = c((i+1):10))
    }
    # Specify continuous & categorical attributes
    conDf <- data.frame(scale(mixdt1df[,(i+1):10]))
    catDf <- dummyCodeFactorDf(data.frame(mixdt1df[,1:i]))
    # Number of dimensions to be retained using GCV
    howmany <- estim_ncp(data.matrix(mixdt1df))$ncp
    if (howmany<=1){
      howmany <- 2
    }
    outpcamix <- FAMD(mixdt1df, ncp = howmany, graph=FALSE)
    nlpca <- homals(mixdt1df, active=TRUE, verbose=1, rank=1,
                    level=c(rep('nominal', i), rep('numerical', (10-i))), ndim = howmany)
    for (j in 1:niter){
      # PAM with Gower's
      pam_fit <- pam(gower_dist, diss = TRUE, k = 3, do.swap = FALSE, cluster.only = TRUE,
                     medoids=sample(1:nrow(mixdt1df), 3))
      mixsimresmat2[j,(72*(l-1)+(i-1)*8+1)] <- adjustedRandIndex(pam_fit, mixdt1$id)
      cat('PAM done for iteration',j,'of',i,'\n')
      # K-prototypes
      outk <- kproto(mixdt1df, 3)
      mixsimresmat2[j,(72*(l-1)+(i-1)*8+2)] <- adjustedRandIndex(outk$cluster, mixdt1$id)
      cat('K-prototypes done for iteration',j,'of',i,'\n')
      # Mixed K-means
      if (i>1){
        kmedres <- fastkmed(mix_sim, 3, iterate = 100, init = sample(1:nrow(mixdt1df), 3))
        mixsimresmat2[j,(72*(l-1)+(i-1)*8+3)] <- adjustedRandIndex(kmedres$cluster, mixdt1$id)
        cat('Mixed K-means done for iteration',j,'of',i,'\n')
      }
      # Modha-Spangler
      if (i<9){
        msRes <- gmsClust(conDf, catDf, nclust = 3, searchDensity = 5)
        mixsimresmat2[j,(72*(l-1)+(i-1)*8+4)] <- adjustedRandIndex(msRes$results$cluster, mixdt1$id)
        cat('Modha-Spangler done for iteration',j,'of',i,'\n')
      }
      # FAMD + K-means
      outkm <- kmeans(outpcamix$ind$coord, 3, iter.max=100000, algorithm = "MacQueen")
      mixsimresmat2[j,(72*(l-1)+(i-1)*8+5)] <- adjustedRandIndex(outkm$cluster, mixdt1$id)
      cat('FAMD done for iteration',j,'of',i,'\n')
      # NLPCA + K-means
      nlpcakmeans <- kmeans(nlpca$objscores, 3, iter.max=100000, algorithm = "MacQueen")
      mixsimresmat2[j,(72*(l-1)+(i-1)*8+6)] <- adjustedRandIndex(nlpcakmeans$cluster, mixdt1$id)
      cat('NLPCA done for iteration',j,'of',i,'\n')
      # Mixed RKM
      outmix = cluspcamix(mixdt1df, nclus = 3, ndim = 2)
      mixsimresmat2[j,(72*(l-1)+(i-1)*8+7)] <- adjustedRandIndex(outmix$cluster, mixdt1$id)
      cat('Mixed RKM done for iteration',j,'of',i,'\n')
      # Mixed FKM
      outfmix = cluspcamix(mixdt1df, nclus = 3, ndim = 2, method = 'mixedFKM')
      mixsimresmat2[j,(72*(l-1)+(i-1)*8+8)] <- adjustedRandIndex(outfmix$cluster, mixdt1$id)
      cat('Mixed FKM done for iteration',j,'of',i,'\n')
      cat('Iteration',j,'of',i,'for data set',l+1,'complete \n')
    }
    save(mixsimresmat2, file='mixsimresmat2.RData')
  }
}