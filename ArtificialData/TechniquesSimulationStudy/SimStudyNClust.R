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

# Simulation study on effect of number of clusters

niter <- 100

# Empty matrices to store simulation results
mixsimresmatnclust2 <- matrix(NA, nrow=100, ncol=32)
mixsimresmatnclust2 <- matrix(NA, nrow=100, ncol=128)

# Empty matrices to store artificial data sets constructed & their labels
fulldtmatnclust <- matrix(NA, nrow=480, ncol=40)
fullidmatnclust <- matrix(NA, nrow=480, ncol=4)

fulldtmatnclust2 <- matrix(NA, nrow=480, ncol=160)
fullidmatnclust2 <- matrix(NA, nrow=480, ncol=16)

# Simulation for 1 artificial data set
for (i in 1:4){
  set.seed(1234)
  # Construction of artificial data set with:
    # Average overlap: 0.05
    # Maximum Overlap: 0.07
    # Number of clusters: i+2
    # Number of observations: 480
    # Number of variables: 10
  mixsimaux <- MixSim(BarOmega = 0.05, MaxOmega = 0.07, PiLow=1.0, K = i+2, p = 10, resN = 100000)
  mixdtaux <- simdataset(n = 480, Pi = mixsimaux$Pi, Mu = mixsimaux$Mu, S = mixsimaux$S)
  fullidmatnclust[,i] <- mixdtaux$id
  # Store data set labels
  save(fullidmatnclust, file='fullidmatnclust.RData')
  # Discretise first 2 attributes using 3 levels
  for (k in 1:2){
    mixdtaux$X[,k] <- discretise(mixdtaux$X[,k], nlevels=3)
  }
  mixdt1df <- as.data.frame(mixdtaux$X)
  for (k in 1:2){
    mixdt1df[,k] <- as.factor(mixdt1df[,k])
  }
  # Store final data set
  fulldtmatnclust[, ((i-1)*10+1):(i*10)] <- as.matrix(mixdt1df)
  save(fulldtmatnclust, file='fulldtmatnclust.RData')
  # Calculate Gower's dissimilarities and Ahmad distances
  gower_dist <- daisy(mixdt1df, metric = "gower")
  mix_sim <- distmix(mixdt1df, method = "ahmad", idcat = c(1:2), idnum = c(3:10))
  # Specify continuous & categorical attributes
  conDf <- data.frame(scale(mixdt1df[,3:10]))
  catDf <- dummyCodeFactorDf(data.frame(mixdt1df[,1:2]))
  # Number of dimensions to be retained using GCV
  howmany <- estim_ncp(data.matrix(mixdt1df))$ncp
  if (howmany<=1){
    howmany <- 2
  }
  outpcamix <- FAMD(mixdt1df, ncp = howmany, graph=FALSE)
  nlpca <- homals(mixdt1df, active=TRUE, verbose=1, rank=1,
                  level=c(rep('nominal', 2), rep('numerical', 8)), ndim = howmany)
  for (j in 1:niter){
    # PAM with Gower's
    pam_fit <- pam(gower_dist, diss = TRUE, k = i+2, do.swap = FALSE, cluster.only = TRUE,
                   medoids=sample(1:480, i+2))
    mixsimresmatnclust[j,((i-1)*8+1)] <- adjustedRandIndex(pam_fit, mixdtaux$id)
    cat('PAM done for iteration',j,'of',i,'\n')
    # K-prototypes
    outk <- kproto(mixdt1df, i+2)
    mixsimresmatnclust[j,((i-1)*8+2)] <- adjustedRandIndex(outk$cluster, mixdtaux$id)
    cat('K-prototypes done for iteration',j,'of',i,'\n')
    # Mixed K-means
    kmedres <- fastkmed(mix_sim, i+2, iterate = 100, init = sample(1:nrow(mixdt1df), i+2))
    mixsimresmatnclust[j,((i-1)*8+3)] <- adjustedRandIndex(kmedres$cluster, mixdtaux$id)
    cat('Mixed K-means done for iteration',j,'of',i,'\n')
    # Modha-Spangler
    msRes <- gmsClust(conDf, catDf, nclust = i+2, searchDensity = 5)
    mixsimresmatnclust[j,((i-1)*8+4)] <- adjustedRandIndex(msRes$results$cluster, mixdtaux$id)
    cat('Modha-Spangler done for iteration',j,'of',i,'\n')
    # FAMD + K-means
    outkm <- kmeans(outpcamix$ind$coord, i+2, iter.max=100000, algorithm = "MacQueen")
    mixsimresmatnclust[j,((i-1)*8+5)] <- adjustedRandIndex(outkm$cluster, mixdtaux$id)
    cat('FAMD done for iteration',j,'of',i,'\n')
    # NLPCA + K-means
    nlpcakmeans <- kmeans(nlpca$objscores, i+2, iter.max=100000, algorithm = "MacQueen")
    mixsimresmatnclust[j,((i-1)*8+6)] <- adjustedRandIndex(nlpcakmeans$cluster, mixdtaux$id)
    cat('NLPCA done for iteration',j,'of',i,'\n')
    # Mixed RKM
    outmix = cluspcamix(mixdt1df, nclus = i+2, ndim = 2)
    mixsimresmatnclust[j,((i-1)*8+7)] <- adjustedRandIndex(outmix$cluster, mixdtaux$id)
    cat('Mixed RKM done for iteration',j,'of',i,'\n')
    # Mixed FKM
    outfmix = cluspcamix(mixdt1df, nclus = i+2, ndim = 2, method = 'mixedFKM')
    mixsimresmatnclust[j,((i-1)*8+8)] <- adjustedRandIndex(outfmix$cluster, mixdtaux$id)
    cat('Mixed FKM done for iteration',j,'of',i,'\n')
    cat('Iteration',j,'of',i,'complete \n')
  }
  save(mixsimresmatnclust, file='mixsimresmatnclust.RData')
}


# Simulation for 4 different artificial data sets
for (l in 1:4){
  set.seed(1234+l)
  for (i in 1:4){
    # Construct artificial data set
    mixsimaux <- MixSim(BarOmega = 0.05, MaxOmega = 0.07, PiLow=1.0, K = i+2, p = 10, resN = 100000)
    mixdtaux <- simdataset(n = 480, Pi = mixsimaux$Pi, Mu = mixsimaux$Mu, S = mixsimaux$S)
    fullidmatnclust2[,(4*(l-1)+i)] <- mixdtaux$id
    # Store data set labels
    save(fullidmatnclust2, file='fullidmatnclust2.RData')
    # Discretise first 2 attributes using 3 levels
    for (k in 1:2){
      mixdtaux$X[,k] <- discretise(mixdtaux$X[,k], nlevels=3)
    }
    mixdt1df <- as.data.frame(mixdtaux$X)
    for (k in 1:2){
      mixdt1df[,k] <- as.factor(mixdt1df[,k])
    }
    # Store final data set
    fulldtmatnclust2[, (40*(l-1)+(i-1)*10+1):(40*(l-1)+(i*10))] <- as.matrix(mixdt1df)
    save(fulldtmatnclust2, file='fulldtmatnclust2.RData')
    # Calculate Gower's dissimilarities and Ahmad distances
    gower_dist <- daisy(mixdt1df, metric = "gower")
    mix_sim <- distmix(mixdt1df, method = "ahmad", idcat = c(1:2), idnum = c(3:10))
    # Specify continuous & categorical attributes
    conDf <- data.frame(scale(mixdt1df[,3:10]))
    catDf <- dummyCodeFactorDf(data.frame(mixdt1df[,1:2]))
    # Number of dimensions to be retained using GCV
    howmany <- estim_ncp(data.matrix(mixdt1df))$ncp
    if (howmany<=1){
      howmany <- 2
    }
    outpcamix <- FAMD(mixdt1df, ncp = howmany, graph=FALSE)
    nlpca <- homals(mixdt1df, active=TRUE, verbose=1, rank=1,
                    level=c(rep('nominal', 2), rep('numerical', 8)), ndim = howmany)
    for (j in 1:niter){
      # PAM with Gower's
      pam_fit <- pam(gower_dist, diss = TRUE, k = i+2, do.swap = FALSE, cluster.only = TRUE,
                     medoids=sample(1:480, i+2))
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+1)] <- adjustedRandIndex(pam_fit, mixdtaux$id)
      cat('PAM done for iteration',j,'of',i,'\n')
      # K-prototypes
      outk <- kproto(mixdt1df, i+2)
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+2)] <- adjustedRandIndex(outk$cluster, mixdtaux$id)
      cat('K-prototypes done for iteration',j,'of',i,'\n')
      # Mixed K-means
      kmedres <- fastkmed(mix_sim, i+2, iterate = 100, init = sample(1:nrow(mixdt1df), i+2))
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+3)] <- adjustedRandIndex(kmedres$cluster, mixdtaux$id)
      cat('Mixed K-means done for iteration',j,'of',i,'\n')
      # Modha-Spangler
      msRes <- gmsClust(conDf, catDf, nclust = i+2, searchDensity = 5)
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+4)] <- adjustedRandIndex(msRes$results$cluster, mixdtaux$id)
      cat('Modha-Spangler done for iteration',j,'of',i,'\n')
      # FAMD + K-means
      outkm <- kmeans(outpcamix$ind$coord, i+2, iter.max=100000, algorithm = "MacQueen")
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+5)] <- adjustedRandIndex(outkm$cluster, mixdtaux$id)
      cat('FAMD done for iteration',j,'of',i,'\n')
      # NLPCA + K-means
      nlpcakmeans <- kmeans(nlpca$objscores, i+2, iter.max=100000, algorithm = "MacQueen")
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+6)] <- adjustedRandIndex(nlpcakmeans$cluster, mixdtaux$id)
      cat('NLPCA done for iteration',j,'of',i,'\n')
      # Mixed RKM
      outmix = cluspcamix(mixdt1df, nclus = i+2, ndim = 2)
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+7)] <- adjustedRandIndex(outmix$cluster, mixdtaux$id)
      cat('Mixed RKM done for iteration',j,'of',i,'\n')
      # Mixed FKM
      outfmix = cluspcamix(mixdt1df, nclus = i+2, ndim = 2, method = 'mixedFKM')
      mixsimresmatnclust2[j,(32*(l-1)+(i-1)*8+8)] <- adjustedRandIndex(outfmix$cluster, mixdtaux$id)
      cat('Mixed FKM done for iteration',j,'of',i,'\n')
      cat('Iteration',j,'of',i,'for data set',l+1,'complete \n')
    }
    save(mixsimresmatnclust2, file='mixsimresmatnclust2.RData')
  }
}