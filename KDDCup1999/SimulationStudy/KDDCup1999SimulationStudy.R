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
require(caret)

# Load data set and convert to data frame
load('kddcup99.Rdata')
kddcup99 <- as.data.frame(kddcup99)

# Remove duplicates
kddcup99 <- unique(kddcup99)

# Remove redundant columns with just 1 value
kddcup99 <- kddcup99[,c(1:19, 22:ncol(kddcup99))]

# Define continuous & categorical columns
numericcols <- c(1, 5, 6, 10, 13, 16, 17, 19, 21:39)
categoricalcols <- setdiff(1:ncol(kddcup99sample), numericcols)

# Convert numeric columns to numeric and categorical to factors
for (i in numericcols){
  kddcup99[,i] <- as.numeric(kddcup99[,i])
}

for (i in categoricalcols){
  kddcup99[,i] <- as.factor(kddcup99[,i])
}

# Preview data set variables
str(kddcup99)

# Obtain data set sample
set.seed(2)
train_idx <- createDataPartition(kddcup99$label, times = 1, p = 0.1, list=F)
kddcup99sample <- kddcup99[train_idx,]

# Define labels for network intrusion types
label <- c()
u2r <- c('buffer_overflow','loadmodule','perl','rootkit')
dos <- c('back','land','neptune','pod','smurf','teardrop')
r21 <- c('ftp_write','guess_passwd','imap','multihop','phf','spy','warezclient', 'warezmaster')
probe <- c('ipsweep','portsweep','satan','nmap')

# Convert intrusion type to its label
for (i in 1:nrow(kddcup99sample)){
  if (kddcup99sample[i,40]=="normal"){
    label <- c(label,"normal")
  }
  else if (kddcup99sample[i,40] %in% u2r){
    label <- c(label,"u2r")
  }
  else if (kddcup99sample[i,40] %in% dos){
    label <- c(label,"dos")
  }
  else if (kddcup99sample[i,40] %in% r21){
    label <- c(label,"r21")
  }
  else{
    label <- c(label,"probe")
  }
}

# Replace label column & convert to factor
kddcup99sample$label <- label
kddcup99sample$label <- as.factor(kddcup99sample$label)

# Preview final data set
str(kddcup99sample)

# You may load the data set sample from here if you want to skip the parts above
load('kddcup99sample.Rdata')

# Compute Gower distance
gower_dist_kddcup9 <- daisy(kddcup99sample[,-40], metric = "gower")
gower_mat_kddcup9 <- as.matrix(gower_dist_kddcup9)

# Implement PAM
pam_fit_kddcup9 <- pam(gower_dist_kddcup9, diss = TRUE, k = 5)

### Extra: Calculate Average Silhouette Width
sil_width_kddcup9 <- pam_fit_kddcup9$silinfo$avg.width  
sil_width_kddcup9

# Print cluster size
table(pam_fit_kddcup9$clustering)

# Calculate ARI
adjustedRandIndex(pam_fit_kddcup9$clustering, kddcup99sample[,40])

# K-prototypes
outk_kddcup9 <- kproto(kddcup99sample[,-40], 5)

# Print cluster size
table(outk_kddcup9$cluster)

# Calculate ARI
adjustedRandIndex(outk_kddcup9$cluster, kddcup99sample[,40])

# Mixed K-means
# Define categorical columns
categoricalcols <- setdiff(1:39, numericcols)

# Distance options for mixed variables data set: gower, wishart, podani, huang, harikumar, ahmad 
# We use the Ahmad distance
mix_kddcup9 <- distmix(kddcup99sample[,-40], method = "ahmad", idcat = categoricalcols, idnum = numericcols)

# Implement K-Medoids on resulting distances
kmedres_kddcup9 <- fastkmed(mix_kddcup9, 5, iterate = 50, init = NULL)

# Print cluster size
table(kmedres_kddcup9$cluster)

# Calculate ARI
adjustedRandIndex(kmedres_kddcup9$cluster, kddcup99sample[,40])

# Modha-Sprangler K-means

# Construct data frames with scaled contnuous and dummy coded categorical attributes
conDf_kddcup9 <- data.frame(scale(kddcup99sample[,numericcols]))
catDf_kddcup9 <- dummyCodeFactorDf(data.frame(kddcup99sample[,categoricalcols]))

# Apply Modha-Spangler K-means clustering
msRes_kddcup9 <- gmsClust(conDf_kddcup9, catDf_kddcup9, nclust = 5, searchDensity = 2)

# Calculate ARI
adjustedRandIndex(msRes_kddcup9$results$cluster, kddcup99sample[,40])

# FAMD + K-means (Tandem approach)

# Use GCV to calculate number of dimensions to be retained
howmany <- estim_ncp(data.matrix(kddcup99sample[,-40]))$ncp

# Apply FAMD
outpcamix_kddcup9 <- FAMD(kddcup99sample[,-40], ncp = howmany, graph=FALSE)
outkm_kddcup9 <- kmeans(outpcamix_kddcup9$ind$coord, 5)

# Calculate ARI
adjustedRandIndex(outkm_kddcup9$cluster, kddcup99sample[,40])


# NLPCA + K-Means
nlpcakddcup9 <- homals(kddcup99sample, active=c(rep(TRUE, 39), FALSE), verbose=1, rank=1,
                       level=c('numerical', rep('nominal',3), rep('numerical',2),
                               'nominal', rep('ordinal',2), 'numerical',
                               'ordinal', 'nominal', 'numerical', rep('nominal',2),
                               rep('numerical',2), 'ordinal', 'numerical',
                               'nominal', rep('numerical',19), 'nominal'), ndim = howmany)

nlpcakddcup9kmeans <- kmeans(nlpcakddcup9$objscores, 5)

# Calculate ARI
adjustedRandIndex(nlpcakddcup9kmeans$cluster, kddcup99sample[,40])

# Mixed Reduced K-means

outmix_kddcup9 <- cluspcamix(kddcup99sample[,-40], nclus = 5, ndim = 2)

# Calculate ARI
adjustedRandIndex(outmix_kddcup9$cluster, kddcup99sample[,40])

# Mixed Factorial K-means

outfmix_kddcup9 = cluspcamix(kddcup99sample[,-40], nclus = 5, ndim = 2, method = 'mixedFKM')

#Calculate ARI
adjustedRandIndex(outfmix_kddcup9$cluster, kddcup99sample[,40])

#### SIMULATION STUDY- 50 ITERATIONS

niter <- 50
kddcup99aris <- matrix(NA, nrow = niter, ncol = 8)

for (i in 1:niter){
  # PAM with Gower's
  pam_fit_kddcup9 <- pam(gower_dist_kddcup9, diss = TRUE, k = 5)
  kddcup99aris[i,1] <- adjustedRandIndex(pam_fit_kddcup9$clustering, kddcup99sample[,40])
  cat('PAM done for iteration',i,'\n')
  # K-prototypes
  outk_kddcup9 <- kproto(kddcup99sample[,-40], 5)
  kddcup99aris[i,2] <- adjustedRandIndex(outk_kddcup9$cluster, kddcup99sample[,40])
  cat('K-prototypes done for iteration',i,'\n')
  # Mixed K-means
  kmedres_kddcup9 <- fastkmed(mix_kddcup9, 5, iterate = 50, init = NULL)
  kddcup99aris[i,3] <- adjustedRandIndex(kmedres_kddcup9$cluster, kddcup99sample[,40])
  cat('Mixed K-means done for iteration',i,'\n')
  # Modha-Spangler K-Means
  msRes_kddcup9 <- gmsClust(conDf_kddcup9, catDf_kddcup9, nclust = 5, searchDensity = 2)
  kddcup99aris[i,4] <- adjustedRandIndex(msRes_kddcup9$results$cluster, kddcup99sample[,40])
  cat('Modha-Spangler done for iteration',i,'\n')
  # FAMD + K-means
  outkm_kddcup9 <- kmeans(outpcamix_kddcup9$ind$coord, 5)
  kddcup99aris[i,5] <- adjustedRandIndex(outkm_kddcup9$cluster, kddcup99sample[,40])
  cat('FAMD done for iteration',i,'\n')
  # NLPCA + K-means
  nlpcakddcup9kmeans <- kmeans(nlpcakddcup9$objscores, 5)
  kddcup99aris[i,6] <- adjustedRandIndex(nlpcakddcup9kmeans$cluster, kddcup99sample[,40])
  cat('NLPCA done for iteration',i,'\n')
  # Mixed RKM
  outmix_kddcup9 = cluspcamix(kddcup99sample[,-40], nclus = 5, ndim = 2)
  kddcup99aris[i,7] <- adjustedRandIndex(outmix_kddcup9$cluster, kddcup99sample[,40])
  cat('Mixed RKM done for iteration',i,'\n')
  # Mixed FKM
  outfmix_kddcup9 = cluspcamix(kddcup99sample[,-40], nclus = 5, ndim = 2, method = 'mixedFKM')
  kddcup99aris[i,8] <- adjustedRandIndex(outfmix_kddcup9$cluster, kddcup99sample[,40])
  cat('Mixed FKM done for iteration',i,'\n')
  cat('Iteration',i,'complete \n')
  # Uncomment line below if you want to save the file after every iteration
  #save(kddcup99aris, file='kddcup99aris.RData')
}
# If you want to skip all simulations above, simulation results file can be found below
load('kddcup99aris.Rdata')
# Box plot of simulations results
boxplot(kddcup99aris,
        names=c('PAM', 'K-PROT', 'MIX', 'MS', 'FAMD', 'NLPCA', 'RKM', 'FKM'),
        ylab='Adjusted Rand Index (ARI)',
        col=c(rep('khaki1',4), rep('darkolivegreen3', 2), rep('deepskyblue1',2)),
        main='ARI Values for clustering methods on KDD Cup 1999 data set')
