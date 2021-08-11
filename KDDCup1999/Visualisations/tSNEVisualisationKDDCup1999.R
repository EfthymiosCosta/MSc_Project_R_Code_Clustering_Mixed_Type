require(Rtsne)
require(rgl)

load('kddcup99sample.Rdata')

# Set colours
cols<- c()

for (i in 1:nrow(kddcup99sample2)){
  if (kddcup99sample[i,40]=="normal"){
    cols <- c(cols,"red")
  }
  else if (kddcup99sample[i,40]=="u2r"){
    cols <- c(cols,"blue")
  }
  else if (kddcup99sample[i,40]=="dos"){
    cols <- c(cols,"forestgreen")
  }
  else if (kddcup99sample[i,40]=="r21"){
    cols <- c(cols,"black")
  }
  else{
    cols <- c(cols,"purple")
  }
}

# t-SNE applied on KDD Cup 1999 data set in 2 dimensions
tsne_kdd2d <- Rtsne(kddcup99sample, is_distance=FALSE, dims=2, theta=0.1)

# Plot in 2 dimensions
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=cols, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main='2D t-SNE plot of KDD Cup 1999 data set')
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)


# t-SNE applied on KDD Cup 1999 data set in 3 dimensions
tsne_kdd3d <- Rtsne(kddcup99sample, is_distance=FALSE, dims=3, theta=0.1)
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=cols, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# Gower's dissimilarity & PAM

load('pam_fit_kddcup9.RData')

colspam <- c()

for (i in 1:nrow(kddcup99sample)){
  if (pam_fit_kddcup9[i]==1){
    colspam <- c(colspam,"red")
  }
  else if (pam_fit_kddcup9[i]==3){
    colspam <- c(colspam,"blue")
  }
  else if (pam_fit_kddcup9[i]==4){
    colspam <- c(colspam,"forestgreen")
  }
  else if (pam_fit_kddcup9[i]==5){
    colspam <- c(colspam,"black")
  }
  else{
    colspam <- c(colspam,"purple")
  }
}

# 2D t-SNE Plot for Gower's dissimilarity with PAM clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colspam, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of PAM with Gower's dissimilarity clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for Gower's dissimilarity with PAM clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colspam, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# K-Prototypes

load('outk_kddcup9.RData')

colskprot <- c()

for (i in 1:nrow(kddcup99sample)){
  if (outk_kddcup9$cluster[i]==3){
    colskprot <- c(colskprot,"red")
  }
  else if (outk_kddcup9$cluster[i]==2){
    colskprot <- c(colskprot,"blue")
  }
  else if (outk_kddcup9$cluster[i]==4){
    colskprot <- c(colskprot,"forestgreen")
  }
  else if (outk_kddcup9$cluster[i]==1){
    colskprot <- c(colskprot,"black")
  }
  else{
    colskprot <- c(colskprot,"purple")
  }
}

# 2D t-SNE Plot for K-Prototypes clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colskprot, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of K-prototypes clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for K-Prototypes clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colskprot, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# Mixed K-Means

load('kmedres_kddcup9.RData')

colsmix <- c()

for (i in 1:nrow(kddcup99sample)){
  if (kmedres_kddcup9$cluster[i]==5){
    colsmix <- c(colsmix,"red")
  }
  else if (kmedres_kddcup9$cluster[i]==2){
    colsmix <- c(colsmix,"blue")
  }
  else if (kmedres_kddcup9$cluster[i]==3){
    colsmix <- c(colsmix,"forestgreen")
  }
  else if (kmedres_kddcup9$cluster[i]==1){
    colsmix <- c(colsmix,"black")
  }
  else{
    colsmix <- c(colsmix,"purple")
  }
}

# 2D t-SNE Plot for Mixed K-Means clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colsmix, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of Mixed K-means clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for Mixed K-Means clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colsmix, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# Modha-Spangler K-Means

load('msRes_kddcup9.RData')

colsms <- c()

for (i in 1:nrow(kddcup99sample)){
  if (msRes_kddcup9$results$cluster[i]==1){
    colsms <- c(colsms,"red")
  }
  else if (msRes_kddcup9$results$cluster[i]==2){
    colsms <- c(colsms,"blue")
  }
  else if (msRes_kddcup9$results$cluster[i]==3){
    colsms <- c(colsms,"forestgreen")
  }
  else if (msRes_kddcup9$results$cluster[i]==4){
    colsms <- c(colsms,"black")
  }
  else{
    colsms <- c(colsms,"purple")
  }
}

# 2D t-SNE Plot for Modha-Spangler K-Means clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colsms, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of Modha-Spangler K-means clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for Modha-Spangler K-Means clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colsms, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# FAMD & K-Means

load('outkm_kddcup9.RData')

colsfamd <- c()

for (i in 1:nrow(kddcup99sample)){
  if (outkm_kddcup9$cluster[i]==4){
    colsfamd <- c(colsfamd,"red")
  }
  else if (outkm_kddcup9$cluster[i]==2){
    colsfamd <- c(colsfamd,"blue")
  }
  else if (outkm_kddcup9$cluster[i]==1){
    colsfamd <- c(colsfamd,"forestgreen")
  }
  else if (outkm_kddcup9$cluster[i]==3){
    colsfamd <- c(colsfamd,"black")
  }
  else{
    colsfamd <- c(colsfamd,"purple")
  }
}

# 2D t-SNE Plot for FAMD & K-Means clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colsfamd, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of FAMD & K-means clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for FAMD & K-Means clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colsfamd, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# NLPCA & K-Means

load('nlpcakddcup9kmeans.RData')

colsnlpca <- c()

for (i in 1:nrow(kddcup99sample)){
  if (nlpcakddcup9kmeans$cluster[i]==5){
    colsnlpca <- c(colsnlpca,"forestgreen")
  }
  else if (nlpcakddcup9kmeans$cluster[i]==2){
    colsnlpca <- c(colsnlpca,"blue")
  }
  else if (nlpcakddcup9kmeans$cluster[i]==4){
    colsnlpca <- c(colsnlpca,"red")
  }
  else if (nlpcakddcup9kmeans$cluster[i]==1){
    colsnlpca <- c(colsnlpca,"black")
  }
  else{
    colsnlpca <- c(colsnlpca,"purple")
  }
}

# 2D t-SNE Plot for NLPCA & K-Means clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colsnlpca, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of NLPCA & K-means clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for NLPCA & K-Means clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colsnlpca, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# Mixed RKM

load('outmix_kddcup9.RData')

colsrkm <- c()

for (i in 1:nrow(kddcup99sample)){
  if (outmix_kddcup9$cluster[i]==1){
    colsrkm <- c(colsrkm,"red")
  }
  else if (outmix_kddcup9$cluster[i]==5){
    colsrkm <- c(colsrkm,"blue")
  }
  else if (outmix_kddcup9$cluster[i]==2){
    colsrkm <- c(colsrkm,"forestgreen")
  }
  else if (outmix_kddcup9$cluster[i]==4){
    colsrkm <- c(colsrkm,"black")
  }
  else{
    colsrkm <- c(colsrkm,"purple")
  }
}

# 2D t-SNE Plot for Mixed Reduced K-Means clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colsrkm, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of Mixed RKM clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for Mixed Reduced K-Means clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colsrkm, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')

# Mixed FKM

load('outfmix_kddcup9.RData')

colsfkm <- c()

for (i in 1:nrow(kddcup99sample)){
  if (outfmix_kddcup9$cluster[i]==1){
    colsfkm <- c(colsfkm,"red")
  }
  else if (outfmix_kddcup9$cluster[i]==5){
    colsfkm <- c(colsfkm,"blue")
  }
  else if (outfmix_kddcup9$cluster[i]==2){
    colsfkm <- c(colsfkm,"forestgreen")
  }
  else if (outfmix_kddcup9$cluster[i]==4){
    colsfkm <- c(colsfkm,"black")
  }
  else{
    colsfkm <- c(colsfkm,"purple")
  }
}

# 2D t-SNE Plot for Mixed Factorial K-Means clusters
plot(x=tsne_kdd2d$Y[,1], y=tsne_kdd2d$Y[,2], col=colsfkm, pch=16, cex=0.5,
     xlab='Dimension 1', ylab='Dimension 2', main="2D t-SNE plot of Mixed FKM clusters")
legend('topright', legend=c('Normal', 'U2R','DOS','R21','Probe'),
       col=c('red','blue','forestgreen','black','purple'), pch=16)

# 3D t-SNE Plot for Mixed Factorial K-Means clusters
plot3d(x = tsne_kdd3d$Y[,1], y = tsne_kdd3d$Y[,2], z=tsne_kdd3d$Y[,3],
       col=colsfkm, pch=16, size=1,
       xlab='Dimension 1', ylab='Dimension 2',
       zlab='Dimension 3', main='3D t-SNE plot of labels')