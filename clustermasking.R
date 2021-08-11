require(RColorBrewer)
require(mvtnorm)
require(factoextra)

set.seed(1)
n <- 50

# Create clusters in 2 dimensions
sigma <- matrix(c(1,0,0,1),2,2)
c1 <- rmvnorm(n, mean=c(-5,0), sigma=sigma)
c2 <- rmvnorm(n, mean=c(0,5), sigma=sigma)
c3 <- rmvnorm(n, mean=c(5,0), sigma=sigma)

# Visualise clusters
plot(c1, ylim=c(-2,7), xlim=c(-7, 7), pch=16, col='#E41A1C', xlab='Dimension 1', ylab='Dimension 2',
     main='Underlying Clustering Structure in 2 dimensions')
points(c2, pch=16, col='#377EB8')
points(c3, pch=16, col='#4DAF4A')
legend("topleft", legend=c("Cluster 1", "Cluster 2", "Cluster 3"),
       col=c("#E41A1C", "#377EB8","#4DAF4A"), pch=16, text.font=2)

# Define data frames with first 2 dimensions
c1df <- as.data.frame(c1)
c2df <- as.data.frame(c2)
c3df <- as.data.frame(c3)
clustersdf <- rbind(c1df, c2df)
clustersdf <- rbind(clustersdf, c3df)

# Additional columns/dimensions
for (i in 1:4){
  set.seed(i)
  newcol <- rnorm(3*n, mean=0, sd=5)
  newcol <- as.data.frame(newcol)
  clustersdf <- cbind(clustersdf, newcol)
}

# PCA
clusterspca <- prcomp(clustersdf, scale. = TRUE)

# scree plot
barplot(summary(clusterspca)$importance[c(2,5,8,11,14,17)],
        main='Proportion of variance explained by principal components',
        xlab='Principal Components',
        ylab='Proportion of variance (%) explained',
        col='cyan3',
        names.arg = c('PC1','PC2','PC3','PC4','PC5','PC6'))

# eigenvalues
get_eigenvalue(clusterspca)

# SVD
X <- as.matrix(clustersdf)
svdX <- svd(scale(X))

# Matrix Recovered from SVD
svdX$u%*%diag(svdX$d)%*%t(svdX$v)

# Plot using first 2 PCs
u2 <- svdX$u[,1:2]
d2 <- diag(svdX$d[1:2])
v2 <- svdX$v[,1:2]
X2 <- u2%*%d2%*%t(v2)
plot(X2[,1],X2[,2])

# Plot using first 3 PCs
u3 <- svdX$u[,1:3]
d3 <- diag(svdX$d[1:3])
v3 <- svdX$v[,1:3]
X3 <- u3%*%d3%*%t(v3)
plot(X3[,1],X3[,2])

# Plot using first 4 PCs
u4 <- svdX$u[,1:4]
d4 <- diag(svdX$d[1:4])
v4 <- svdX$v[,1:4]
X4 <- u4%*%d4%*%t(v4)
plot(X4[,1],X4[,2])

# Plot using first 5 PCs
u5 <- svdX$u[,1:5]
d5 <- diag(svdX$d[1:5])
v5 <- svdX$v[,1:5]
X5 <- u5%*%d5%*%t(v5)
plot(X5[,1],X5[,2])

# Plot using all 6 PCs
u6 <- svdX$u[,1:6]
d6 <- diag(svdX$d[1:6])
v6 <- svdX$v[,1:6]
X6 <- u6%*%d6%*%t(v6)
plot(X6[,1],X6[,2])

# Visualise clusters in first 2 dimensions when retaining 2 PCs
palette(brewer.pal(9,'Set1'))
k2 <- kmeans(X2, 3, nstart = 25, iter.max=1000)
plot(X2[,1],X2[,2], col=k2$cluster, pch=16)
plot(x=X[,1], y=X[,2], col=k2$cluster, pch=16,
     xlab='Dimension 1', ylab='Dimension 2',
     main='k-Means Clustering based on first 2 principal component scores')
legend("topleft", legend=c("Cluster 1", "Cluster 2", "Cluster 3"),
       col=c("#E41A1C", "#377EB8","#4DAF4A"), pch=16, text.font=2)

# Visualise clusters in first 2 dimensions when retaining 3 PCs
k3 <- kmeans(X3, 3, nstart = 25, iter.max=1000)
plot(X3[,1],X3[,2], col=k3$cluster, pch=16)
plot(x=X[,1], y=X[,2], col=k3$cluster, pch=16,
     xlab='Dimension 1', ylab='Dimension 2',
     main='k-Means Clustering based on first 3 principal component scores')
legend("topleft", legend=c("Cluster 1", "Cluster 2", "Cluster 3"),
       col=c("#E41A1C", "#377EB8","#4DAF4A"), pch=16, text.font=2)

# Visualise clusters in first 2 dimensions when retaining 4 PCs
k4 <- kmeans(X4, 3, nstart = 25, iter.max=1000)
plot(X4[,1],X4[,2], col=k4$cluster, pch=16)
plot(x=X[,1], y=X[,2], col=k4$cluster, pch=16,
     xlab='Dimension 1', ylab='Dimension 2',
     main='k-Means Clustering based on first 4 principal component scores')
legend("topleft", legend=c("Cluster 1", "Cluster 2", "Cluster 3"),
       col=c("#E41A1C", "#377EB8","#4DAF4A"), pch=16, text.font=2)

# Visualise clusters in first 2 dimensions when retaining 5 PCs
k5 <- kmeans(X5, 3, nstart = 25, iter.max=1000)
plot(X5[,1],X5[,2], col=k5$cluster, pch=16)
plot(x=X[,1], y=X[,2], col=k5$cluster, pch=16,
     xlab='Dimension 1', ylab='Dimension 2',
     main='k-Means Clustering based on first 5 principal component scores')
legend("topleft", legend=c("Cluster 1", "Cluster 2", "Cluster 3"),
       col=c("#E41A1C", "#377EB8","#4DAF4A"), pch=16, text.font=2)

# Visualise clusters in first 2 dimensions when retaining all 6 PCs
k6 <- kmeans(X6, 3, nstart = 25, iter.max=1000)
plot(X6[,1],X6[,2], col=k6$cluster, pch=16)
plot(x=X[,1], y=X[,2], col=k6$cluster, pch=16)