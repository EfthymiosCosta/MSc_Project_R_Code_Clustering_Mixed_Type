# Discretisation function

discretise <- function(vec, nlevels){
  range <- max(vec)-min(vec)
  width <- range/nlevels
  vecnew <- c()
  for (i in 1:length(vec)){
    for (j in 1:nlevels){
      if (vec[i]<= min(vec)+j*width & vec[i]>min(vec)+(j-1)*width){
        val <- j
      }
      if (vec[i]==min(vec)){
        val <- 1
      }
      if (vec[i]==max(vec)){
        val <- nlevels
      }
    }
    vecnew <- c(vecnew, val)
  }
  vecnew <- as.factor(vecnew)
}