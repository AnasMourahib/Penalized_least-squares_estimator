starting_point <-function(data, nrcol, quant = 0.9){
  N <- nrow(data)

  dataP <- apply(data, 2, function(i) N/(N + 0.5 - rank(i))) 
  U <- quantile(rowSums(dataP),quant) 
  dataU <-dataP[rowSums(dataP)>U,] 
  ndata <- t(apply(dataU,1, function(i) i/sum(i))) 
  
  kmean <- kmeans(ndata,centers=nrcol,nstart=5) 
  startk <- sapply(c(1:nrcol), function(j) kmean$centers[j,]*kmean$size[j])
  resk <- t(apply(startk, 1, function(x) x/sum(x)))
  return(resk)
}
#starting_point<-function(N, r,k, d, A, alpha){
#  data<-N_generate_Mix_log(N,A,alpha)
#  dataP <- apply(data, 2, function(i) N/(N + 0.5 - rank(i))) #transform to unit Pareto
#  U1 <- quantile(rowSums(dataP),0.9) #pick a high quantile
#  dataU1 <-dataP[rowSums(dataP)>U1,] #select datapoints exceeding this quantile
#  ndata1 <- t(apply(dataU1,1, function(i) i/sum(i))) #normalized to unit simplex
#  kmean1 <- kmeans(ndata1,centers=k,nstart=5) #for an r-column max-linear parameter matrix
#  startall1<-matrix(0, nrow=d, ncol=k)
#  for (j in 1:k){
#    startall1[,j]=kmean1$centers[j,]*(d/nrow(ndata1))*kmean1$size[j]
#  }
#  return(t(apply(startall1, 1, function(x) x/sum(x))))
#}



