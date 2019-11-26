uptri <- function(network, diag = F)  {
  network[upper.tri(network, diag = diag)]
}
to_matrix <- function(vec) {
  Nodes = (sqrt(8*length(vec) +1) + 1)/2
  network = matrix(0,nrow = Nodes, ncol = Nodes)
  network[upper.tri(network)] = as.double(vec)
  network = network + t(network)
  return(network)
}


spectral_clustering <- function(A, K=4) {
  require(irlba)
  #browser()
  #eig <- eigs_sym(A, K)$vectors
  eig <- irlba(A, nv = K)$u
  clust <- kmeans(eig, centers = K, nstart = 2000)
  Z <- t(sapply(clust$cluster, function(k) {
    u = rep(0, K)
    u[k] = 1
    return(u)}))
  return(Z)
}


spectral_clustering2 <- function(A, K=4) {
  require(rARPACK)
  eig <- eigs_sym(A, K)$vectors
  require(mclust)
  clust <- Mclust(eig, G = K)
    #kmeans(eig, centers = K, nstart = 20)
  Z <- t(sapply(clust$classification, function(k) {
    u = rep(0, K)
    u[k] = 1
    return(u)}))
  return(Z)
}
