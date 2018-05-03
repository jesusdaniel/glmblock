##

sbm <- function(p,q, K=4, nk=10) {
  Psi = (p-q)*diag(K)  + q
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)), W = A))
}



sbm_11 <- function(p,q, K=4, nk=10) {
  Psi = (p-q)*diag(K)  + q
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  2*(Aobs + t(Aobs))-1
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)), W = A))
}

erdos_renyi_1_2 <- function(N) {
  edges = 2*rbinom(N*(N-1)/2, 1, prob = 0.5) -1
  A = matrix(0, nrow = N, ncol = N)
  A[upper.tri(A)] = edges
  A = A + t(A)
  return(A)
}

spectral_clustering <- function(A, K=4) {
  #browser()
  eig <- eigen(A)
  largest = order(abs(eig$values),decreasing = T)
  lat_pos <- eig$vectors[,largest[1:K]] %*% diag(sqrt(abs(eig$values[largest[1:K]])))
  clust <- kmeans(lat_pos, centers = K, nstart = 100)
  Z <- t(sapply(clust$cluster, function(k) {
    u = rep(0, K)
    u[k] = 1
    return(u)}))
  return(Z)
}

sbm_normal <- function(p,q, K=4, nk=10) {
  Psi = (p-q)*diag(K)  + q
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rnorm(n = 1, mean = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)), W = A))
}



generate_gaussian_blockmodel <- function(Nc = 30, K = 2, rho = 0, U = NULL) {
  A <- matrix(0, nrow = Nc*K, ncol = Nc*K)
  Asub <- matrix(0, nrow = Nc, ncol = Nc)
  #Generate intra-community edges
  Ne = (Nc^2 - Nc)/2
  Sigma = (1-rho) * diag(Ne) + rho
  if(is.null(U))
    U = chol(Sigma)
  for(i in 1:K) {
    Asub[upper.tri(Asub)] = rnorm(Ne)%*%U
    A[((i-1)*Nc+1):(i*Nc), ((i-1)*Nc+1):(i*Nc)] = Asub
  }
  # generate inter-commuinities edges
  if(K > 1) {
    Ne = Nc^2
    Sigma = (1-rho) * diag(Ne) + rho
    U = chol(Sigma)
    for(i in 1:(K-1)) {
      for(j in (i+1):K) {
        Asub = rnorm(Ne) %*% U
        A[((i-1)*Nc+1):(i*Nc), ((j-1)*Nc+1):(j*Nc)] = Asub
      }
    }  
  }
  A = A + t(A)
  return(A)
}

gaussian_blockmodel_cov_inv <- function(Nc = 30, K = 2, rho = 0) {
  Sigma_AA = array(0, dim =  c(Nc*K, Nc*K, Nc*K, Nc*K))
  Sigma_AA_inv = Sigma_AA
  for(i in 1:K) {
    for(j in 1:K) {
      Sigma_AA[((i-1)*Nc+1):(i*Nc), ((j-1)*Nc+1):(j*Nc), ((i-1)*Nc+1):(i*Nc), ((j-1)*Nc+1):(j*Nc)] = rho
      Sigma_AA_inv[((i-1)*Nc+1):(i*Nc), ((j-1)*Nc+1):(j*Nc), ((i-1)*Nc+1):(i*Nc), ((j-1)*Nc+1):(j*Nc)] = 
        -rho/((1-rho)^2 + rho*(1-rho)*(Nc^2- (i==j)*Nc))
    }
  }
  require(gtools)
  Sigma_AA[cbind(permutations(K*Nc, 2,repeats.allowed=T), permutations(K*Nc, 2,repeats.allowed=T))] = 1
  for(i in 1:(K*Nc)) {
    for(j in 1:(K*Nc)) {
      Sigma_AA_inv[i,j,i,j] = 1/(1-rho)-rho/((1-rho)^2 + rho*(1-rho)*(Nc^2-Nc*(floor((i-1)/Nc)==floor((j-1)/Nc))))
    }
  }
  for(i in 1:(Nc*K)) {
    Sigma_AA[i,i,,] = 0
    Sigma_AA[,,i,i] = 0
    Sigma_AA_inv[i,i,,] = 0
    Sigma_AA_inv[,,i,i] = 0
  }
  return(list(Sigma_AA = Sigma_AA, Sigma_AA_inv = Sigma_AA_inv))
}


eigenvector_clustering <- function(A,K=2) {
  eig <- eigen(A)
  lat_pos <- eig$vectors[,1]
  clust <- kmeans(lat_pos, centers = 2, nstart = 100)
  Z <- t(sapply(clust$cluster, function(k) {
    u = rep(0, 2)
    u[k] = 1
    return(u)}))
  return(Z)
}