
glmblock.fit <- function(Adj_list, y, 
                         family = "gaussian", intercept = T,
                         Z, gamma = 0, lambda = 0) {
  #' 
  #' Function to perform block penalized regression using a specified community membership
  #' 
  #' @param Adj_list a list with adjacency matrices of the same size n x n
  #' @param y a vector of responses
  #' @param family response type, either "gaussian" for continuous responses, or "binomial" for classification
  #' @param intercept indicate whether intercept should be fitted
  #' @param Z membership matrix of size n by K specifying the community assignments
  #' @param gamma ridge penalty parameter
  #' @param lambda lasso penalty parameter
  #' 
  #' @return a list containing the fitted coefficents matrix \code{B}, intercept \code{b}, 
  #' and the value of the objective function and classification error (if binomial).
  #' 
  #' @references 
  #' 
  #' @author Jes\'us Arroyo <jarroyor@umich.edu>
  # fits B and b in the model
  # l(B,b) + lambda* ||weights*B||_1 + gamma/2 * ||weights*B||_F^2
  # subject to B = ZCZ^T
  Z = Matrix(Z)
  n = ncol(Adj_list[[1]])
  m = length(Adj_list)
  ZAZ = lapply(Adj_list, function(A) Matrix::crossprod(Z, Matrix::crossprod(A,Z)))
  ZAZr <- as.matrix(t(sapply(ZAZ, function(U) {
    U = 2*U
    diag(U) = diag(U)/2
    uptri(U, diag =  T)
  })))  
  
  weights = 2*(tcrossprod(diag(crossprod(Z)))  - crossprod(Z))
  diag(weights) = diag(weights)/2
  weights = uptri(weights, diag = T)
  weights = ifelse(weights==0, 1, weights)
  
  
  Alpha = lambda / (gamma + lambda + 1*((lambda+gamma)==0))
  Lambda = lambda + gamma#* sum(weights) /ncol(ZAZr)
  glmfit <- glmnet(x= ZAZr, y = y, family = family, 
                   alpha = Alpha,
                   lambda = Lambda*sum(weights)/length(weights),
                   penalty.factor = weights*length(weights)/sum(weights), 
                   standardize = FALSE)
  glmfitbeta = glmfit$beta
  b <- glmfit$a0
  
  P_hatr <- glmfitbeta
  
  Phat <- matrix(0, ncol = ncol(Z), nrow = ncol(Z))
  Phat[upper.tri(Phat, diag=T)] <- as.double(P_hatr)
  Phat <- Phat + t(Phat)
  diag(Phat) = diag(Phat)/2
  Class_error = NULL
  if(family=="gaussian") {
    responses = predict(glmfit, newx = ZAZr)
    loss = sum((y - responses)^2)/(2*length(y))
  }
  if(family=="binomial") {
    responses = as.numeric(predict(glmfit, newx = ZAZr, type = "class"))
    loss = -sum((y==1)*responses-log(1+exp(responses)))/(length(y))
    Class_hat = max(y)*(responses >=0) + min(y)*(responses < 0)
    Class_error = sum(Class_hat != y) / length(y)
  }
  weights <- tcrossprod(diag(crossprod(Z))) - crossprod(Z)
  penalty <- lambda * sum(abs(weights*Phat)) + gamma/2 * sum(weights*(Phat)^2)
  objective <- loss + penalty
  #diag(Phat) <- diag(Phat) / 2
  return(list(B= Z%*%(Phat%*%t(Z)), b = b, 
              objective = objective, Class_error = Class_error))
}


Z_membership_from_list <- function(community_list) {
  require(Matrix)
  K = length(community_list)
  n = max(unlist(community_list))
  Z = Matrix(0, n, K)
  for(k in 1:length(community_list)) {
    Z[community_list[[k]], k] = 1
  }
  return(Z)
}

Z_membership_from_vector <- function(community_vec) {
  require(Matrix)
  n = length(community_vec)
  com_labels = unique(community_vec)
  K = length(com_labels)
  Z = Matrix(0, n, K)
  for(k in 1:K) {
    Z[which(community_vec == com_labels[k]), k] = 1
  }
  return(Z)
}


Z_membership_to_list <- function(Z) {
  require(Matrix)
  n = nrow(Z)
  K = ncol(Z)
  communities = list()
  for(k in 1:K) {
    communities[[k]] = which(Z[,k]>0)
  }
  return(communities)
}