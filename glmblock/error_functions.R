MSE <- function(Atrain, Ytrain, Z,  Atest, Ytest) {
  training <-  glmblock.fit(Atrain, y = Ytrain, Z = Z, intercept = TRUE)
  Bhat <- training$B
  bhat <- training$b
  Yhat <- apply(Atest, 3, function(A) sum(diag(crossprod(A, Bhat)))) + bhat
  return(sum((Yhat-Ytest)^2)/length(Ytest))
}

class_error <- function(Atrain, Ytrain, Z,  Atest, Ytest, lambda = 0, gamma = 0, returnLossValue = FALSE) {
  browser()
  training <-  glmblock.fit(A = Atrain, family = "binomial", y = Ytrain, Z = Z, intercept = TRUE, 
                            lambda = lambda, gamma = gamma)
  Bhat <- training$B
  bhat <- training$b
  Yhat <- 2*(apply(Atest, 3, function(A) sum(diag(crossprod(A, Bhat)))) + bhat > 0)-1
  result <- sum((Yhat != Ytest))/length(Ytest)
  if(returnLossValue) result <- c(result, training$objective)
  return(result)
}

class_error_Bb <- function(B, b,  Atest, Ytest) {
  Bhat <- B
  bhat <- b
  Yhat <- 2*(apply(Atest, 3, function(A) sum(diag(crossprod(A, Bhat)))) + bhat > 0)-1
  return(sum((Yhat != Ytest))/length(Ytest))
}



beta_l2error <- function(A, Y, Z, B) {
  #browser()
  sqrt(sum((glmblock.fit(A, y = Y, Z = Z, intercept = TRUE)$B - B)^2))
}

#install.packages("combinat")
percentage_error_nodes <- function(Z, Zhat) {
  K = ncol(Z)
  if(K<=6) {
    require(combinat)
    1-max(sapply(permn(1:ncol(Z)), function(perm) sum(diag(crossprod(Z, Zhat[,perm]))))/nrow(Z)) 
  }else{
    cpZZ = crossprod(Z, Zhat)
    maxs = apply(cpZZ, 2, which.max)
    if(length(unique(maxs)) == ncol(Z)) {
      diag(cpZZ[maxs,]) = 0
      sum(cpZZ) / nrow(Z)
    }else{
      return(NULL)
    }
  }
}

## coclustering error
coclustering_error <- function(Z, Zhat) {
  sum(abs(tcrossprod(Z) - tcrossprod((Zhat))))/nrow(Z)^2
}
