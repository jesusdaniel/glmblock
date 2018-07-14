# Full ADMM algorithm
library(tensor)
library(glmnet)
library(Matrix)
# min l(B,b) + gamma/2*||B||_F^2 + lambda*||B||_1
#   s.t. B=ZPZ^T
glmblock <- function(A, y, K, gamma = 0, lambda = 0,
                     family = "gaussian", intercept = TRUE, Z0 = NULL,
                     fclust = "SC", MAX_STEPS = 30, rho = 1, EPS = 1e-4, rho_factor = 1, autotuning = FALSE, 
                     return_path = FALSE,
                     verbose = FALSE,
                     ...) {
  # constants
  N = ncol(A)
  n = length(y)
  #Initialize functions
  if(fclust=="SC") { fclust<- spectral_clustering} #user specified clustering
  #Initialization
  if(is.null(Z0)) {
    sigma_AY <- apply(A, c(1,2), function(a) cov(a,y)) # or do t-test for binomial
    sigma_AY <- ifelse(is.na(sigma_AY),0,sigma_AY)
    Z0 <- fclust(sigma_AY,K)              # specify clustering function
  }
  Z <- Z0
  Ar <- 2*t(apply(A, 3, uptri))
  
  # Initialize variables
  V = matrix(0, nrow = N, ncol = N)
  init = glmblock.fit(A = A, y = y, Z = Z0,
                      gamma = gamma, lambda = lambda, family = family, intercept = intercept)
  W = init$B
  b = init$b
  B = V
  converge = FALSE
  steps <- 0
  primal_crit <- c();  dual_crit <- c(); obj_value <- c();  Z_crit <- c()
  if(return_path){    Zpath = list()
  }else{    Zpath = NULL  }
  if(autotuning)
    rho <- tune_rho(V,W,b, Ar, y, intercept, K, family, gamma, lambda, fclust)
  rholist <- c()
  #################################################################################
  while(!converge & steps < MAX_STEPS) {
    Bnew_step = loss_step_offset(V, W, rho, Ar, y, intercept, family, 
                                 gamma, lambda)
    Bnew = Bnew_step$Bhat
    b = Bnew_step$b
    W_Z = clustering_step(Bnew, V, rho, K, fclust)
    Wnew = W_Z$W;    Znew <- W_Z$Z
    Vnew = V_step(V, Bnew, Wnew, rho)
    #Stopping criteria
    # 1) Primal
    crit_1 <- sqrt(sum((Bnew-Wnew)^2)) / N
    # 2) Dual
    crit_2 <- rho*sqrt(sum((Wnew - W)^2)) / N
    # 3) Change in communities
    crit_3 <- sum(abs(tcrossprod(Znew)-tcrossprod(Z))) / N^2
    t4 <- proc.time()[3]
    converge <- (crit_1 < EPS & crit_2 < EPS)
    if(return_path) Zpath[[steps + 1]] = Matrix(Znew)
    #update values
    V = Vnew; B = Bnew; Z = Znew; W = Wnew;
    primal_crit <- c(primal_crit, crit_1); dual_crit <- c(dual_crit, crit_2); Z_crit <- c(Z_crit, crit_3)
    steps <- steps + 1
    if(verbose) cat(paste("-- Rho:", rho, " Iteration:", steps, "Convergence criteria:", crit_1, crit_2, "\n"))
    ## Rho updating rules
    if(crit_1 < 10*crit_2) {  rho = rho*rho_factor
    }else{  if(crit_2 < 10*crit_1)  rho = rho/rho_factor  }
    rholist = c(rholist, rho)
  }
  result <- list(Bhat = B, b = b, What = W, Vhat = V, Zhat = Z, rholist = rholist,
                 primal = primal_crit, dual = dual_crit, primmal_dual = primal_crit/dual_crit, 
                 comm_change = Z_crit, Zpath = Zpath)
  return(result)
}

#implement parallel...
glmblock_grid <- function(A, y, K, gamma = 0, lambda = 0, family = "gaussian", intercept = TRUE, Z0 = NULL,
                          fclust = "SC", MAX_STEPS = 100, rho_grid = 1, EPS = 1e-4, 
                          return_full = FALSE, 
                          return_path = FALSE,
                          verbose = FALSE,
                          ...) {
  fit_rho <- lapply(rho_grid, function(rho) tryCatch({
    print(paste("Fitting rho = ", rho))
    glmblock(A, y, K, gamma = gamma, lambda = lambda,
             family = family,
             intercept = intercept, Z0 = Z0,
             fclust = fclust,
             MAX_STEPS = MAX_STEPS, rho = rho,
             EPS = EPS, rho_factor = 1, autotuning = FALSE,
             return_path = return_path, verbose = verbose)
    
  },
  error = function(e) {
    warning("Parameter ",rho, " terminated with an error!")
    return(NULL)
  }))
  if(family=="gaussian")
    loss_function <- sapply(fit_rho, function(glmb) {
      if(is.null(glmb))
        return(Inf)
      glmblock.fit(A,y,Z = glmb$Zhat, gamma = gamma, lambda = lambda,family = "gaussian", intercept=intercept)$objective
    })
  
  if(family=="binomial")
    loss_function <- sapply(fit_rho, function(glmb) {
      if(is.null(glmb))
        return(Inf)
      glmblock.fit(A, y, glmb$Zhat,  gamma = gamma, lambda = lambda,
                   family = "binomial", intercept=intercept)$objective
    })
  if(return_full)
    return(list(grid_obj = fit_rho, grid_loss = loss_function))
  return(fit_rho[[which.min(loss_function)]])
}


#implement parallel...
#only gaussian implemented....
cv.glmblock <- function(A, y, K_values = 1:6, gamma = 0, lambda = 0, num_folds = 5,
                             family = "gaussian", intercept = TRUE, 
                          fclust = "SC", MAX_STEPS = 100, rho_grid = 1, EPS = 1e-4, return_full = FALSE, return_path = FALSE,
                          ...) {
  #browser()
  require(cvTools)
  fold_partition = cvFolds(n = length(y), K = num_folds)
  fold_list <- lapply(1:fold_partition$K, function(j) 
    which(fold_partition$which==j))
  
  CV_LOSS <- matrix(0, nrow = length(K_values), ncol = num_folds)
  INITIAL_CV_LOSS <- matrix(0, nrow = length(K_values), ncol = num_folds)
  
  j = 1
  for(K in K_values) {
    for(i in 1:num_folds) {
      # Train test
      Atrain = A[,,-fold_list[[i]]]
      Ytrain = y[-fold_list[[i]]]
      Atest = A[,,fold_list[[i]]]
      Ytest = y[fold_list[[i]]]
      #Atrain_tensor <- array(unlist(Atrain), dim = c(N, N, length(Atrain)))
      #Atest_tensor <- array(unlist(Atest), dim = c(N, N, length(Atest)))
      
      # Initial value-----------------------------------------------------
      #Set initial value
      sigma_AY <- apply(Atrain, c(1,2),function(a) cov(Ytrain,a))
      sigma_AY <- ifelse(is.na(sigma_AY), 0, sigma_AY)
      Zhat = spectral_clustering(sigma_AY, K)
      
      # fit method--------------------------------------------------------
      fitted_Kfold <- glmblock_grid(Atrain, Ytrain, K, intercept = intercept, 
                                    gamma = gamma, lambda = lambda, 
                                    family = family,  Z0 = Zhat,
                                    MAX_STEPS = MAX_STEPS, rho_grid = rho_grid, 
                                    EPS = EPS, return_full = FALSE, return_path = FALSE)
      ### Evaluate initial test error
      INITIAL_CV_LOSS[j,i] <- MSE(Atrain, Ytrain, Zhat, 
                                  Atest, Ytest)
      ## Evaluate test error
      CV_LOSS[j,i] <- MSE(Atrain, Ytrain, 
                          fitted_Kfold$Zhat, Atest, Ytest)
      
    }
    print(K)
    j = j+1
  }
  
  Kmin = K_values[which.min(rowMeans(CV_LOSS))]
  #Set initial value
  #Atensor <- array(unlist(A), dim = c(N, N, length(A)))
  sigma_AY <- apply(A, c(1,2),function(a) cov(y,a))
  sigma_AY <- ifelse(is.na(sigma_AY), 0, sigma_AY)
  Zhat = spectral_clustering(sigma_AY, Kmin)
  
  best_K <- glmblock_grid(A, y, Kmin, intercept = intercept, 
                          gamma = gamma, lambda = lambda, 
                          family = family, Z0 = Zhat,
                                MAX_STEPS = MAX_STEPS, rho_grid = rho_grid, 
                                EPS = EPS, return_full = FALSE, return_path = FALSE)
  return(list(fit_bestK = best_K, fit_initial_bestK = Zhat, Kmin = Kmin, 
              INITIAL_CV_LOSS = INITIAL_CV_LOSS, CV_LOSS = CV_LOSS))
}


#not checked
tune_rho <- function(V,W,b, Ar, y, intercept, K, family, gamma, lambda, fclust) {
  candidate_rho = 10^seq(-2,3,by = 0.1)
  crit_ratio = c()
  for(rho in candidate_rho) {
    Bnew_step = loss_step_offset(V, W, rho, Ar, y, intercept, family, gamma, lambda)
    Bnew = Bnew_step$Bhat
    b = Bnew_step$b
    W_Z = clustering_step(Bnew, V, rho, K, fclust)
    Wnew = W_Z$W;    Znew <- W_Z$Z
    Vnew = V_step(V, Bnew, Wnew, rho)
    #Stopping criteria
    # 1) Primal
    crit_1 <- sqrt(sum((Bnew-Wnew)^2)) / ncol(W)
    # 2) Dual
    crit_2 <- rho*sqrt(sum((Wnew - W)^2)) / ncol(W)
    crit_ratio = c(crit_ratio, crit_1/crit_2)
  }
  return(candidate_rho[which.min(abs(1-crit_ratio))])
}

#lambda > 0 only implemented for logistic
loss_step_offset <- function(V, W, rho, Ar, y, intercept, family, gamma = 0, lambda = 0) {
  if(family == "binomial" & lambda > 0){
    res <- loss_step_offset_logistic(V, W, rho, Ar, y, intercept, family, gamma, lambda)
    return(res)
  }
  Vr = uptri(V)
  Wr = uptri(W)
  U <- as.vector((rho / (gamma + rho)) * (Wr - (1/rho)*Vr))
  offset = c(Ar %*% U)  
  if(lambda == 0) {
    bhat = glm(y~1, family = family, offset = offset)$coefficients
    offset = offset + bhat
    glmfit = glmnet(x = Ar, y = y, alpha = 0, 
                    lambda = 2*(rho + gamma),  
                    intercept = intercept, 
                    family = family, offset = offset, 
                    standardize = FALSE)
    b = glmfit$a0  + bhat
    Uhat = to_matrix(glmfit$beta)
    Bhat = Uhat + (rho/(gamma + rho)) *  W - (1/(gamma + rho))*V
    return(list(Bhat = Bhat, b = b))
  }
  glmfit1 = gaussian_ridge(X = Ar, Y = y, R = U, 
                           gamma = 2*(rho + gamma), lambda = 2*lambda, 
                           MAX_ITER = 50, CONV_CRIT = 1e-5)
  Bhat <- get_matrix(glmfit1$best_beta)
  b <- glmfit1$best_b
  return(list(Bhat = Bhat, b = b))
}

loss_step_offset_logistic <- function(V, W, rho, Ar, y, intercept, family, gamma = 0, lambda = 0) {
  y01 = sapply(y, function(y_i) ifelse(y_i==max(y),1, 0))
  Vr = uptri(V)
  Wr = uptri(W)
  U = rho / (gamma + rho) * (Wr - (1/rho) * Vr)
  offset = c(Ar%*%as.vector(U))
  bhat = glm(y01~1, family = family, offset = offset)$coefficients
  offset = offset + bhat
  
  glmfit = logistic_ridge(X = Ar, Y = y, 
                          gamma = 2*(rho + gamma), 
                          offset = offset, lambda = 2*lambda, 
                          R = U)
  b = glmfit$b  + bhat
  Uhat = to_matrix(glmfit$x)
  Bhat = Uhat + (rho/(gamma + rho)) *  W - (1/(gamma + rho))*V
  return(list(Bhat = Bhat, b = b))
}

clustering_step <- function(B, V, rho, K, fclust) {
  BV <- B + (1/rho)*V
  Z <- fclust(BV, K)
  P <- solve(crossprod(Z)) %*% (crossprod(Z, BV %*% Z) %*% solve(crossprod(Z)))
  W <- Z %*% (P %*% t(Z))
  return(list(W = W, Z = Z))
}

V_step <- function(V, B, W, rho) {
  V <- V + rho*(B-W)
  return(V)
}



# fits P and b in the model
# l_Z(P,b) + lambda* ||weights*P||_1 + gamma/2 * ||weights*P||_F^2
glmblock.fit <- function(A, y, Z, gamma = 0, lambda = 0,
                         family = "gaussian", intercept = T) {

  Z = Matrix(Z)
  N = ncol(A)
  ZAZ = apply(A, 3, function(a) Matrix::crossprod(Z, Matrix::crossprod(a,Z)))
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
    responses = predict(glmfit, newx = ZAZr, type = "class")
    loss = -sum((y==1)*responses-log(1+exp(responses)))/(length(y))
    Class_hat = max(y)*(responses >=0) + min(y)*(responses < 0)
    Class_error = sum(Class_hat != y) / length(y)
  }
  weights <- tcrossprod(diag(crossprod(Z))) - crossprod(Z)
  penalty <- lambda * sum(abs(weights*Phat)) + gamma/2 * sum(weights*(Phat)^2)
  objective <- loss + penalty
  #diag(Phat) <- diag(Phat) / 2
  return(list(B= Z%*%(Phat%*%t(Z)), b = b, objective = objective, Class_error = Class_error))
}


