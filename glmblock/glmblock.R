# Full ADMM algorithm
# min l(B,b) + gamma/2*||B||_F^2 + lambda*||B||_1
#   s.t. B=ZPZ^T


#' 
#' Function to perform block penalized regression with supervised community detection
#' 
#' @param Adj_list a list with adjacency matrices of the same size n x n
#' @param y a vector of responses
#' @param family response type, either "gaussian" for continuous responses, or "binomial" for classification
#' @param intercept indicate whether intercept should be fitted
#' @param K number of communities in the block constraint. If \code{initial.communities} is given, then this parameter is ignored.
#' @param gamma ridge penalty parameter
#' @param lambda lasso penalty parameter
#' @param method.communities method for supervised community detection. The options available are "SC" for
#' spectral clustering or "ADMM" for alternating direction method of multipliers. Default is \code{"SC"}.
#' @param rho.grid ADMM penalty parameters to try in the optimization problem.
#' @param Z0 membership matrix of size n by K, with the initial value for the community membership search of ADMM.
#' @param fclust clustering function for approximating the solution of the block constraint. Currently only spectral clustering ("SC") is implemented.
#' @param MAX_STEPS maximum number of iterations of the ADMM algorithm
#' @param EPS stopping criteria od ADMM algorithm
#' @param rho.factor change of rho on each iteration of ADMM algorithm
#' @param verbose print output on each iteration of ADMM
#' 
#' @return 
#' 
#' @references 
#' 
#' @author Jes\'us Arroyo <jarroyor@umich.edu>
glmblock <- function(Adj_list, y, 
                     family = "gaussian", intercept = TRUE, # regression parameters
                     K, gamma = 0, lambda = 0,              # tuning parameters
                     method.communities = c("SC", "ADMM"),
                     rho.grid = 10^seq(-2, 3),                          # controls the performance of ADMM
                     Z0 = NULL,                             # ADMM parameters
                     fclust = "SC", MAX_STEPS = 30, 
                     EPS = 1e-4, rho.factor = 1, 
                     verbose = FALSE,
                     ...) {
  # constants
  n = ncol(Adj_list[[1]])
  m = length(y)
  
  
  #Initialize functions
  if(fclust=="SC") { fclust<- spectral_clustering} #user specified clustering
  # Methods
  method.communities <-  match.arg(method.communities, c("SC", "ADMM", "user"))
  ############################################################################
  ############################################################################
  # 1) Spectral clustering
  if(method.communities == "SC") {
    # make networks into 3-dimensional array
    A = array(0, dim = c(n, n, m))
    for(i in 1:m) {
      A[, , i] = Adj_list[[i]]
    }
    sigma_AY <- apply(A, c(1,2), function(a) cov(a,y)) # or do t-test for binomial
    sigma_AY <- ifelse(is.na(sigma_AY),0,sigma_AY)
    Z <- fclust(sigma_AY,K)              # specify clustering function
    glm.result = glmblock.fit(Adj_list, y, family, intercept, Z, gamma, lambda)
  }
  ############################################################################
  ############################################################################
  if(method.communities == "ADMM") {
    # check initial
    if(is.null(Z0)) {
      if(verbose) cat(paste(c("Calculating initialization...\n")))
      # make networks into 3-dimensional array
      A = array(0, dim = c(n, n, m))
      for(i in 1:m) {
        A[, , i] = Adj_list[[i]]
      }
      sigma_AY <- apply(A, c(1,2), function(a) cov(a,y)) # or do t-test for binomial
      sigma_AY <- ifelse(is.na(sigma_AY),0,sigma_AY)
      Z0 <- fclust(sigma_AY,K)              # specify clustering function
    } else{
      K.z0 <- ncol(Z0)
      if(K.z0 != K) {
        stop("Number of communities in Z0 and K are different.")
      }
    }
    glm_solutions_grid <- lapply(rho.grid, function(rho) tryCatch({
      if(verbose) cat(paste(c("Running rho =", rho, "...\n")))
      admm.rho = glmblock.admm(Adj_list, y, Z0, rho,
                    gamma, lambda, family, intercept,
                    EPS, MAX_STEPS, fclust, rho.factor, autotuning = FALSE, 
                    verbose = verbose, return_path = TRUE)
    },
    error = function(e) {
      warning("Parameter ",rho, " terminated with an error!")
      return(NULL)
    }))
    loss_function <- sapply(glm_solutions_grid, function(glmb) {
        if(is.null(glmb))
          return(Inf)
        # check path of solutions
        grid.losses = sapply(glmb$Zpath, function(Z) {
          glmblock.fit(Adj_list, y, family, intercept, Z, gamma, lambda)$objective
        })
        min(grid.losses)
      })
    
    if(min(loss_function) == Inf) { stop("All rho values terminated with an error.")}
    
    # search minimum
    which.rho.min = which.min(loss_function)
    rho.min.losses = sapply(glm_solutions_grid[[which.rho.min]]$Zpath, function(Z) {
      glmblock.fit(Adj_list, y, family, intercept, Z, gamma, lambda)$objective
    })
    Z = glm_solutions_grid[[which.min(loss_function)]]$Zpath[[which.min(rho.min.losses)]]
    glm.result = glmblock.fit(Adj_list, y, family, intercept, 
                              Z, gamma, lambda)
  }
  result = list(glm.result = glm.result, Z = Z)
  
  return(result)
}





