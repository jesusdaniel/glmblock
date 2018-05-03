####
# Solves the problem
# min l(U + R) +  gamma/2 |U||_F^2  +  lambda ||U+R||_1
my_logistic_ridge <- function(X, Y, offset, gamma, lambda, R, beta_start = NULL, b_start = NA,
                               MAX_ITER = 100, verbose = FALSE,
                               CONV_CRIT = 1e-04, MAX_TIME = Inf,
                               Xtest = NULL, Ytest= NULL) {
  #browser()
  if(min(Y)==0) # remove this later
    Y = 2*Y-1
  n = length(Y); p = dim(X)[2]
  Xbar <- colMeans(X)
  Xsds <- apply(X,2, sd)
  Xsdsd_inv <- 1/Xsds
  X <- scale(X)
  # Define functions wrt X,Y,D---------------------------------------------------------------
  # B derivative, B hessian------------------------------------------------------------------
  #Derivative
  b_derivative = function(Xbeta,b)
    sum(-Y/(1+exp(Y*(Xbeta+b + offset))))/n
  #Hessian
  b_hessian = function(Xbeta,b)
    sum(1/(exp(-Y*(Xbeta+b + offset))+exp(Y*(Xbeta+b + offset))+2))/n
  
  # Beta derivative -------------------------------------------------------------------------
  grad_f <- function(Xbeta,b, beta)
    -crossprod(X,Y/(1+exp(Y*(Xbeta+b + offset))))/n + gamma*(Xsdsd_inv^2*beta)
  # F evaluation-----------------------------------------------------------------------------
  f = function(Xbeta,b, beta)
    sum(log(1+exp(-Y*(Xbeta+b + offset))))/n + gamma/2*crossprod(Xsdsd_inv*beta)
  penalty = function(beta,b)
    lambda*sum(abs(Xsdsd_inv*beta + R))
  
  #Proximal and b step-----------------------------------------------------------------------
  proximal <- function(u,alpha) {
    if(alpha>0){
      #return(list(x = Xsds*sign(R + Xsdsd_inv*u)*pmax(abs(R+Xsdsd_inv*u)- alpha*lambda,0) - Xsds*R))
      XsdsR = Xsds*R
      return(list(x = sign(XsdsR + u)*pmax(abs(XsdsR + u)- alpha*lambda*Xsdsd_inv,0) - XsdsR))
    }else{
      return(list(x = u))
    }
  }
  b_step <- function(Xbeta, b_start = 0) {
    TOL_B = 1e-04; MAX_S_B = 100
    b_n = b_start
    i = 0
    b_deriv = Inf
    while(abs(b_deriv)>TOL_B & i < MAX_S_B) {
      b_deriv = b_derivative(Xbeta,b_n)
      b_n = b_n - b_deriv/(b_hessian(Xbeta,b_n)+1*(abs(b_deriv/b_hessian(Xbeta,b_n))>100))
      if(abs(b_n)>1000){
        warning("The intercept is too big: ",b_n,"\n  Derivative = ",b_deriv,
                "\n  Hessian = ", b_hessian(Xbeta,b_n))
      }
      i = i+1
    }
    return(b_n)
  }
  
  if(is.null(beta_start))  beta_start = rep(0,p)
  if(is.na(b_start))  b_start = 0
  optimal = fista_line_search_ridge(proximal,b_step,f,grad_f,penalty,beta_start,b_start,X,Y,verbose,
                                    MAX_STEPS = MAX_ITER, TOLERANCE = CONV_CRIT,
                                    Xtest = Xtest, Ytest = Ytest,offset =  offset, Xsds = Xsdsd_inv, gamma = gamma,
                                    lambda = lambda, R = R)
  optimal$x <- optimal$x * Xsdsd_inv
  optimal$best_beta <- optimal$best_beta * Xsdsd_inv
  optimal$b <- as.double(optimal$b - crossprod(Xbar, optimal$x))
  optimal$best_b <- as.double(optimal$best_b - crossprod(Xbar, optimal$best_beta))
  return(optimal)
}

norm_s = function(vec) sqrt(sum(vec^2))
# STATUS
fista_line_search_ridge <- function(proximal_f,b_step,f,grad_f,penalty, x_start, b_start,
                                    X, Y, verbose = F,
                                    MAX_STEPS = 300, TOLERANCE = 1e-06,
                                    MAX_TIME = 10800,
                                    Xtest = NULL, Ytest = NULL, offset, Xsdsd_inv, gamma, lambda, R) {
  crits = c();  fs = c();  tks = c()
  crits_f = c()
  ID = "NULL"
  # Initialize -------------------------------------------------------------------
  xk1 = x_start;  xk = x_start
  criterion = Inf
  crit_f = Inf
  k = 0;  tk = 0.125;  beta_step = 1/2
  Xbeta = X%*%x_start
  b = b_step(Xbeta,b_start)
  best_beta = x_start;  best_b = b;
  best_f = f(Xbeta,b,x_start) + penalty(x_start,b);  best_prox = NULL
  newf = best_f
  if(verbose) cat(sprintf("---%s - Iteration %d. Criterion = -. Crit_f = -.  F = %.4f\n",ID,0,newf))
  time_start = proc.time()[1]
  crit_f_1  = crit_f
  beta_path  = list()
  beta_path[[1]] = xk1
  is_best_end = FALSE
  while(k<=5 | ((crit_f > TOLERANCE  & criterion > TOLERANCE) & (k < MAX_STEPS) &
                proc.time()[1]-time_start < MAX_TIME)){
    crit_f = newf
    is_best_end = FALSE
    k = k+1
    xk_1 = xk;    xk = xk1
    y = xk + (k-2)/(k+1) * (xk-xk_1)
    ######## line search ###################
    repeat {
      Xy = X%*%y
      z=y-tk*grad_f(Xy,b,y)
      prox = proximal_f(z,tk)
      z = prox$x
      Xbeta = X%*%z
      if(f(Xbeta,b,z) <= as.double(f(Xy,b,y)) +t(grad_f(Xy,b,y))%*%(z-y) + 1/(2*tk) * norm_s(z-y)^2)
        break #tk = 1/beta_step*tk
      #cat("  ...trying new tk\tt =",proc.time()[1]-time_start, "\n")
      tk = beta_step*tk
    }
    
    ######## line search ###################
    b = b_step(Xbeta,b)
    tks = c(tks,tk)
    xk1 = z;    criterion = norm_s(xk1-xk)
    crits = c(crits,criterion)
    newf = f(Xbeta,b, z) + penalty(xk1,b)
    
    if(crit_f>0) #crit_f = newf above
      crit_f = abs(newf - crit_f)/crit_f
    if(newf < best_f) {
      is_best_end = TRUE
      best_beta = xk1;     best_b = b;
      best_f = newf
    }
    crits_f = c(crits_f, crit_f)
    fs = c(fs, newf)
    if(verbose)  cat(sprintf("--%s - Iter %d. Sparse = %d.  Crit_f = %.4f. tk = %.2f. F = %.4f. t=%2f. N=%d. ER=%.1f CE=%.2f\n",
                             ID,k,sum(xk1==0),crit_f,tk,newf,proc.time()[1]-time_start,NA, NA,
                             classification_error(Xtest,Ytest,xk1,b)))
    #if(is.na(crit_f) | is.null(crits_f) | crit_f==0) browser()
    if(crit_f==0 | abs(crit_f-crit_f_1)/crit_f<0.1)
      tk = 2*tk
    crit_f_1 = crit_f
    beta_path[[k+1]] = xk1
  }
  if(crit_f < TOLERANCE | criterion < TOLERANCE) {
    status = 1
  }else if(k>=MAX_STEPS){
    status = 2
  }else{
    status = 3
  }
  return(list(x = xk1, b=b, crits = crits_f, fs = fs, L = tk, tks = tks, conv = criterion,
              best_f = best_f, best_beta = best_beta, best_b = best_b,
              is_best_end = is_best_end,
              beta_path = beta_path,
              status = status))
}
