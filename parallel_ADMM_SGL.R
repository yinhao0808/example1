#####################################
### ADMM: Sparse Group Lasso 2    ###
### This is for convergence check ###
#####################################

sgl_admm <- function(X,Y,rho, group, lambda1, lambda2, b0, iter=500, family = "gaussian", e.abs=1E-3,e.rel=1E-3){
  
  
  # Lasso function #
  l1 <- function(u, lambda){
    uhat <- abs(u) - lambda
    prox <- sign(u) * pmax(rep(0, length(u)), uhat)
    return(prox)
  }
  
  # group Lasso function #
  l2 <- function(a, lambda) {
    a * pmax(0, 1 - lambda / sqrt(sum(a ^ 2)))
  }
  
  # norm function
  l2norm <- function(x) sqrt(sum(x^2))
  
  n=nrow(X); p <- ncol(X)
  
  numofgroups <- group
  groupsize <- p/numofgroups
  groupID <- split(1:(groupsize*numofgroups), rep(1:numofgroups, rep(groupsize, numofgroups)))
  
  for (i in seq(numofgroups)){
    X[,groupID[[i]]] <- sqrt(n)*orthonormalization(X[,groupID[[i]]],basis=FALSE)
  }
  
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  

  # Initialize coefficients matrix 
  beta <- matrix(0, nrow=iter, ncol=p) 
  beta[1,] <- b0#rep(1, p)
  
  alpha <- matrix(0, nrow = iter, ncol=q)
  
  # Initialize gamma and v 
  # In our proof, gamma is z and v is u
  gamma <- matrix(0, nrow=iter, ncol=p)
  v <- rep(0, p)    # Lagrangian; 
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  
  # Compute inverse matrix to obtain beta
  if (family == "gaussian"){weight_adjust = 1}
  if (family == "binomial"){weight_adjust = 4}
  invmat <- solve(XtX/(weight_adjust*n) + diag(rho, p) + diag(rho, p))
  
  # Initialize residual vectors.
  s <- 0    # dual residual
  r <- 0    # primal residual
  
  # Initialize iterator
  t <- 0
  
  if (family == "gaussian"){
  # ADMM updates
  for (t in 2:iter){
    
    # Update beta, gamma, v, and obj
    beta[t,] <- invmat %*% (XtY/n + (rho+rho) * (gamma[t-1,]-v) )
    
    for (g in seq(group)){
      g.index <- ((g-1)*(p/group)+1):(g*(p/group))
      M <- (rho + rho) * (beta[t,g.index] + v[g.index])
      den <- l2norm(gamma[t-1,g.index])
      if (den == 0){den <- 10^(-6)}
      G <- rho + lambda1/den
      h <- l1(M/(G+rho), lambda2/(G+rho))*G
      gamma[t,g.index] <- l2(h/rho, lambda1/rho)
    }
    
    v <- v + beta[t,] - gamma[t,]
    
    # Calculate residuals for iteration t
    r <- beta[t,] - gamma[t,]
    s <- -rho * (gamma[t,] - gamma[t-1,])
    
    r.norm <- l2norm(r)
    s.norm <- l2norm(s)
    
    e.primal <- sqrt(p) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
    e.dual <-  sqrt(n) * e.abs + e.rel * l2norm(v)
    
    # Check convergence
    if (r.norm <= e.primal && s.norm <= e.dual){
      # Remove excess beta
      beta <- beta[-((t+1):nrow(beta)),]
      
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      break
    }
  }
  }
  
  if (family == "binomial"){
    
  Y_origin <- Y
  # ADMM updates
  for (t in 2:iter){
    
    ww=1/4;
    P0=1/(1+exp(-E%*%alpha[t-1,]-X%*%beta[t-1,]))
    y.work = E%*%alpha[t-1,] + X%*%beta[t-1,]+ww^(-1)*(Y_origin-P0)
    alpha[t,] <- ww*solve(t(E)%*%E)%*%t(E)%*%(Y_origin-P0)
    Y=y.work - E%*%alpha[t,]
    
    XtY <- crossprod(X, Y)
    
    # Update beta, gamma, v, and obj
    beta[t,] <- invmat %*% (XtY/(4*n) + (rho+rho) * (gamma[t-1,]-v) )
    
    for (g in seq(group)){
      g.index <- ((g-1)*(p/group)+1):(g*(p/group))
      M <- (rho + rho) * (beta[t,g.index] + v[g.index])
      den <- l2norm(gamma[t-1,g.index])
      if (den == 0){den <- 10^(-6)}
      G <- rho + lambda1/den
      h <- l1(M/(G+rho), lambda2/(G+rho))*G
      gamma[t,g.index] <- l2(h/rho, lambda1/rho)
    }
    
    v <- v + beta[t,] - gamma[t,]
    
    # Calculate residuals for iteration t
    r <- beta[t,] - gamma[t,]
    s <- -rho * (gamma[t,] - gamma[t-1,])
    
    r.norm <- l2norm(r)
    s.norm <- l2norm(s)
    
    e.primal <- sqrt(p) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
    e.dual <-  sqrt(n) * e.abs + e.rel * l2norm(v)
    
    # Check convergence
    if (r.norm <= e.primal && s.norm <= e.dual){
      # Remove excess beta
      beta <- beta[-((t+1):nrow(beta)),]
      alpha <- alpha[-((t+1):nrow(alpha)),]
      
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      break
    }
  }
  }
  
  result <- list("beta.hat"=beta[nrow(beta),],"beta"=beta, "alpha.hat" = alpha[nrow(alpha),],
                 "conv"=message, "iter"=t)
  return(result)
}


sgl_admm_sample <- function(X,Y,rho, group, fold=1, lambda1, lambda2, b0, family = "gaussian", iter=500,e.abs=1E-3,e.rel=1E-3){
  
  numofgroups_sample <- fold
  
  start_time <- Sys.time()
  # Lasso function #
  l1 <- function(u, lambda){
    uhat <- abs(u) - lambda
    prox <- sign(u) * pmax(rep(0, length(u)), uhat)
    return(prox)
  }
  
  # group Lasso function #
  l2 <- function(a, lambda) {
    a * pmax(0, 1 - lambda / sqrt(sum(a ^ 2)))
  }
  
  # norm function
  l2norm <- function(x) sqrt(sum(x^2))
  
  n=nrow(X); p <- ncol(X)
  
  numofgroups <- group
  groupsize <- p/numofgroups
  groupID <- split(1:(groupsize*numofgroups), rep(1:numofgroups, rep(groupsize, numofgroups)))
  
  ### Parallel
  for (i in seq(numofgroups)){
    X[,groupID[[i]]] <- sqrt(n)*orthonormalization(X[,groupID[[i]]],basis=FALSE)
  }
  
  # Initialize coefficients matrix 
  beta <- matrix(0, nrow=iter, ncol=p) 
  beta[1,] <- b0#rep(1, p)
  
  # Initialize gamma and v 
  # In our proof, gamma is z and v is u
  gamma <- matrix(0, nrow=iter, ncol=p)
  v <- matrix(0, nrow = numofgroups_sample, p)      # Lagrangian; 
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  
  # Initialize residual vectors.
  s <- 0    # dual residual
  r <- 0    # primal residual
  
  # Initialize iterator
  t <- 0
  
  groupsize_sample <- n/numofgroups_sample
  groupID_sample <- split(1:(groupsize_sample*numofgroups_sample), rep(1:numofgroups_sample, rep(groupsize_sample, numofgroups_sample)))
  
  end_time <- Sys.time()
  time0 <- as.numeric(end_time - start_time)
  
  if (family == "gaussian"){weight_adjust = 1}
  if (family == "binomial"){weight_adjust = 4}
  
  time02 <- NULL
  invmat_group <- list()
  XtY_group <- list()
  for ( i in seq(length(groupID_sample))){
    start_time <- Sys.time()
    tmp_id <- groupID_sample[[i]]
    invmat_group[[i]] <- diag(1/rho, ncol(X[tmp_id,])) - (1/rho)*t(X[tmp_id,])%*%solve(weight_adjust*groupsize_sample*diag(rho, nrow(X[tmp_id,])) + X[tmp_id,]%*%t(X[tmp_id,]))%*%X[tmp_id,]
    XtY_group[[i]] <- t(X[tmp_id,])%*%Y[tmp_id]
    end_time <- Sys.time()
    time02 <- c(time02, as.numeric(end_time - start_time))
  }
  time02 <- max(time02)
  
  time <- rep(0,iter)
  
  if (family == "gaussian"){
  
  # ADMM updates
  for (t in 2:iter){

    ########################################
    # Parallel  ADMM steps
    time2 <- NULL
    beta_i <- list()
    for ( i in seq(length(groupID_sample))){
      start_time <- Sys.time()
      tmp_id <- groupID_sample[[i]]
      #invmat <- diag(1/rho, ncol(X[tmp_id,])) - (1/rho)*t(X[tmp_id,])%*%solve(groupsize*diag(rho, nrow(X[tmp_id,])) + X[tmp_id,]%*%t(X[tmp_id,]))%*%X[tmp_id,]  #solve(t(X[tmp_id,])%*%X[tmp_id,]+rho*diag(p))
      beta_i[[i]] <-  invmat_group[[i]]%*%(XtY_group[[i]]/(weight_adjust*groupsize_sample)+rho*(gamma[t-1,] - v[i,])) #invmat%*%(t(X[tmp_id,])%*%Y[tmp_id]/groupsize+rho*(gamma[t-1,] - v[i,]))
      end_time <- Sys.time()
      time2 <- c(time2, as.numeric(end_time - start_time))
    }
    time2 <- max(time2)
    ########################################
    
    start_time <- Sys.time()
    
    tmp_beta_i <- matrix(unlist(beta_i), ncol = p, byrow = TRUE)
    
    mean_beta_i <- apply(tmp_beta_i,2,mean)
    mean_v <- apply(v,2,mean)
    
    beta[t,] <- mean_beta_i
    
    for (g in seq(group)){
      g.index <- ((g-1)*(p/group)+1):(g*(p/group))
      M <- (rho + rho) * (mean_beta_i[g.index] + mean_v[g.index])
      den <- l2norm(gamma[t-1,g.index])
      if (den == 0){den <- 10^(-6)}
      G <- rho + lambda1/den
      h <- l1(M/(G+rho), lambda2/(G+rho))*G
      gamma[t,g.index] <- l2(h/rho, lambda1/rho)
    }
    
    v <- t(t(v + tmp_beta_i) - gamma[t,])
    
    # Calculate residuals for iteration t
    r <- beta[t,] - gamma[t,]
    s <- -rho * (gamma[t,] - gamma[t-1,])
    
    r.norm <- l2norm(r)
    s.norm <- l2norm(s)
    
    # e.primal <- sqrt(p) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
    # e.dual <-  sqrt(groupsize_sample) * e.abs + e.rel * l2norm(v)
    
    end_time <- Sys.time()
    time3 <- as.numeric(end_time - start_time)
    
    time[t] <- time2+time3
    # Check convergence
    if (r.norm <= e.rel && s.norm <= e.abs){
      # Remove excess beta
      beta <- beta[-((t+1):nrow(beta)),]
      
      time <- time[-((t+1):length(time))]
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      break
    }
  }
  }
  
  if (family == "binomial"){
    
    Y_origin <- Y
  
  for (t in 2:iter){
    
    ww=1/4;
    P0=1/(1+exp(-X%*%beta[t-1,]))
    y.work=X%*%beta[t-1,]+ww^(-1)*(Y_origin-P0)
    Y=y.work
    
    ########################################
    # Parallel  ADMM steps
    time2 <- NULL
    beta_i <- list()
    for ( i in seq(length(groupID_sample))){
      start_time <- Sys.time()
      tmp_id <- groupID_sample[[i]]
      #invmat <- diag(1/rho, ncol(X[tmp_id,])) - (1/rho)*t(X[tmp_id,])%*%solve(groupsize*diag(rho, nrow(X[tmp_id,])) + X[tmp_id,]%*%t(X[tmp_id,]))%*%X[tmp_id,]  #solve(t(X[tmp_id,])%*%X[tmp_id,]+rho*diag(p))
      beta_i[[i]] <-  invmat_group[[i]]%*%(XtY_group[[i]]/(weight_adjust*groupsize_sample)+rho*(gamma[t-1,] - v[i,])) #invmat%*%(t(X[tmp_id,])%*%Y[tmp_id]/groupsize+rho*(gamma[t-1,] - v[i,]))
      end_time <- Sys.time()
      time2 <- c(time2, as.numeric(end_time - start_time))
    }
    time2 <- max(time2)
    ########################################
    
    start_time <- Sys.time()
    
    tmp_beta_i <- matrix(unlist(beta_i), ncol = p, byrow = TRUE)
    
    mean_beta_i <- apply(tmp_beta_i,2,mean)
    mean_v <- apply(v,2,mean)
    
    beta[t,] <- mean_beta_i
    
    for (g in seq(group)){
      g.index <- ((g-1)*(p/group)+1):(g*(p/group))
      M <- (rho + rho) * (mean_beta_i[g.index] + mean_v[g.index])
      den <- l2norm(gamma[t-1,g.index])
      if (den == 0){den <- 10^(-6)}
      G <- rho + lambda1/den
      h <- l1(M/(G+rho), lambda2/(G+rho))*G
      gamma[t,g.index] <- l2(h/rho, lambda1/rho)
    }
    
    v <- t(t(v + tmp_beta_i) - gamma[t,])
    
    # Calculate residuals for iteration t
    r <- beta[t,] - gamma[t,]
    s <- -rho * (gamma[t,] - gamma[t-1,])
    
    r.norm <- l2norm(r)
    s.norm <- l2norm(s)
    
    # e.primal <- sqrt(p) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
    # e.dual <-  sqrt(groupsize_sample) * e.abs + e.rel * l2norm(v)
    
    end_time <- Sys.time()
    time3 <- as.numeric(end_time - start_time)
    
    time[t] <- time2+time3
    # Check convergence
    if (r.norm <= e.rel && s.norm <= e.abs){
      # Remove excess beta
      beta <- beta[-((t+1):nrow(beta)),]
      
      time <- time[-((t+1):length(time))]
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      break
    }
  }
  }
  
  
  result <- list("beta.hat"=beta[nrow(beta),],"beta"=beta,
                 "conv"=message, "iter"=t, "time" = time, "time0" = time02)
  return(result)
}



sgl_admm_par <- function(X, Y, rho, b0=b0, lambda1, lambda2, numofgroups, iter=50000,e.abs=1E-3,e.rel=1E-6){
  
  start_time <- Sys.time()
  # Lasso function #
  l1 <- function(u, lambda){
    uhat <- abs(u) - lambda
    prox <- sign(u) * pmax(rep(0, length(u)), uhat)
    return(prox)
  }
  
  # group Lasso function #
  l2 <- function(a, lambda) {
    a * pmax(0, 1 - lambda / sqrt(sum(a ^ 2)))
  }
  
  
  # norm function
  l2norm <- function(x) sqrt(sum(x^2))
  
  n=nrow(X); p <- ncol(X)
  
  groupsize <- p/numofgroups
  groupID <- split(1:(groupsize*numofgroups), rep(1:numofgroups, rep(groupsize, numofgroups)))
  
  for (i in seq(numofgroups)){
    X[,groupID[[i]]] <- orthonormalization(X[,groupID[[i]]],basis=FALSE)
  }
  
  # Initialize coefficients matrix 
  beta <- matrix(0, nrow=iter, ncol=p) 
  beta[1,] <- b0#rep(0, p)
  
  # Initialize gamma and v 
  # In our proof, gamma is z and v is u
  gamma <- matrix(0, nrow=iter, ncol=n)
  v <- rep(0, n)    # Lagrangian; 
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  
  # Initialize residual vectors.
  s <- 0    # dual residual
  r <- 0    # primal residual
  
  # Initialize iterator
  t <- 0
  
  end_time <- Sys.time()
  time0 <- as.numeric(end_time - start_time)
  
  time <- rep(0,iter)
  # ADMM updates
  for (t in 2:iter){
    
    start_time <- Sys.time()
    X_i_beta_i_mean <- 1/numofgroups*X%*%beta[t-1,]
    end_time <- Sys.time()
    time1 <- as.numeric(end_time - start_time)
    
    #################################################
    # Parallel ADMM steps
    beta_i <- NULL
    time2 <- NULL
    for (i in seq(numofgroups)){
      start_time <- Sys.time()
      M <- X[,groupID[[i]]]%*%beta[t-1,groupID[[i]]] + gamma[t-1,] - X_i_beta_i_mean - v
      den <- l2norm(beta[t-1,groupID[[i]]])
      if (den == 0){den=10^(-6)}
      g <- 1+(lambda2/rho)/den
      beta_i_tmp <- l1(1/g*t(X[,groupID[[i]]])%*%M,lambda1/(rho*g))
      beta_i <- rbind(beta_i,  l2(beta_i_tmp*g, lambda2/rho))# 1/n
      end_time <- Sys.time()
      time2 <- c(time2, as.numeric(end_time - start_time))
    }
    beta[t,] <- beta_i
    time2 <- max(time2)
    #################################################
    
    start_time <- Sys.time()
    X_i_beta_i_mean_new <- 1/numofgroups*X%*%beta[t,]
    gamma[t,] <- 1/(numofgroups+rho)*(rho*X_i_beta_i_mean_new + rho*v + Y)
    
    v <- v + X_i_beta_i_mean_new - gamma[t,]
    
    # Calculate residuals for iteration t
    r <- X_i_beta_i_mean_new - gamma[t,]
    s <- -rho * (gamma[t,] - gamma[t-1,])
    
    r.norm <- l2norm(r)
    s.norm <- l2norm(s)
    
    e.primal <- sqrt(n) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
    e.dual <-  sqrt(n) * e.abs + e.rel * l2norm(v)
    
    end_time <- Sys.time()
    time3 <- as.numeric(end_time - start_time)
    
    time[t] <- time1+time2+time3
    
    # Check convergence
    if (r.norm <= e.primal && s.norm <= e.dual){
      # Remove excess beta
      beta <- beta[-((t+1):nrow(beta)),]
      
      time <- time[-((t+1):length(time))] #+ time0
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      break
    }
  }
  result <- list("beta.hat"=beta[nrow(beta),],"beta"=beta,
                 "conv"=message, "iter"=t, "time" = time, "time0" = time0)
  return(result)
}


# sgl_admm_par_v2 <- function(X, Y, rho, b0=b0, lambda1, lambda2, numofgroups, iter=50000,e.abs=1E-3,e.rel=1E-6){
#   
#   start_time <- Sys.time()
#   # Lasso function #
#   l1 <- function(u, lambda){
#     uhat <- abs(u) - lambda
#     prox <- sign(u) * pmax(rep(0, length(u)), uhat)
#     return(prox)
#   }
#   
#   # group Lasso function #
#   l2 <- function(a, lambda) {
#     a * pmax(0, 1 - lambda / sqrt(sum(a ^ 2)))
#   }
#   
#   
#   # norm function
#   l2norm <- function(x) sqrt(sum(x^2))
#   
#   n=nrow(X); p <- ncol(X)
#   
#   groupsize <- p/numofgroups
#   groupID <- split(1:(groupsize*numofgroups), rep(1:numofgroups, rep(groupsize, numofgroups)))
#   
#   for (i in seq(numofgroups)){
#     X[,groupID[[i]]] <- orthonormalization(X[,groupID[[i]]],basis=FALSE)
#   }
#   
#   # Initialize coefficients matrix 
#   beta <- matrix(0, nrow=iter, ncol=p) 
#   beta[1,] <- b0#rep(0, p)
#   
#   # Initialize gamma and v 
#   # In our proof, gamma is z and v is u
#   gamma <- matrix(0, nrow=iter, ncol=n)
#   v <- rep(0, n)    # Lagrangian; 
#   
#   # Initialize convergence message in case convergence not reached
#   message <- "Convergence not reached..."
#   
#   # Initialize residual vectors.
#   s <- 0    # dual residual
#   r <- 0    # primal residual
#   
#   # Initialize iterator
#   t <- 0
#   
#   end_time <- Sys.time()
#   time0 <- as.numeric(end_time - start_time)
#   
#   time <- rep(0,iter)
#   # ADMM updates
#   for (t in 2:iter){
#     
#     start_time <- Sys.time()
#     X_i_beta_i_mean <- 1/numofgroups*X%*%beta[t-1,]
#     end_time <- Sys.time()
#     time1 <- as.numeric(end_time - start_time)
#     
#     #################################################
#     # Parallel ADMM steps
#     beta_i <- NULL
#     time2 <- NULL
#     for (i in seq(numofgroups)){
#       start_time <- Sys.time()
#       M <- X[,groupID[[i]]]%*%beta[t-1,groupID[[i]]] + gamma[t-1,] - X_i_beta_i_mean - v
#       den <- l2norm(beta[t-1,groupID[[i]]])
#       if (den == 0){den=10^(-6)}
#       g <- 1+(lambda2/rho)/den
#       beta_i_tmp <- l1(1/g*t(X[,groupID[[i]]])%*%M,lambda1/(rho*g))
#       beta_i <- rbind(beta_i,  l2(beta_i_tmp*g, lambda2/rho))# 1/n
#       end_time <- Sys.time()
#       time2 <- c(time2, as.numeric(end_time - start_time))
#     }
#     beta[t,] <- beta_i
#     time2 <- max(time2)
#     #################################################
#     
#     start_time <- Sys.time()
#     X_i_beta_i_mean_new <- 1/numofgroups*X%*%beta[t,]
#     gamma[t,] <- 1/(numofgroups+rho)*(rho*X_i_beta_i_mean_new + rho*v + Y)
#     
#     v <- v + X_i_beta_i_mean_new - gamma[t,]
#     
#     # Calculate residuals for iteration t
#     r <- X_i_beta_i_mean_new - gamma[t,]
#     s <- -rho * (gamma[t,] - gamma[t-1,])
#     
#     r.norm <- l2norm(r)
#     s.norm <- l2norm(s)
#     
#     e.primal <- sqrt(n) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
#     e.dual <-  sqrt(n) * e.abs + e.rel * l2norm(v)
#     
#     end_time <- Sys.time()
#     time3 <- as.numeric(end_time - start_time)
#     
#     time[t] <- time1+time2+time3
#     
#     # Check convergence
#     if (r.norm <= e.primal && s.norm <= e.dual){
#       # Remove excess beta
#       beta <- beta[-((t+1):nrow(beta)),]
#       
#       time <- time[-((t+1):length(time))] #+ time0
#       # Update convergence message
#       message <- sprintf("Convergence reached after %i iterations", (t))
#       break
#     }
#   }
#   result <- list("beta.hat"=beta[nrow(beta),],"beta"=beta,
#                  "conv"=message, "iter"=t, "time" = time, "time0" = time0)
#   return(result)
# }



sgl_admm_par_v2 <- function(X, Y, rho, b0=b0, lambda1, lambda2, group, fold=1, family = "gaussian", iter=500,e.abs=1E-3,e.rel=1E-3){
  
  start_time <- Sys.time()
  # Lasso function #
  l1 <- function(u, lambda){
    uhat <- abs(u) - lambda
    prox <- sign(u) * pmax(rep(0, length(u)), uhat)
    return(prox)
  }
  
  # group Lasso function #
  l2 <- function(a, lambda) {
    a * pmax(0, 1 - lambda / sqrt(sum(a ^ 2)))
  }
  
  
  # norm function
  l2norm <- function(x) sqrt(sum(x^2))
  
  n=nrow(X); p <- ncol(X)
  
  groupsize <- p/fold
  group_fold <- group/fold
  groupsize_fold <- groupsize/group_fold
  
  groupID <- split(1:(groupsize*fold), rep(1:fold, rep(groupsize, fold)))
  
  for (i in seq(fold)){
    X_temp <- X[,groupID[[i]]]
    for (ii in seq(group_fold) ){
      temp_id <- ((ii-1)*groupsize_fold+1):(ii*groupsize_fold)
      X_temp[, temp_id] <- orthonormalization(X_temp[, temp_id], basis=FALSE)
    }
    X[,groupID[[i]]] <- X_temp
  }
  
  # Initialize coefficients matrix 
  beta <- matrix(0, nrow=iter, ncol=p) 
  beta[1,] <- b0#rep(0, p)
  
  alpha <- matrix(0, nrow = iter, ncol=q)
  
  # Initialize gamma and v 
  # In our proof, gamma is z and v is u
  gamma <- matrix(0, nrow=iter, ncol=n)
  v <- rep(0, n)    # Lagrangian; 
  
  # Initialize convergence message in case convergence not reached
  message <- "Convergence not reached..."
  
  # Initialize residual vectors.
  s <- 0    # dual residual
  r <- 0    # primal residual
  
  # Initialize iterator
  t <- 0
  
  end_time <- Sys.time()
  time0 <- as.numeric(end_time - start_time)
  
  time <- rep(0,iter)
  
  if (family == "gaussian"){
  # ADMM updates
  for (t in 2:iter){
    
    start_time <- Sys.time()
    X_i_beta_i_mean <- 1/fold*X%*%beta[t-1,]
    end_time <- Sys.time()
    time1 <- as.numeric(end_time - start_time)
    
    #################################################
    # Parallel ADMM steps
    beta_i <- NULL
    time2 <- NULL
    for (i in seq(fold)){
      start_time <- Sys.time()
      X_temp = X[,groupID[[i]]]
      beta_temp = beta[t-1, groupID[[i]]]
      
      beta_i_j <- NULL
      for(jj in seq(group_fold)){
        j.index = ((jj-1)*groupsize_fold+1):(jj*groupsize_fold)
        M = X_temp[,j.index]%*%beta_temp[j.index] + gamma[t-1,] - X_i_beta_i_mean - v
        den <- l2norm(beta_temp[j.index])
        if (den == 0){den=10^(-6)}
        g <- 1+(lambda2/rho)/den
        beta_i_j_tmp <- l1(1/g*t(X_temp[,j.index])%*%M,lambda1/(rho*g))
        beta_i_j <- rbind(beta_i_j,  l2(beta_i_j_tmp*g, lambda2/rho))
      }
      
      beta_i <- rbind(beta_i, beta_i_j)
      end_time <- Sys.time()
      time2 <- c(time2, as.numeric(end_time - start_time))
    }
    beta[t,] <- beta_i
    time2 <- max(time2)
    #################################################
    
    start_time <- Sys.time()
    X_i_beta_i_mean_new <- 1/fold*X%*%beta[t,]
    gamma[t,] <- 1/(fold+rho)*(rho*X_i_beta_i_mean_new + rho*v + Y)
    
    v <- v + X_i_beta_i_mean_new - gamma[t,]
    
    # Calculate residuals for iteration t
    r <- X_i_beta_i_mean_new - gamma[t,]
    s <- -rho * (gamma[t,] - gamma[t-1,])
    
    r.norm <- l2norm(r)
    s.norm <- l2norm(s)
    
    e.primal <- sqrt(n) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
    e.dual <-  sqrt(n) * e.abs + e.rel * l2norm(v)
    
    end_time <- Sys.time()
    time3 <- as.numeric(end_time - start_time)
    
    time[t] <- time1+time2+time3
    
    # Check convergence
    if (r.norm <= e.primal && s.norm <= e.dual){
      # Remove excess beta
      beta <- beta[-((t+1):nrow(beta)),]
      
      time <- time[-((t+1):length(time))] #+ time0
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      break
    }
  }
  }
  
  if (family == "binomial"){
    
  Y_origin <- Y
  # ADMM updates
  for (t in 2:iter){
    
    ww=1/4;
    P0=1/(1+exp(-E%*%alpha[t-1,]-X%*%beta[t-1,]))
    y.work = E%*%alpha[t-1,] + X%*%beta[t-1,]+ww^(-1)*(Y_origin-P0)
    alpha[t,] <- ww*solve(t(E)%*%E)%*%t(E)%*%(Y_origin-P0)
    Y=y.work - E%*%alpha[t,]
    
    start_time <- Sys.time()
    X_i_beta_i_mean <- 1/fold*X%*%beta[t-1,]
    end_time <- Sys.time()
    time1 <- as.numeric(end_time - start_time)
    
    #################################################
    # Parallel ADMM steps
    beta_i <- NULL
    time2 <- NULL
    for (i in seq(fold)){
      start_time <- Sys.time()
      X_temp = X[,groupID[[i]]]
      beta_temp = beta[t-1, groupID[[i]]]
      
      beta_i_j <- NULL
      for(jj in seq(group_fold)){
        j.index = ((jj-1)*groupsize_fold+1):(jj*groupsize_fold)
        M = X_temp[,j.index]%*%beta_temp[j.index] + gamma[t-1,] - X_i_beta_i_mean - v
        den <- l2norm(beta_temp[j.index])
        if (den == 0){den=10^(-6)}
        g <- 1+(lambda2/rho)/den
        beta_i_j_tmp <- l1(1/g*t(X_temp[,j.index])%*%M,lambda1/(rho*g))
        beta_i_j <- rbind(beta_i_j,  l2(beta_i_j_tmp*g, lambda2/rho))
      }
      
      beta_i <- rbind(beta_i, beta_i_j)
      end_time <- Sys.time()
      time2 <- c(time2, as.numeric(end_time - start_time))
    }
    beta[t,] <- beta_i
    time2 <- max(time2)
    #################################################
    
    start_time <- Sys.time()
    X_i_beta_i_mean_new <- 1/fold*X%*%beta[t,]
    gamma[t,] <- 1/(fold+rho)*(rho*X_i_beta_i_mean_new + rho*v + Y)
    
    v <- v + X_i_beta_i_mean_new - gamma[t,]
    
    # Calculate residuals for iteration t
    r <- X_i_beta_i_mean_new - gamma[t,]
    s <- -rho * (gamma[t,] - gamma[t-1,])
    
    r.norm <- l2norm(r)
    s.norm <- l2norm(s)
    
    e.primal <- sqrt(n) * e.abs + e.rel * max(l2norm(beta[t,]), l2norm(gamma[t,])) 
    e.dual <-  sqrt(n) * e.abs + e.rel * l2norm(v)
    
    end_time <- Sys.time()
    time3 <- as.numeric(end_time - start_time)
    
    time[t] <- time1+time2+time3
    
    print(t)
    # Check convergence
    if (r.norm <= e.primal && s.norm <= e.dual){
      # Remove excess beta
      beta <- beta[-((t+1):nrow(beta)),]
      alpha <- alpha[-((t+1):nrow(alpha)),]
      
      time <- time[-((t+1):length(time))] #+ time0
      # Update convergence message
      message <- sprintf("Convergence reached after %i iterations", (t))
      break
    }
  }
  }
  
  
  result <- list("beta.hat"=beta[nrow(beta),],"beta"=beta, "alpha.hat" = alpha[nrow(alpha),],
                 "conv"=message, "iter"=t, "time" = time, "time0" = time0)
  return(result)
}


