VarAdjNetBen <- function(in_w, u_vec, s_df, A_df, Vmat, Vmat_inv,
                           C_grad, Cst, lambda_C, discnt_df)
{
  #  if(abs(mean(in_w)) > 10^2){in_w <- runif(n_proj, 0, 10)}
  if(sum(in_w) != 1){in_w <- in_w / sum(in_w)}
  V <- t(in_w) %*% Vmat %*% in_w
  u <- u_vec
  n_dims <- ncol(b_likely_df)
  n_proj <- length(u)
  q <- c()
  PB <- c()
  Pb_grad_list <- list()
  for(i in 1:n_dims)
  {
    A <- A_df[, i]
    s <- s_df[, i]
    discnt <- discnt_df[, i]
    b_targ <- unlist(A * exp(-u / in_w))
    m <- log(b_targ) + s^2
    dmdw <- as.vector(unlist(log(b_targ / A)^2 / u))
    Eb <- unlist(exp(m + s^2 / 2))
    Pb <- diag(discnt) %*% Eb
    #print(diag(dmdw) %*% Pb)
    Pb_grad <- as.vector(diag(dmdw) %*% Pb)
    Pb_grad_list[[i]] <- Pb_grad
    q[i] <- t(Pb) %*% -log(b_targ / A)
    #q[i] <- t(Pb_grad) %*% in_w
    #    q[i] <- t(Pb_grad) %*% (-u / log(b_targ / A))
    #    q[i] <- beta * (t(Pb) %*% c(1, 1, 1))
    PB[i] <- sum(Pb)
  }
  Pb_grad_mat <- do.call(cbind, Pb_grad_list)
  Grad_mat <- as.matrix(cbind(Pb_grad_mat, C_grad))
  Merton_mat <- round(t(Grad_mat) %*% Vmat_inv %*% Grad_mat, 4)
  qCV <- c(q, Cst, 2 * V)
  augMerton_mat <- cbind(Merton_mat, rep(0, (n_dims + 1)))
  augMerton_mat <- rbind(augMerton_mat, qCV)
  augMerton_mat <- round(augMerton_mat,4)
  augMerton_mat_inv <- solve(augMerton_mat)
  qC0 <- c(q, Cst, 0)
  xx <- -2 * augMerton_mat_inv %*% qC0
  
  
  o1 <- (t(Grad_mat) %*% Vmat_inv %*%  Grad_mat %*% xx[1:(n_dims + 1)]) + 2 * (t(Grad_mat) %*% in_w)
  
  alt_w <- -1/2 * (Vmat_inv %*%  Grad_mat %*% xx[1:(n_dims + 1)])
  
  o2 <- (t(Grad_mat) %*% Vmat_inv %*%  Grad_mat %*% xx[1:(n_dims + 1)]) + 2 * (t(Grad_mat) %*% alt_w)
  
  dif_w <- in_w - alt_w
  
  # print(o1);print(o2)
  # print(data.frame(in_w, alt_w, dif_w))
  # print("Lambdas div'd by lambda_V. Last entry should = 1.")
  # print(xx)
  
  
  if(abs(xx[n_dims + 2] - 1) > 10^-5){solflag <- 0; outlist <- list(solflag, alt_w)}
  else
  {
    print("Convergence!")
    solflag <- 1
    lambdas_divby_lambdaV <- xx
    #Derive lambda_V by applying convention lambda_C = -1
    lambdaC_divby_lambdaV <- lambdas_divby_lambdaV[n_dims + 1]
    lambda_V <- - 1 / lambdaC_divby_lambdaV
    lambdas <- round(lambda_V * lambdas_divby_lambdaV, 4)
    raNB <- round(c(PB, Cst, V) %*% lambdas, 4)
    NB <- raNB - lambda_V * V
    Benefit <- round(PB %*% lambdas[1:n_dims], 4)
    frontier_check <- round(qCV %*% lambdas, 4)
    print("Frontier check (should = 0)")
    print(frontier_check)
    # w_check / in_w
    # print("w check (should = 0)")
    # print(w_check)
    # M_check <- augMerton_mat %*% lambda_vec + 2 * lambda_V * qC0
    # print("Merton mat check (should = 0)")
    # print(M_check)
    outlist <- list(solflag, alt_w, lambdas, PB, raNB, Benefit, V, NB)
    
  }
  #---
return(outlist)
  #---
}