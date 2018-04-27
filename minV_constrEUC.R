minV_constrEUC <- function(in_vec, Vend_targ, C_targ, m_h, S_h, rootfn = T, quietly = T){
  #in_vec <- in_star
  n_proj <- length(m_h)
  lw <- in_vec[1:n_proj]
  w <- exp(lw)
  #print(in_vec)
  l_EUr <- in_vec[n_proj + 1]
  l_C <- in_vec[n_proj + 2]
  C <- sum(w)
  #--------------------------------
  nablw_EUr <- m_h
  nablw_C <- w
  nablw_V <- 2 * S_h %*% lw
  nablw_mat <- cbind(nablw_V, nablw_EUr, nablw_C)
  l_V <- 1
  l_vec <- c(l_V, l_EUr, l_C)
  nablw_L <- nablw_mat %*% l_vec
  Vend <- t(lw) %*% S_h %*% lw
  slack_Vend <- Vend - Vend_targ
  slack_C <- C - C_targ
  slack <- c(nablw_L, slack_Vend, slack_C)
  #--------------------------------
  if(quietly == F){
    Sh_inv <- round(solve(S_h), 6)
    #--lambdas slack
    Mert <- t(nablw_mat[, c(2:3)]) %*% Sh_inv %*% nablw_mat[, c(2:3)]
    Mert_inv <- round(solve(Mert), 6)
    MD <- t(l_vec[2:3]) %*% Mert_inv %*% l_vec[2:3]
    EUr_end <- m_h %*% lw
    Uc <- w %*% lw
    MD_test <- -2 * (l_vec[2:3] %*% t(nablw_mat[, c(2:3)]) %*% lw)
    MD_check <- MD - MD_test
    MD <- MD_test
    l_effective <- -2 * Mert_inv %*% c(EUr_end, Uc)
    slack_l <- 1 / l_V * l_vec[2:3] - l_effective
    #--wgts slack
    w_effective <- exp(-1 / 2 * Sh_inv %*% (1 / l_V * (nablw_mat[, c(2:3)] %*% l_vec[2:3])))
    slack_w <- w - w_effective
    slack_w_and_l <- c(slack_w, slack_l)
    w_and_l <- c(w, l_EUr, l_C)
    #--Second order condition (must be < 0 for a max(EUr))
    SOC <- t(lw) %*% (2 * S_h + l_C * diag(w)) %*% lw
  }
  #--------------------------------
  if(rootfn == T){
    if(quietly == F){
      print(data.frame(slack, w_and_l, slack_w_and_l))
    }
    return(slack)
  }else{
    if(quietly == F){
      frontier_root <- t(lw) %*% nablw_mat %*% l_vec
      print(data.frame(slack, w_and_l, slack_w_and_l))
      print(data.frame(frontierRoot = frontier_root, MD_check = MD_check))
      print(data.frame(EUr_end, Vend, Uc, MD, SOC))
      #---
    }else{
      Sh_inv <- round(solve(S_h), 6)
      Mert <- t(nablw_mat[, c(2:3)]) %*% Sh_inv %*% nablw_mat[, c(2:3)]
      Mert_inv <- round(solve(Mert), 6)
      #MD <- t(l_vec[2:3]) %*% Mert_inv %*% l_vec[2:3]
      MD <- -2 * (l_vec[2:3] %*% t(nablw_mat[, c(2:3)]) %*% lw)
      EUr_end <- t(lw) %*% m_h
      SOC <- t(lw) %*% (2 * S_h + l_C * diag(w)) %*% lw
    }
    outlist <- list(Vend, EUr_end, SOC, MD)
    return(outlist)
  }
  #--------------------------------
}
