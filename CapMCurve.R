CapMCurve <- function(in_vec, EU_targ, m_b, S_b, r_f, Budg, E_lA, rootfn = T, quietly = T){
  n_proj <- length(m_b)
  w <- in_vec[1:n_proj]
  l_EU <- in_vec[n_proj + 1]
  #--------------------------------
  naba_Uf <- -w^2 * r_f
  naba_Ur <- m_b
  naba_EU <- naba_Uf + naba_Ur
  naba_V <- 2 * S_b %*% (1 / w)
  naba_mat <- cbind(naba_V, naba_EU)
  l_V <- 1
  l_vec <- c(l_V, l_EU)
  naba_L <- naba_mat %*% l_vec
  E_Uf <- (Budg - t(w) %*% rep(1, n_proj)) * r_f
  E_Ur <- E_lA + t(1 / w) %*% m_b
  E_U <- E_Ur + E_Uf
  slack_EU <- E_U - EU_targ
  slack <- c(naba_L, slack_EU)
  #--------------------------------
  if(rootfn == T){
    if(quietly == F){
      w_and_l <- c(w, l_EU)
      print(data.frame(slack, w_and_l))
    }
    return(slack)
  }else{
    Sb_inv <- round(solve(S_b), 6)
    MD <- t(naba_EU) %*% Sb_inv %*% naba_EU
    MD_wl <- l_EU * t(naba_EU) %*% Sb_inv %*% naba_EU * l_EU
    V <- t(1 / w) %*% S_b %*% (1 / w)
    if(quietly == F){
      w_and_l <- c(w, l_EU)
      print(data.frame(slack, w_and_l))
      print(data.frame(E_U, V, MD, MD_wl))
    }
    outlist <- list(V, E_U, MD_wl, w)
    return(outlist)
  }
  #--------------------------------
}