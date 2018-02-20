ENetU <- function(in_vec, Vb_targ, m_b, S_b, E_lA, rootfn = T, quietly = T){
  n_proj <- length(m_b)
  w <- in_vec[1:n_proj]
  wi <- 1 / w
  l_Enu <- in_vec[n_proj + 1]
  C_r <- sum(w)
  #--------------------------------
  nabwi_C <- w^2
  nabwi_EUr <- m_b
  nabwi_Enu <- nabwi_EUr + nabwi_C
  nabwi_V <- 2 * S_b %*% wi
  nabwi_mat <- cbind(nabwi_V, nabwi_Enu)
  l_V <- 1
  l_vec <- c(l_V, l_Enu)
  nabwi_L <- nabwi_mat %*% l_vec
  E_Ur <- E_lA + t(wi) %*% m_b
  E_NU <- E_Ur - C_r
  V_b <- t(wi) %*% S_b %*% wi
  slack_Vb <- Vb_targ - V_b
  slack <- c(nabwi_L, slack_Vb)
  #--------------------------------
  if(rootfn == T){
    if(quietly == F){
      w_and_l <- c(w, l_Enu)
      Sb_inv <- round(solve(S_b), 6)
      MD <- t(nabwi_Enu) %*% Sb_inv %*% nabwi_Enu
      term <- E_Ur - E_lA + C_r
      l_Enu_test1 <- -2 * term / MD
      slack_wi <- wi + 1 / 2 * l_Enu * Sb_inv %*% nabwi_Enu
      slack_lEnu <- l_Enu_test1 - l_Enu
      slack_wi_and_l <- c(slack_wi, slack_lEnu)
      print(data.frame(slack, w_and_l, slack_wi_and_l))
    }
    return(slack)
  }else{
    Sb_inv <- round(solve(S_b), 6)
    MD <- t(nabwi_Enu) %*% Sb_inv %*% nabwi_Enu
    if(quietly == F){
      w_and_l <- c(w, l_Enu)
      term <- E_Ur - E_lA + C_r
      l_Enu_test1 <- -2 * V_b / term
      l_Enu_test2 <- -2 * term / MD
      slack_wi <- wi + 1 / 2 * l_Enu * Sb_inv %*% nabwi_Enu
      slack_lEnu <- l_Enu_test1 - l_Enu
      slack_wi_and_l <- c(slack_wi, slack_lEnu)
      print(data.frame(slack, w_and_l, slack_wi_and_l))
      print(data.frame(E_NU, V_b, MD))
      #---
      frontier_root <- term + sign(l_Enu) * sqrt(V_b * MD)
      print(data.frame(l_Enu = l_Enu, lEnutest1 = l_Enu_test1, lEnutest2 = l_Enu_test2, frontierRoot = frontier_root))
      #---
    }
    outlist <- list(V_b, E_NU, C_r, MD)
    return(outlist)
  }
  #--------------------------------
}