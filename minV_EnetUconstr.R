minV_EnetUconstr <- function(in_vec, Vb_targ, m_b, S_b, rootfn = T, quietly = T){
  n_proj <- length(m_b)
  w <- in_vec[1:n_proj]
  wi <- 1 / w
  l_NEU <- in_vec[n_proj + 1]
  C <- sum(w)
  #--------------------------------
  nabwi_NEU <- m_b + 1 / C * w^2
  #nabwi_C <- -w^2
  nabwi_V <- 2 * S_b %*% wi
  nabwi_mat <- cbind(nabwi_V, nabwi_NEU)
  l_V <- 1
  l_vec <- c(l_V, l_NEU)
  nabwi_L <- nabwi_mat %*% l_vec
  V_b <- t(wi) %*% S_b %*% wi
  slack_Vb <- Vb_targ - V_b
  slack <- c(nabwi_L, slack_Vb)
  #--------------------------------
  if(quietly == F){
    w_and_l <- c(w, l_NEU)
    Sb_inv <- round(solve(S_b), 6)
    #--lambdas slack
    EU_b <- t(wi) %*% m_b
    q <- EU_b + 1
    MD <- t(nabwi_mat[, 2]) %*% Sb_inv %*% nabwi_mat[, 2]
    #l_endog <- -2 * V_b / q
    l_endog <- as.numeric(-2 * q / MD)
    slack_l <- l_vec[2] - l_endog
    #--wgts slack
    wi_endog <- -l_endog / 2 * Sb_inv %*% nabwi_mat[, 2]
    slack_wi <- wi - wi_endog
    slack_wi_and_l <- c(slack_wi, slack_l)
  }
  #--------------------------------
  if(rootfn == T){
    if(quietly == F){
      print(data.frame(slack, w_and_l, slack_wi_and_l))
    }
    return(slack)
  }else{
    if(quietly == F){
      frontier_root <- 2 * V_b + l_NEU * q
      frontier_root2 <- q^2 - V_b * MD
      print(data.frame(slack, w_and_l, slack_wi_and_l))
      print(data.frame(frontierRoot = frontier_root, frontierRoot2 = frontier_root2))
      print(data.frame(EU_b, V_b, C, MD))
      #---
    }else{
      MD <- t(nabwi_mat[, 2]) %*% Sb_inv %*% nabwi_mat[, 2]
      EU_b <- t(wi) %*% m_b
    }
    outlist <- list(V_b, EU_b, MD)
    return(outlist)
  }
  #--------------------------------
}
