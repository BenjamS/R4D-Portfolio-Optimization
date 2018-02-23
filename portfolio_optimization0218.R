library(ggplot2)
library(gridExtra)
library(stats)
library(pracma)
library(nleqslv)
source('D:/OneDrive - CGIAR/Documents/EUwCconstr.R')
source('D:/OneDrive - CGIAR/Documents/coloredNoise.R')
# source('~/EUwCconstr.R')
# source('~/coloredNoise.R')
#============================================
#N risky assets
#============================================
m_b <- c(-0.37, -0.55, -0.78, -0.83)
cv_b <- c(0.65, 0.56, 0.79, 0.43) * sign(m_b)
s_b <- cv_b * m_b
m_lA <- c(-.65, -0.75, -0.53, -0.3)
cv_lA <- c(0.43, 0.81, 0.51, 0.3) * sign(m_lA)
s_lA <- cv_lA * m_lA
#--------------------------------------------
n_proj <- length(m_b)
onevec <- rep(1, n_proj)
#--------------------------------------------
s2_lA <- s_lA^2
E_lA <- sum(m_lA)
#--------------------------------------------
noisevec <- coloredNoise(n_proj^2, a = -2, normalize = T, graph = T)
Q <- matrix(exp(noisevec), n_proj, n_proj)
Corr <- Q %*% t(Q)
Corr <- Corr / max(Corr) * .5
diag(Corr) <- rep(1, n_proj)
#print(Corr)
S_b <- diag(s_b) %*% Corr %*% diag(s_b)
Sb_inv <- round(solve(S_b), 5)
#============================================
C_targ <- 1
Vb_targ <- 0.6
#--------------------------------------------
rootfn <- T
quietly <- F
in_vec <- runif(n_proj + 1)
out <- nleqslv(in_vec, minV_EnetUconstr, jac = NULL, Vb_targ, m_b, S_b, rootfn, quietly)
rootfn <- F
in_star <- out$x
outstar <- minV_EnetUconstr(in_star, Vb_targ, m_b, S_b, rootfn, quietly)
wstar <- in_star[1:n_proj]
Cost <- sum(wstar)
Cost
#============================================
pwr10 <- 2
range_V <- c(0, 2.5)
first_i <- range_V[1] * 10^pwr10
last_i <- range_V[2] * 10^pwr10
#---------------
EUb_vec <- c()
Vb_vec <- c()
MD_vec <- c()
lNEU_vec <- c()
#lC_vec <- c()
C_vec <- c()
wstar_list <- list()
slackvec_list <- list()
negpos <- c()
t <- 0
for(i in first_i:last_i){
  #----
  Vb_targ <- i / 10^pwr10
  #----
  for(j in 1:40){
  rootfn <- T
  quietly <- T
  in_vec <- runif(n_proj + 1)
  out <- try(nleqslv(in_vec, minV_EnetUconstr, jac = NULL, Vb_targ, m_b, S_b, rootfn, quietly))
  if(inherits(out, "try-error")){print("gotta skip"); next()}
  if(max(abs(out$fvec)) > 10^-5){print("max iter with no convergence");next()}
  #----
  instar <- out$x
  #----
  rootfn <- F
  quietly <- T
  outstar <- minV_EnetUconstr(instar, Vb_targ, m_b, S_b, rootfn, quietly)
  #----
  EUb_try <- round(outstar[[2]], 5)
  if(EUb_try %in% EUb_vec){print("already have"); next()}
  #----
  t <- t + 1
  #----
  wstar <- instar[1:n_proj]
  l_NEU <- instar[n_proj + 1]
  #l_C <- instar[n_proj + 2]  
  slackvec_list[[t]] <- outstar$fvec
  wstar_list[[t]] <- wstar
  Vb_vec[t] <- Vb_targ
  EUb_vec[t] <- EUb_try
  lEU_vec[t] <- l_EU
  #lC_vec[t] <- l_C
  MD_vec[t] <- outstar[[3]]
  C_vec[t] <- sum(wstar)
  #---
  if(sum(wstar < 0) > 0){negpos[t] <- "negative"}else{negpos[t] <- "positive"}
  }
  
}

risk_vec <- sqrt(Vb_vec)
wstars <- wstar_list
df <- data.frame(EUb = EUb_vec, Risk = risk_vec, lNEU = lEU_vec, MD = MD_vec, C = C_vec, negpos, Vb = Vb_vec)
df$NEUb <- 1 / 9.624643 * (df$EUb - log(df$C))
ind_negw <- which(df$negpos == "negative")
n_neg <- length(ind_negw)
n_tot <- nrow(df)
frac_neg <- n_neg / n_tot
print(frac_neg)
df <- df[-ind_negw, ]
wstars <- wstars[-ind_negw]
ind_MD <- which(df$MD > 10)
length(ind_MD)
df <- df[-ind_MD, ]
wstars <- wstars[-ind_MD]
ind_rm <- which(df$lrat < -10)
length(ind_rm)
df <- df[-ind_rm, ]
wstars <- wstars[-ind_rm]
#------------------
ind <- which(df$NEUb == max(df$NEUb))
data.frame(opt_NEUb = df$NEUb[ind], opt_Risk = df$Risk[ind], opt_NEUbRiskRat = df$NEUb[ind] / df$Risk[ind], opt_lNEU = df$lNEU[ind])
ggplot(df, aes(x = Risk, y = NEUb, color = MD)) + geom_point() + geom_vline(xintercept = df$Risk[ind])

df_wStack <- data.frame(Risk = NA, EUb= NA, project = NA, w = NA)
for(i in 1:nrow(df)){df_wStack <- rbind(df_wStack, data.frame(Risk = df$Risk[i], EUb = df$EUb[i], project = c(1:n_proj), w = c(wstars[[i]][1:n_proj])))}
#for(i in 1:nrow(df)){df_wStack <- rbind(df_wStack, data.frame(Risk = df$Risk[i], EUb = df$EUb[i], project = c(1:n_proj), w = c(wstars[[i]][1:n_proj] / df$C[i])))}
df_wStack <- df_wStack[-1, ]

df_wStack$project <- as.character(df_wStack$project)
df_wStack$w <- df_wStack$w / 9.624643
gg <- ggplot(df_wStack, aes(x = Risk, y = w)) + geom_area(aes(fill = project))
gg <- gg + geom_vline(xintercept = df$Risk[ind])
gg

wstar_opt <- wstars[[ind]]
C_opt <- sum(wstar_opt)
































































#============================================
#============================================
#============================================
#N risky assets with return and budget constraint
#============================================
m_u <- c(0.43, 0.25, 0.38)
m_lA <- c(-.65, -0.75, -0.53)
#m_lA <- c(.02, 0.03, 0.05)
cv_u <- c(0.15, 0.46, 0.29)
cv_lA <- c(0.23, 0.81, 0.51) * sign(m_lA)
s_u <- cv_u * m_u
#s2_u <- s_u^2
s_lA <- cv_lA * m_lA
#s2_lA <- s_lA^2
#E_A <- exp(m_lA + 1 / 2 * s2_lA)
n_proj <- length(m_u)
onevec <- rep(1, n_proj)
#-----------------
noisevec <- coloredNoise(n_proj^2, a = -2, normalize = T, graph = T)
Q <- matrix(exp(noisevec), n_proj, n_proj)
Corr <- Q %*% t(Q)
Corr <- Corr / max(Corr) * .98
diag(Corr) <- rep(1, n_proj)
#print(Corr)
S_u <- diag(s_u) %*% Corr %*% diag(s_u)
Su_inv <- round(solve(S_u), 5)
#-----------------
noisevec <- coloredNoise(n^2, a = -2, normalize = T, graph = T)
Q <- matrix(exp(noisevec), n, n)
Corr <- Q %*% t(Q)
Corr <- Corr / max(Corr) * .98
diag(Corr) <- rep(1, n)
S_lA <- diag(s_lA) %*% Corr %*% diag(s_lA)
V_lA <- t(onevec) %*% S_lA %*% onevec
#-----------------


# l <- eigen(Corr)$value
# e <- eigen(Corr)$vector
# print(l)
# print(e)
# #colMeans(e)
# Y <- diag(1 / s) %*% e
# W <- apply(Y, 2, function(x) x / sum(x))
# colSums(W)
# w_mkt <- W[, 1]
# round(t(w_mkt) %*% S %*% w_mkt, 4)
# round(t(W[, 2]) %*% S %*% w_mkt, 4)
# W_inv <- solve(W)
# m_effective <- W_inv %*% c(1, rep(0, (n - 1)))
# m_effective
# V_vec <- diag(t(W) %*% S %*% W)
# ER_vec <- W %*% m_effective
# V_mkt <- V_vec[1]
# ER_mkt <- ER_vec[1]
#-----
