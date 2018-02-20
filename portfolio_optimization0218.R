library(ggplot2)
library(gridExtra)
library(stats)
library(pracma)
library(nleqslv)
source('~/ENetU.R')
source('~/coloredNoise.R')
source('~/ENetU.R')
#============================================
#N risky assets and 1 risk free asset
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
#--------------------------------------------
r_f <- 0.05
Budg <- 1
#--------------------------------------------
Vb_targ <- 0.2
#----
rootfn <- T
quietly <- F
in_vec <- runif(n_proj + 1)
out <- nleqslv(in_vec, ENetU, jac = NULL, Vb_targ, m_b, S_b, E_lA, rootfn, quietly)
rootfn <- F
in_star <- out$x
outstar <- ENetU(in_star, Vb_targ, m_b, S_b, E_lA, rootfn, quietly)
wstar <- in_star[1:n_proj]
Cost_r <- sum(wstar)
Cost_r
#--------------------------------------------
pwr10 <- 2
range_EU <- c(0, 2)
first_i <- range_EU[1] * 10^pwr10
last_i <- range_EU[2] * 10^pwr10
#---------------
Enu_vec <- c()
Vb_vec <- c()
MD_vec <- c()
lEnu_vec <- c()
Costr_vec <- c()
RtoR_vec <- c()
EUr_vec <- c()
wstar_list <- list()
slackvec_list <- list()
negpos <- c()
t <- 0
for(i in first_i:last_i){
  #----
  Vb_targ <- i / 10^pwr10
  #----
  rootfn <- T
  quietly <- T
  in_vec <- runif(n_proj + 1)
  out <- try(nleqslv(in_vec, ENetU, jac = NULL, Vb_targ, m_b, S_b, E_lA, rootfn, quietly))
  if(inherits(out, "try-error")){print("gotta skip"); next()}
  if(max(abs(out$fvec)) > 10^-5){print("max iter with no convergence");next()}
  #----
  instar <- out$x
  #----
  rootfn <- F
  quietly <- T
  outstar <- ENetU(instar, Vb_targ, m_b, S_b, E_lA, rootfn, quietly)
  #----
  Enu_try <- round(outstar[[2]], 5)
  if(Enu_try %in% Enu_vec){print("already have"); next()}
  #----
  t <- t + 1
  #----
  wstar <- instar[1:n_proj]
  l_Enu <- instar[n_proj + 1]
  slackvec_list[[t]] <- outstar$fvec
  wstar_list[[t]] <- wstar
  Vb_vec[t] <- Vb_targ
  Enu_vec[t] <- Enu_try
  lEnu_vec[t] <- l_Enu
  MD_vec[t] <- outstar[[3]]
  #---
  Costr <- sum(wstar)
  Costr_vec[t] <- Costr
  RtoR_vec[t] <- (Enu_try - Costr) / Vb_targ
  EUr_vec[t] <- Enu_try + Costr
  #---
  if(sum(wstar < 0) > 0){negpos[t] <- "negative"}else{negpos[t] <- "positive"}
  
}

risk_vec <- sqrt(Vb_vec)
wstars <- wstar_list
Sb_inv <- round(solve(S_b), 6)
MD_unconst <- as.numeric(t(m_b) %*% Sb_inv %*% m_b)
EU_unconst <- -risk_vec * sqrt(MD_unconst) + E_lA
df <- data.frame(Enu = Enu_vec, EUr = EUr_vec, Risk = risk_vec, lEnu = lEnu_vec, MD = MD_vec, RtoR = RtoR_vec, Costr = Costr_vec, EUunc = EU_unconst, negpos, Vb = Vb_vec)
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
#------------------
gg <- ggplot(df, aes(x = Risk, y = Enu, color = MD)) + geom_point()
gg <- gg + geom_point(aes(x = Risk, y = EUr, color = MD))
gg

# df$diff <- df$EUr - df$EUunc
# ggplot(df, aes(x = EUr, y = diff, color = MD)) + geom_point()

# ind_sort <- order(df$V)
# w_sort <- wvec_list[ind_sort]
# V_sort <- df$V[ind_sort]
df_wStack <- data.frame(Risk = NA, Enu= NA, project = NA, w = NA)
for(i in 1:nrow(df)){df_wStack <- rbind(df_wStack, data.frame(Risk = df$Risk[i], Enu = df$Enu[i], project = c(1:n_proj), w = c(wstars[[i]][1:n_proj])))}
#for(i in 1:nrow(df)){df_wStack <- rbind(df_wStack, data.frame(Risk = df$Risk[i], Enu = df$Enu[i], project = c(1:n_proj), w = c(wstars[[i]][1:n_proj] / df$Costr[i])))}
df_wStack <- df_wStack[-1, ]

df_wStack$project <- as.character(df_wStack$project)
gg <- ggplot(df_wStack, aes(x = Risk, y = w)) + geom_area(aes(fill = project))
gg




































































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
