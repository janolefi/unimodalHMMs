## packages and utility functions
library(LaMa)
source("utils.R")


## reading in the data
data = readRDS("./narwal_case_study/data/MaxDepth_TimeSeries.RData")
# subtracting the threshold
data$maxdep = data$maxdep - 20
nrow(data)
range(data$maxdep)


## defining color vector
color = c("orange", "deepskyblue", "seagreen2")


## exploring initial values
hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", ylim = c(0, 0.015))
curve(0.8 * dgamma2(x, 20, 30), add = TRUE)
curve(0.1 * dgamma2(x, 360, 50), add = TRUE)
curve(0.1 * dgamma2(x, 490, 40), add = TRUE)



# Parametric model --------------------------------------------------------

# simple model with gamma state-dependent distributions
nll = function(par){
  getAll(par, dat) # makes everything accessible w/o $
  ## state process
  Gamma = tpm(beta0) # mlogit for transitiom probabilities
  delta = stationary(Gamma) # stationary distribution of the state process
  ## state-dependent process
  # parameter transformations
  mu = exp(logmu); REPORT(mu) # reporting to make accessible later
  sigma = exp(logsigma); REPORT(sigma)
  # state-dependent densities
  allprobs = matrix(1, length(maxdep), N)
  ind = which(!is.na(maxdep)) # NA handling
  for(j in 1:N) { # loop over states
    allprobs[ind,j] = dgamma2(maxdep[ind], mu[j], sigma[j])
  }
  # forward algorithm (separately for different tracks and then summed)
  -forward(delta, Gamma, allprobs, ID)
}


## fitting the model

# Define base values (based on initial exploration)
beta0_base <- rep(-2, 6)
logmu_base <- log(c(20, 360, 490))
logsigma_base <- log(c(30, 50, 50))

# Define noise levels
beta0_sd <- 0.5       # Standard deviation for beta0
logmu_sd <- 0.5       # Standard deviation for logmu
logsigma_sd <- 0.5    # Standard deviation for logsigma

# number of random initial values
n_random = 100 # 500

# Seed for RNG
set.seed(123)

# Generate n_random random sets
init_values <- lapply(1:n_random, function(i) {
  list(
    beta0 = rnorm(6, mean = beta0_base, sd = beta0_sd),
    logmu = rnorm(3, mean = logmu_base, sd = logmu_sd),
    logsigma = rnorm(3, mean = logsigma_base, sd = logsigma_sd)
  )
})

# data and hyperparamter object
dat = list(maxdep = data$maxdep, ID = data$ID, N = 3)

# fitting models with different initial values
llks = rep(NA, n_random)
for(i in 1:n_random){
  cat("Iteration", i, "\n")
  thispar = init_values[[i]]
  obj = MakeADFun(nll, thispar, silent = TRUE) # creating automatically diff'able objective
  opt = tryCatch(
    nlminb(obj$par, obj$fn, obj$gr), # optimising it
    error = function(e) list(value = NA)) # handle errors
  llks[[i]] = - opt$objective # saving log-likelihood
  gc() # cleaning up
}

plot(llks)

# refit best model
par = init_values[[which.max(llks)]]
obj = MakeADFun(nll, par)
opt = nlminb(obj$par, obj$fn, obj$gr)

# extract estimated parameters
mod_par = obj$report()

# save estimated model
# saveRDS(mod_par, "./narwal_case_study/mod_normal.rds")
# mod_par = readRDS("./narwal_case_study/mod_normal.rds")

mu = mod_par$mu
sigma = mod_par$sigma
delta_par = mod_par$delta


## plotting estimated state-dependent distributions
par(mfrow = c(1,1))
N = 3
# zoomed out
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.002), xlim = c(0, 1100), main = "", xlab = "maxdepth", ylab = "density")
for(i in 1:N) curve(delta_par[i] * dgamma2(x, mu[i], sigma[i]), add = TRUE, col = color[i], lwd = 2)
curve(delta_par[1] * dgamma2(x, mu[1], sigma[1]) +
        delta_par[2] * dgamma2(x, mu[2], sigma[2]) +
        delta_par[3] * dgamma2(x, mu[3], sigma[3]), add = TRUE, lwd = 2, lty = 2)




# More complex parametric model -------------------------------------------

# # non-standard model:
# # state 1 follows a gamma distribution
# # states 2 and 3 follow skew-normal distributions
# # this likelihood now only works for 3 states
# nll2 = function(par){
#   getAll(par, dat) # makes everything accessible w/o $
#   ## state process
#   Gamma = tpm(beta0) # mlogit for transitiom probabilities
#   delta = stationary(Gamma) # stationary distribution of the state process
#   ## state-dependent process
#   # parameter transformations
#   mu = exp(logmu); REPORT(mu)
#   sigma = exp(logsigma); REPORT(sigma)
#   REPORT(alpha) # these are the 2 skewness parameters
#   # state-dependent densities
#   allprobs = matrix(1, length(maxdep), 3)
#   ind = which(!is.na(maxdep)) # NA handling
#   allprobs[ind,1] = dgamma2(maxdep[ind], mu[1], sigma[1]) # state 1: gamma
#   # states 2 and 3: skew-normal
#   for(j in 2:3) allprobs[ind,j] = dskewnorm(maxdep[ind], mu[j], sigma[j], alpha[j-1])
#   # forward algorithm
#   -forward(delta, Gamma, allprobs, ID)
# }
# 
# 
# ## fitting the model
# 
# # Define base values (based on initial exploration)
# alpha_base = 1
# # Define noise levels
# alpha_sd = 1
# 
# # number of random initial values
# n_random = 500
# 
# # Seed for RNG
# set.seed(12)
# 
# # Generate n_random random sets
# init_values <- lapply(1:n_random, function(i) {
#   list(
#     beta0 = rnorm(6, mean = beta0_base, sd = beta0_sd),
#     logmu = rnorm(3, mean = logmu_base, sd = logmu_sd),
#     logsigma = rnorm(3, mean = logsigma_base, sd = logsigma_sd),
#     alpha = rnorm(2, mean = alpha_base, sd = alpha_sd)
#   )
# })
# 
# # data and hyperparameter object
# dat = list(maxdep = data$maxdep, ID = data$ID)
# 
# # # fitting models with different initial values
# # llks = rep(NA, n_random)
# # for(i in 1:n_random){
# #   cat("Iteration", i, "\n")
# #   thispar = init_values[[i]]
# #   obj = MakeADFun(nll2, thispar, silent = TRUE) # creating automatically diff'able objective
# #   opt = optim(obj$par, obj$fn, obj$gr, method = "BFGS") # optimising it
# #   llks[[i]] = - opt$value # saving log-likelihood
# #   gc() # cleaning up
# # }
# # 
# # plot(llks)
# # 
# # # refit best model
# # par = init_values[[which.max(llks)]]
# # obj2 = MakeADFun(nll2, par)
# # opt2 = optim(obj2$par, obj2$fn, obj2$gr, method = "BFGS")
# # 
# # # extract estimated parameters
# # mod_par2 = obj2$report()
# 
# # save estimated model
# # saveRDS(mod_par2, "./narwal_case_study/mod_skew.rds")
# mod_par2 = readRDS("./narwal_case_study/mod_skew.rds")
# 
# mu2 = mod_par2$mu
# sigma2 = mod_par2$sigma
# delta_par2 = mod_par2$delta
# alpha = mod_par2$alpha
# 
# 
# ## plotting estimated state-dependent distributions
# par(mfrow = c(1,1))
# hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", 
#      ylim = c(0, 0.005), xlim = c(0, 800), main = "", xlab = "maxdepth", ylab = "density")
# curve(delta_par2[1] * dgamma2(x, mu2[1], sigma2[1]), add = TRUE, col = color[1], lwd = 2)
# for(j in 2:3) curve(delta_par2[j] * dskewnorm(x, mu2[j], sigma2[j], alpha[j-1]), add = TRUE, col = color[j], lwd = 2)
# curve(delta_par2[1] * dgamma2(x, mu2[1], sigma2[1]) +
#         delta_par2[2] * dskewnorm(x, mu2[2], sigma2[2], alpha[1]) +
#         delta_par2[3] * dskewnorm(x, mu2[3], sigma2[3], alpha[2]), add = TRUE, lwd = 2, lty = 2)
# 
# 
# # plotting both parametric ones
# # pdf("narwal_case_study/parametric_models.pdf", width = 8, height = 4)
# par(mfrow = c(1,2))
# N = 3
# color = c("seagreen2", "orange", "deepskyblue")
# hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", 
#      ylim = c(0, 0.002), xlim = c(0, 800), xlab = "maxdepth", ylab = "density",
#      main = "gamma")
# for(i in 1:N) curve(delta_par[i] * dgamma2(x, mu[i], sigma[i]), add = TRUE, col = color[i], lwd = 2)
# curve(delta_par[1] * dgamma2(x, mu[1], sigma[1]) +
#         delta_par[2] * dgamma2(x, mu[2], sigma[2]) +
#         delta_par[3] * dgamma2(x, mu[3], sigma[3]), add = TRUE, lwd = 2, lty = 2)
# 
# color = c("orange", "deepskyblue", "seagreen2")
# hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", 
#      ylim = c(0, 0.002), xlim = c(0, 800), xlab = "maxdepth", ylab = "density",
#      main = "gamma-skewnormal")
# curve(delta_par2[1] * dgamma2(x, mu2[1], sigma2[1]), add = TRUE, col = color[1], lwd = 2)
# for(j in 2:3) curve(delta_par2[j] * dskewnorm(x, mu2[j], sigma2[j], alpha[j-1]), add = TRUE, col = color[j], lwd = 2)
# curve(delta_par2[1] * dgamma2(x, mu2[1], sigma2[1]) +
#         delta_par2[2] * dskewnorm(x, mu2[2], sigma2[2], alpha[1]) +
#         delta_par2[3] * dskewnorm(x, mu2[3], sigma2[3], alpha[2]), add = TRUE, lwd = 2, lty = 2)
# 
# # dev.off()




# Unconstrained nonparametric fit -----------------------------------------

# likelihood function for the unconstrained and constrained model 
# for the unconstrained model, we just set kappa = 0
pnll = function(par){
  getAll(par, dat)
  
  ## state process model
  Gamma = tpm(beta0)
  delta = stationary(Gamma)
  
  ## state-dependent process model
  # multinomial logit for density weights
  alpha = exp(cbind(beta, 0)) # last column is fixed at zero
  alpha = alpha / rowSums(alpha); REPORT(alpha) # reporting for easy access
  # state-dependent densities 
  allprobs = matrix(1, nrow(Z), N)
  ind = which(!is.na(Z[,1])) # NA handling
  # now just linear combinations of precomputed basis functions in Z
  allprobs[ind,] = Z[ind,] %*% t(alpha) # Z %*% alpha[i,] for all states at once
  
  # forward algorithm + smoothness penalty + potential unimodality penalty (for kappa > 0)
  -forward(delta, Gamma, allprobs) + # negative log likelihood
    penalty(beta, S, lambda)# + # smoothness P-spline penalty
    # penalty_uni(beta + logweights, m, kappa, rho = 10) # unimodality penalty
}

## need to do this for random initial values
# number of random initial values
n_random = 100

# set.seed(123)
set.seed(12)

# Generate n_random random sets
init_values <- lapply(1:n_random, function(i) {
  list(
    beta0 = rnorm(6, mean = beta0_base, sd = 0.5),
    logmu = rnorm(3, mean = logmu_base, sd = 0.5),
    logsigma = rnorm(3, mean = logsigma_base, sd = 0.5)
  )
})

# setting up the basis and penalty
k = 35 # number of basis functions
degree = 3 # degree of the B-spline basis
pow = 1 # equidistant spacing
# custom knots
### increase the right knot location to 1.2 or sth
knots = seq(0, (1.1 * max(data$maxdep, na.rm = T))^pow, length = k - degree + 1)^(1/pow)

llks = rep(NA, n_random)
for(i in 1:n_random){
  thispar = init_values[[i]]
  cat("\n\nInitial value set", i, "\n\n")
  sDens = smooth_dens_construct(data["maxdep"],
                                list(maxdep = list(mean = exp(thispar$logmu), 
                                                   sd = exp(thispar$logsigma))),
                                type = "positive",
                                k = k,
                                knots = list(maxdep = knots))
  
  Z = sDens$Z$maxdep # design matrix for the state-dependent densities
  S = sDens$S$maxdep # penalty matrix for the state-dependent densities
  Z_p = sDens$Z_predict$maxdep # prediction design matrix for the state-dependent densities
  xseq = sDens$xseq$maxdep # prediction grid for the state-dependent densities
  w = sDens$basis$maxdep$w
  beta = sDens$coef$maxdep
  
  # initial parameter list
  thispar = list(beta0 = rep(-2, 6),
                 beta = beta) # initial coefficients returned by smooth_dens_construct()
  
  # data and hyperparameter list
  thisdat = list(Z = Z, S = S, # design and penalty matrix
                 N = 3, # number of states
                 lambda = rep(1, 3), # initial smoothenss penalty parameters
                 m = 1:3,
                 # C = construct_C(1:3, 3, k), 
                 kappa = 0,# no constraint (but we need to put some C matrices)
                 logweights = matrix(log(w[-length(w)]) - log(w[length(w)]), nrow(beta), ncol(beta), byrow = TRUE)) 
  
  ## fitting the model using the qREML / extended Fellner-Schall method
  thismod = tryCatch(
    qreml(pnll, thispar, thisdat, 
          random = "beta", # telling the model what are spline coefficients/ random effects
          silent = 1, # more printing
          alpha = 0.3, 
          conv_crit = "relchange"), # exponential smoothing of the outer optimisation
    error = function(e) list(llk = NA))
  
  llks[i] = thismod$llk
  
  gc() # cleaning up
}

plot(llks)

# optpar = init_values[[which.max(llks)]]
# # optpar = init_values[[1]]
# sDens = smooth_dens_construct(data["maxdep"],
#                               list(maxdep = list(mean = exp(optpar$logmu), 
#                                                  sd = exp(optpar$logsigma))),
#                               type = "positive",
#                               k = k,
#                               knots = list(maxdep = knots))

# saveRDS(sDens, "./narwal_case_study/mod_spline_sDens.rds")
sDens = readRDS("./narwal_case_study/mod_spline_sDens.rds")


Z = sDens$Z$maxdep # design matrix for the state-dependent densities
S = sDens$S$maxdep # penalty matrix for the state-dependent densities
Z_p = sDens$Z_predict$maxdep # prediction design matrix for the state-dependent densities
xseq = sDens$xseq$maxdep # prediction grid for the state-dependent densities
w = sDens$basis$maxdep$w
beta = sDens$coef$maxdep

# plot(x[ord], Z[ord,1], type = "l", xlim = c(0, 1100))
# for(i in 2:ncol(Z)) lines(x[ord], Z[ord,i], col = i)
# for(i in 1:ncol(Z_p)) lines(xseq, Z_p[,i], col = i, lty = 2)


# initial parameter list
par = list(beta0 = rep(-2, 6),
           beta = beta) # initial coefficients returned by smooth_dens_construct()

# data and hyperparameter list
dat = list(Z = Z, S = S, # design and penalty matrix
           N = 3, # number of states
           ID = data$ID, # track ID
           lambda = rep(25, 3), # initial smoothenss penalty parameters
           C = construct_C(1:3, 3, k), kappa = 0,
           logweights = matrix(log(w[-length(w)]) - log(w[length(w)]), nrow(beta), ncol(beta), byrow = TRUE)) # no constraint (but we need to put some C matrices)

# ## fitting the model using the qREML / extended Fellner-Schall method
# mod = qreml(pnll, par, dat, 
#             random = "beta", # telling the model what are spline coefficients/ random effects
#             silent = 0, # more printing
#             alpha = 0.3) # exponential smoothing of the outer optimisation

# saveRDS(mod, "./narwal_case_study/mod_spline.rds")
mod = readRDS("./narwal_case_study/mod_spline.rds")


# extracting estimated coefficients
alpha = mod$alpha # coefficients of the state-dependent densities
delta = mod$delta # stationary distribution of the state process
# predict state-dependent densities on fine grid
dens = t(delta * t(Z_p %*% t(alpha))) # all (weighted by delta) at once
# saveRDS(dens, "./narwal_case_study/mod_spline_dens.rds")

## plotting estimated state-dependent distributions
par(mfrow = c(1,1))
hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", 
     ylim = c(0, 0.005), main = "unconstrained", 
     xlab = "max depth", ylab = "density", xlim = c(0, 800))
for(i in 1:3) lines(xseq, dens[,i], col = color[i], lwd = 2)
lines(xseq, rowSums(dens), lty = 2, lwd = 1.5)
abline(v = sDens$basis$maxdep$basis_pos)

hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", 
     ylim = c(0, 0.0005), main = "unconstrained", 
     xlab = "max depth", ylab = "density", xlim = c(500, 1100))
for(i in 1:3) lines(xseq, dens[,i], col = color[i], lwd = 2)
lines(xseq, rowSums(dens), lty = 2, lwd = 1.5)


beta = mod$par$beta
# label switching
par(mfrow = c(2,3))
for(i in c(1,3,2)) {
  plot(beta[i,])
  # abline(h = 0)
}
for(i in c(1,3,2)) {
  plot(alpha[i,])
  # abline(h = 0)
}

# set_m_grid(beta)
# does not work here yet
# pick candidates manually
m1 = 1 # fix this at one
# base the potenial locations for 2 and 3 on unconstrained fit
m2 = 11:15
m3 = 15:19
m_grid = expand.grid(m1 = m1, m2 = m2, m3 = m3)
dat$kappa = 10000
dat$lambda = mod$lambda * 3
# dat$w = sDens$basis$maxdep$w

# llks = rep(NA, nrow(m_grid))
# for(i in 1:nrow(m_grid)){
#   cat("\nPosition", i, "of", nrow(m_grid), "\n")
#   thism = m_grid[i,]
#   dat$C = construct_C(thism, 3, k)
#   thismod = qreml(pnll, par, dat, random = "beta", silent = 1, alpha = 0.4)
#   llks[i] = thismod$llk
#   
#   gb() # cleaning up
# }
# plot(llks)

# m_opt = m_grid[which.max(llks), ]
m_opt = c(1, 13, 17)
dat$C = construct_C(m_opt, 3, k, exclude_last = T)
dat$kappa = 1e3
# par$beta = beta
# par$beta = sDens$coef$maxdep
par$beta = mod$par$beta
dat$lambda = mod$lambda * 2
# map = list(lambda = rep(NA, 3))
map = NULL

# mod_c = qreml(pnll, par, dat, random = "beta",
#               silent = 1, alpha = 0.4, map = map)

# saveRDS(mod_c, "./narwal_case_study/mod_spline_constrained.rds")
mod_c = readRDS("./narwal_case_study/mod_spline_constrained.rds")

alpha_c = mod_c$alpha
delta_c = mod_c$delta
dens_c = t(delta_c * t(Z_p %*% t(alpha_c)))
# saveRDS(dens_c, "./narwal_case_study/mod_spline_constrained_dens.rds")

# compare stationary distributions
delta
delta_c

states = viterbi(mod = mod)
states_c = viterbi(mod = mod_c)


# pdf("narwal_case_study/spline_models.pdf", width = 8, height = 3.5)
# state-dependent and marginal distributions
par(mfrow = c(1,4))
ord = c(1,3,2)
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.002), main = "unconstrained", 
     xlab = "max depth", ylab = "density", xlim = c(0, 1100))
for(i in 1:3) lines(xseq, dens[,i], col = color[ord[i]], lwd = 2)
lines(xseq, rowSums(dens), lty = 2, lwd = 1.5)

hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.002), main = "constrained", 
     xlab = "max depth", ylab = "density", xlim = c(0, 1100))
for(i in 1:3) lines(xseq, dens_c[,i], col = color[i], lwd = 2)
lines(xseq, rowSums(dens_c), lty = 2, lwd = 1.5)


# observed sequence with decoded states
# par(mfrow = c(2,1))
# ind = 1:200
ind = 1:nrow(data)
# plot(data$maxdep[ind], col = color[states[ind]], type = "h", bty = "n", ylab = "max depth")
# plot(data$maxdep[ind], col = color[states_c[ind]], type = "h", bty = "n", ylab = "max depth")
# difference in allocation of states 1 or 2 for very large observations
# what's going on here?:


# state-dependent and marginal distributions, tail only
# par(mfrow = c(1,2))
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.0002), main = "unconstrained", 
     xlab = "max depth", ylab = "density", xlim = c(500, 1100))
for(i in 1:3) lines(xseq, dens[,i], col = color[ord[i]], lwd = 2)
lines(xseq, rowSums(dens), lty = 2, lwd = 1.5)

hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.0002), main = "constrained", 
     xlab = "max depth", ylab = "density", xlim = c(500, 1100))
for(i in 1:3) lines(xseq, dens_c[,i], col = color[i], lwd = 2)
lines(xseq, rowSums(dens_c), lty = 2, lwd = 1.5)
# the unconstrained one actually has a second mode in state 2 in the tail, which leads to bad state allocation
# very very large observations are classified as state 2 not state 3. The unimodality constraint fixes this, yeyy!

legend("topright", col = color, lwd = 2,
       legend = paste("state", 1:3), bty = "n")
# dev.off()




## separate plots

# zoomed out
par(mfrow = c(1,1), xpd = F)
ord = c(1,3,2)

mod_par = readRDS("./narwal_case_study/mod_normal.rds")
mu = mod_par$mu
sigma = mod_par$sigma
delta_par = mod_par$delta

# gamma model
# pdf("./narwal_case_study/joint_figure/parametric.pdf", width = 6, height = 5)
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.002), main = "parametric", 
     xlab = "max depth", ylab = "density", xlim = c(0, 1000))
for(i in 1:N) {
  curve(delta_par[i] * dgamma2(x, mu[i], sigma[i]), add = TRUE, col = color[i], lwd = 2,
        from = 0, to = 1000)
}
# curve(delta_par[1] * dgamma2(x, mu[1], sigma[1]) +
#         delta_par[2] * dgamma2(x, mu[2], sigma[2]) +
#         delta_par[3] * dgamma2(x, mu[3], sigma[3]), add = TRUE, lwd = 2, lty = 2)
legend("topright", col = color, lwd = 2,
       legend = c("shallow", "mid", "deep"), bty = "n")
rect(580, -0.0001/5, 1020, 0.0001)
# dev.off()

# unconstrained model
ind = which(xseq < 1000)
# pdf("./narwal_case_study/joint_figure/nonparametric.pdf", width = 6, height = 5)
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.002), main = "unconstrained", 
     xlab = "max depth", ylab = "density", xlim = c(0, 1000))
for(i in 1:3) lines(xseq[ind], dens[ind,i], col = color[ord[i]], lwd = 2)
# lines(xseq, rowSums(dens), lty = 2, lwd = 1.5)
rect(580, -0.0001/5, 1020, 0.0001)
# dev.off()

# unimodal model
# pdf("./narwal_case_study/joint_figure/unimodal.pdf", width = 6, height = 5)
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.002), main = "constrained", 
     xlab = "max depth", ylab = "density", xlim = c(0, 1000))
for(i in 1:3) lines(xseq[ind], dens_c[ind,i], col = color[i], lwd = 2)
# lines(xseq, rowSums(dens_c), lty = 2, lwd = 1.5)
rect(580, -0.0001/5, 1020, 0.0001)
# dev.off()


## zoomed in

# gamma model
# pdf("./narwal_case_study/joint_figure/parametric_zoom.pdf", width = 6, height = 5)
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white",
     ylim = c(0, 0.0001), main = "parametric", 
     xlab = "max depth", ylab = "density", xlim = c(600, 1000), yaxt = "n")
box()
for(i in 1:N){
  curve(delta_par[i] * dgamma2(x, mu[i], sigma[i]), add = TRUE, col = color[i], 
        lwd = 2, from = 550, to = 1000)
} 
# curve(delta_par[1] * dgamma2(x, mu[1], sigma[1]) +
#         delta_par[2] * dgamma2(x, mu[2], sigma[2]) +
#         delta_par[3] * dgamma2(x, mu[3], sigma[3]), add = TRUE, lwd = 2, lty = 2)
axis(2, at = c(0, 0.00005, 0.0001), labels = c(0, 0.00005, 0.0001))
# dev.off()

# unconstrained model
# pdf("./narwal_case_study/joint_figure/nonparametric_zoom.pdf", width = 6, height = 5)
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.0001), main = "unconstrained", 
     xlab = "max depth", ylab = "density", xlim = c(600, 1000), yaxt = "n")
box()
for(i in 1:3) lines(xseq[ind], dens[ind,i], col = color[ord[i]], lwd = 2)
axis(2, at = c(0, 0.00005, 0.0001), labels = c(0, 0.00005, 0.0001))
# lines(xseq, rowSums(dens), lty = 2, lwd = 1.5)
# dev.off()

# unimodal model
# pdf("./narwal_case_study/joint_figure/unimodal_zoom.pdf", width = 6, height = 5)
hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     ylim = c(0, 0.0001), main = "constrained", 
     xlab = "max depth", ylab = "density", xlim = c(600, 1000), yaxt = "n")
box()
for(i in 1:3) lines(xseq[ind], dens_c[ind,i], col = color[i], lwd = 2)
axis(2, at = c(0, 0.00005, 0.0001), labels = c(0, 0.00005, 0.0001))
# lines(xseq, rowSums(dens_c), lty = 2, lwd = 1.5)
# dev.off()




# Decoded time series -----------------------------------------------------

mod_par$states = viterbi(mod = mod_par)
mod$states = viterbi(mod = mod)
mod_c$states = viterbi(mod = mod_c)

# idx = 2750:3000
idx = 2800 + 1:200
  
par(mfrow = c(1,1))

ylab = "max depth"
# pdf("./narwal_case_study/joint_figure/parametric_states.pdf", width = 6, height = 5)
plot(data$maxdep[idx], type = "h", col = color[mod_par$states[idx]], ylab = ylab, bty = "n")
# dev.off()

# pdf("./narwal_case_study/joint_figure/nonparametric_states.pdf", width = 6, height = 5)
plot(data$maxdep[idx], type = "h", col = color[ord[mod$states[idx]]], ylab = ylab, bty = "n") 
# dev.off()

# pdf("./narwal_case_study/joint_figure/unimodal_states.pdf", width = 6, height = 5)
plot(data$maxdep[idx], type = "h", col = color[mod_c$states[idx]], ylab = ylab, bty = "n")
# dev.off()


