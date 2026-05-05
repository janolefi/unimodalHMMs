## packages and utility functions
library(LaMa)
source("./utils.R")

# faster estimation speeds
TapeConfig(matmul = "plain")


## loading the narwhal data
data = readRDS("~/Downloads/MaxDepth_10_TimeSeries.RData")
data = subset(data, ID == "17001")
data$maxdep = data$maxdep - 10 # subtracting the depth threshold to make it easier for gamma model
nrow(data); range(data$maxdep, na.rm = TRUE)


## exploring potential initial values
hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", ylim = c(0, 0.002),
     main = "Distribution of maximum depth", xlab = "Max depth")
curve(0.75 * dgamma2(x, 20, 30), add = TRUE, col = color[1], lwd = 2)
curve(0.15 * dgamma2(x, 350, 70), add = TRUE, col = color[2], lwd = 2)
curve(0.1 * dgamma2(x, 500, 70), add = TRUE, col = color[3], lwd = 2)



# Parametric gamma model -------------------------------------------------

# negative log likelihood function for model with state-depenent gamma distributions
nll = function(par){
  getAll(par, dat) # makes everything accessible w/o $
  
  ## state process
  Gamma = tpm(eta) # mlogit for transitiom probabilities
  delta = stationary(Gamma) # stationary distribution of the state process
  
  ## state-dependent process
  # parameter transformations
  mu = exp(logmu); REPORT(mu) # reporting to make accessible later
  sigma = exp(logsigma); REPORT(sigma)
  # state-dependent densities
  allprobs = matrix(1, length(maxdep), N)
  ind = which(!is.na(maxdep)) # NA handling
  for(j in 1:N) allprobs[ind,j] = dgamma2(maxdep[ind], mu[j], sigma[j])
  
  ## forward algorithm
  -forward(delta, Gamma, allprobs)
}


## fitting the model

# Define base values (based on initial exploration)
eta_base = rep(-2, 6)
logmu_base = log(c(15, 350, 500))
logsigma_base = log(c(30, 70, 70))

# number of random initial values
n_init = 500

# Seed for RNG
set.seed(123)

# Generate n_init sets of random initial values
init_values <- lapply(1:n_init, function(i) {
  list(
    eta = rnorm(6, eta_base, 0.5),
    logmu = rnorm(3, logmu_base, 0.5),
    logsigma = rnorm(3, logsigma_base, 0.5)
  )
})

# data and hyperparamter object
dat = list(maxdep = data$maxdep, N = 3)

# fitting models with different initial values
llks = rep(NA, n_init)
pb = txtProgressBar(min = 0, max = n_init, style = 3)
for(i in 1:n_init){
  setTxtProgressBar(pb, i)
  obj = MakeADFun(nll, init_values[[i]], silent = TRUE) # creating automatically diff'able objective
  opt = tryCatch(
    nlminb(obj$par, obj$fn, obj$gr), # optimising it
    error = function(e) NULL) # handle errors
  if(!is.null(opt$objective)) llks[i] = -opt$objective # saving log-likelihood
  gc() # cleaning up
}
close(pb)

# plotting log likelihoods
plot(llks)

# refit best model
optpar = init_values[[which.max(llks)]]
saveRDS(optpar, "./narwal_case_study/objects/optpar_norm_10m.rds") # saving best initial values
optpar = readRDS("./narwal_case_study/objects/optpar_norm_10m.rds") # reading best initial values
obj = MakeADFun(nll, optpar)
opt = nlminb(obj$par, obj$fn, obj$gr)

# extract estimated parameters
mod_par = obj$report()

# state decoding
mod_par$states = viterbi(mod = mod_par)

# save estimated model
saveRDS(mod_par, "./narwal_case_study/mod_normal_10m.rds")
# mod_par = readRDS("./narwal_case_study/objects/mod_normal.rds")


###### Nonparametric fits -------------------------------------------------

# negative log likelihood function for the unconstrained AND constrained model 
# for the unconstrained model, we kappa = 0
pnll = function(par){
  getAll(par, dat)
  
  ## state process model
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  
  ## state-dependent process model
  # multinomial logit for density weights
  alpha = exp(cbind(beta, 0)) # last column is fixed at zero
  alpha = alpha / rowSums(alpha); REPORT(alpha) # reporting for easy access
  # state-dependent densities 
  allprobs = matrix(1, nrow(Z), N)
  ind = which(!is.na(Z[,1])) # NA handling
  allprobs[ind,] = Z[ind,] %*% t(alpha) # Z %*% alpha[i,] for all states at once
  
  # forward algorithm + smoothness penalty + potential unimodality penalty (for kappa > 0)
  res = -forward(delta, Gamma, log(allprobs), logspace = TRUE) + # negative log likelihood
    penalty(beta, S, lambda) # P-spline smoothness penalty
  
  if(kappa > 0) res = res + penalty_uni(beta + logweights, m, kappa) # optional unimodality penalty
  
  res
}


# Unconstrained nonparametric model ---------------------------------------

## need to do this for random initial values
# number of random initial values
n_init = 200

logmu_base = log(c(15, 350, 500))
logsigma_base = log(c(30, 70, 70))

# set.seed(123)
set.seed(12)

# Generate n_random random sets
init_values <- lapply(1:n_init, function(i) {
  list(
    eta = rnorm(6, eta_base, 0.2),
    logmu = rnorm(3, logmu_base, 0.2),
    logsigma = rnorm(3, logsigma_base, 0.3)
  )
})

## Setting up the P-spline basis and penalty
k = 30 # number of basis functions
degree = 3 # degree of the B-spline basis
# custom knots (including boundary knots)
knots = seq(0, max(data$maxdep, na.rm = TRUE), length = k - degree + 1) * 1.1

llks = rep(NA, n_init)
lambdas = list()
for(i in 1:n_init){
  cat("\n\nInitial value set", i, "\n")
  
  # constructing the smooth density object
  par = init_values[[i]]
  sDens = smooth_dens_construct(data["maxdep"],
                                list(maxdep = list(mean = exp(par$logmu), sd = exp(par$logsigma))),
                                type = "positive",
                                k = k,
                                knots = list(maxdep = knots))
  
  Z = sDens$Z$maxdep           # B-spline design matrix for the state-dependent densities
  S = sDens$S$maxdep           # P-spline penalty matrix for the state-dependent densities
  Z_p = sDens$Z_predict$maxdep # prediction design matrix for the state-dependent densities
  xseq = sDens$xseq$maxdep     # prediction grid for the state-dependent densities
  beta = sDens$coef$maxdep     # initial coefficients based on mean and sd
  
  # initial parameter list
  par = list(eta = par$eta, beta = beta) 
  
  # data and hyperparameter list
  dat = list(Z = Z, # B-spline desing matrix
             S = S, # P-spline penalty matrix
             N = 3, # number of states
             lambda = runif(3, 5, 20), # initial smoothenss penalty parameters
             kappa = 0) # no unimodality penalty here
  
  ## fitting the model using the qREML / extended Fellner-Schall method
  cat("\n")
  mod = tryCatch(
    qreml(pnll, par, dat, 
          random = "beta", # telling the model what are spline coefficients/ random effects
          conv_crit = "relchange", alpha = 0.3, silent = 1,
          maxiter = 50), # if it didn't converge in 35 interations, bad starting values.
    error = function(e) NULL, warning = function(w) NULL)
  
  if(!is.null(mod)) {
    llks[i] = mod$llk # save log likelihood
    lambdas[[i]] = mod$lambda
  }
  gc() # cleaning up
}

# plotting log likelihoods
plot(llks)


## refitting the best model
optpar = init_values[[which.max(llks)]]
optlambda = lambdas[[which.max(llks)]]

saveRDS(optpar, "./narwal_case_study/objects/optpar_np_10m.rds") # saving best initial values
saveRDS(optlambda, "./narwal_case_study/objects/optlambda_np_10m.rds") # saving best initial lambda

optpar = readRDS("./narwal_case_study/objects/optpar_np_10m.rds") # reading best initial values
optlambda = readRDS("./narwal_case_study/objects/optlambda_np_10m.rds") # reading best initial lambda

# constructing the smooth density object
sDens = smooth_dens_construct(data["maxdep"],
                              list(maxdep = list(mean = exp(optpar$logmu), 
                                                 sd = exp(optpar$logsigma))),
                              type = "positive",
                              k = k,
                              knots = list(maxdep = knots))

Z = sDens$Z$maxdep           # B-spline design matrix for the state-dependent densities
S = sDens$S$maxdep           # P-spline penalty matrix for the state-dependent densities
Z_p = sDens$Z_predict$maxdep # prediction design matrix for the state-dependent densities
xseq = sDens$xseq$maxdep     # prediction grid for the state-dependent densities
beta = sDens$coef$maxdep     # initial coefficients based on mean and sd

# initial parameter list
par = list(eta = optpar$eta, beta = beta) 

# data and hyperparameter list
dat = list(Z = Z, # B-spline desing matrix
           S = S, # P-spline penalty matrix
           N = 3, # number of states
           lambda = optlambda, # initial smoothenss penalty parameters
           kappa = 0) # no unimodality penalty

## fitting the best model
system.time(
  mod_np <- qreml(pnll, par, dat,
                  random = "beta",
                  conv_crit = "relchange")
)


## state decoding
mod_np$states = viterbi(mod = mod_np)

# plotting the estimated state-dependent densities
dens_np = Z_p %*% t(mod_np$alpha) # state-dependent densities
mod_np$dens = dens_np # saving the densities in the model object
mod_np$xseq = xseq # saving the prediction grid
delta_np = mod_np$delta # stationary distribution of the state process

saveRDS(mod_np, "./narwal_case_study/objects/mod_unconstrained_opt_10m.rds") # saving the model")
mod_np = readRDS("./narwal_case_study/objects/mod_unconstrained_opt_10m.rds") # reading the model


hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     xlim = c(0, 1000), ylim = c(0, 0.002),
     main = "Distribution of maximum depth", xlab = "Max depth")
for(i in 1:3) lines(mod_np$xseq, mod_np$delta[i] * mod_np$dens[,i], lwd = 2, col = color[i])

hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", 
     xlim = c(600, 1000), ylim = c(0, 0.0002),
     main = "Distribution of maximum depth", xlab = "Max depth")
for(i in 1:3) lines(mod_np$xseq, mod_np$delta[i] * mod_np$dens[,i], lwd = 2, col = color[i])


# Constrained (unimodal) nonparametric model ------------------------------

## defining potential coefficient modes for states 1:3
# state 1:
plot(mod_np$alpha[1,], pch = 16, col = "lightgray")
text(1:k, mod_np$alpha[1, ], labels = 1:k, cex = 0.9)
m1 = 1 # first mode definitely at one

# state 2:
plot(mod_np$alpha[2,], pch = 16, col = "lightgray")
text(1:k, mod_np$alpha[2, ], labels = 1:k, cex = 0.9)
# label switching, this is state 3
m2 = 11 + c(-2:2)

# state 3:
plot(mod_np$alpha[3,], pch = 16, col = "lightgray")
text(1:k, mod_np$alpha[3, ], labels = 1:k, cex = 0.9)
m3 = 14 + c(-2:2)

m_grid = expand.grid(m1 = m1, m2 = m2, m3 = m3)

## initial values
par = mod_np$par # use estimated parameters from the unconstrained model

## use data object from np fit
# add unimodality sharpness/ strength hyperparameter kappa
dat$kappa = 1e4

# updating initial smoothness parameter lambda
dat$lambda = mod_np$lambda * 3 # increase it a bit to avoid overfitting

# adding logweights for unimodality penalty
w = sDens$basis$maxdep$w # 1 / Integral of basis functions -> needed for reweighting in penalty
dat$logweights = t(log(w[-length(w)]) - log(w[length(w)]))[rep(1, 3), ]


## fitting models for all potential mode positions
llks = rep(NA, nrow(m_grid))
for(i in 1:nrow(m_grid)){
  cat("\nPosition", i, "of", nrow(m_grid), "\n")
  
  dat$m = m_grid[i, ] # setting coefficient modes for unimodality penalty
  mod = tryCatch(
    qreml(pnll, par, dat, 
          random = "beta", 
          alpha = 0.4,
          conv_crit = "relchange"), 
    error = function(e) NULL, warning = function(e) NULL)
  
  if(!is.null(mod)) llks[i] = mod$llk # saving log likelihood value
  gc() # cleaning up
}

plot(llks)

m_opt = m_grid[which.max(llks), ] # best modes
dat$m = m_opt
mod_uni = qreml(pnll, par, dat, 
                random = "beta",
                alpha = 0.4,
                conv_crit = "relchange")

## state decoding
mod_uni$states = viterbi(mod = mod_uni)

dens_uni = Z_p %*% t(mod_uni$alpha)
mod_uni$dens = dens_uni
mod_uni$xseq = xseq
delta_uni = mod_uni$delta

# save the model
saveRDS(mod_uni, "./narwal_case_study/objects/mod_constrained_opt_10m.rds")
mod_uni = readRDS("./narwal_case_study/objects/mod_constrained_opt_10m.rds") # reading the model

hist(data$maxdep, breaks = 30, prob = TRUE, bor = "white", 
     xlim = c(0, 1000), ylim = c(0, 0.002),
     main = "Distribution of maximum depth", xlab = "Max depth")
for(i in 1:3) lines(mod_uni$xseq, mod_uni$delta[i] * mod_uni$dens[,i], lwd = 2, col = color[i])

hist(data$maxdep, breaks = 50, prob = TRUE, bor = "white", 
     xlim = c(600, 1000), ylim = c(0, 0.0002),
     main = "Distribution of maximum depth", xlab = "Max depth")
for(i in 1:3) lines(mod_uni$xseq, mod_uni$delta[i] * mod_uni$dens[,i], lwd = 2, col = color[i])



# Plotting big state-dependent density figure -----------------------------

# pdf("./narwal_case_study/figures/narwhal_statedep_10m.pdf", width = 6, height = 5.5)

lwd = 1.5
ltys <- c(1,2,4)
breaks = 30

m = matrix(c(1,1, 2:7), 4, 2, byrow = TRUE)

layout(m, widths = c(1.3, 1), heights = c(0.06, 1, 1, 1.2))

### state legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("shallow", "mid", "deep"), lwd = lwd, col = color,
       horiz = TRUE, bty = "n", xpd = NA, cex = 1, lty = ltys)

### parametric model - zoomed out
par(mar = c(2,5.5,4,1)+0.1)
hist(data$maxdep + 20, breaks = breaks, prob = TRUE, bor = "white", 
     xlim = c(0, 1000), ylim = c(0, 0.002), las = 1, ylab = "",
     main = "", xaxt = "n")
axis(1, at = c(20, 200, 400, 600, 800, 1000), labels = c(20, 200, 400, 600, 800, 1000))
mtext("(a) Parametric model", side = 3, adj = 0, line = 1, font = 1, cex = 0.8)
mtext("Density", side = 2, line = 4.2, cex = 0.7)  # move y-axis label further from axis

for(i in 1:3) curve(mod_par$delta[i] * dgamma2(x, mod_par$mu[i], mod_par$sigma[i]), 
                    col = color[i], lwd = lwd, add = TRUE, from = 0, to = 1000, 
                    lty = ltys[i], n = 500)
rect(580, -0.0001/5, 1020, 0.0001)

### parametric model - zoomed in
par(mar = c(2,1,4,1)+0.1)
hist(data$maxdep + 20, breaks = breaks, prob = TRUE, bor = "white", 
     xlim = c(600, 1000), ylim = c(0, 0.00007), yaxt = "n",
     main = "", ylab = "")
mtext("(b) Parametric model - tail", side = 3, adj = 0, line = 1, font = 1, cex = 0.8)
for(i in 1:3) curve(mod_par$delta[i] * dgamma2(x, mod_par$mu[i], mod_par$sigma[i]), 
                    col = color[i], lwd = lwd, add = TRUE, from = 550, to = 1000, 
                    lty = ltys[i])
box()


### nonparametric model - zoomed out
par(mar = c(2,5.5,4,1)+0.1)
hist(data$maxdep + 20, breaks = breaks, prob = TRUE, bor = "white", 
     xlim = c(0, 1000), ylim = c(0, 0.002), las = 1, ylab = "",
     main = "", xaxt = "n")
axis(1, at = c(20, 200, 400, 600, 800, 1000), labels = c(20, 200, 400, 600, 800, 1000))

mtext("(c) Nonparametric model", side = 3, adj = 0, line = 1, font = 1, cex = 0.8)
mtext("Density", side = 2, line = 4.2, cex = 0.7)  # move y-axis label further from axis

idx = which(mod_np$xseq <= 1000)
for(i in 1:3) lines(mod_np$xseq[idx], mod_np$delta[i] * mod_np$dens[idx,i], 
                    lwd = lwd, col = color[i], lty = ltys[i])
rect(580, -0.0001/5, 1020, 0.0001)

### nonparametric model - zoomed in
par(mar = c(2,1,4,1)+0.1)
hist(data$maxdep + 20, breaks = breaks, prob = TRUE, bor = "white", 
     xlim = c(600, 1000), ylim = c(0, 0.00007), yaxt = "n",
     main = "", ylab = "")
mtext("(d) Nonparametric model - tail", side = 3, adj = 0, line = 1, font = 1, cex = 0.8)
for(i in 1:3) lines(mod_np$xseq[idx], mod_np$delta[i] * mod_np$dens[idx,i], 
                    lwd = lwd, col = color[i], lty = ltys[i])
box()


### unimodal model - zoomed out
par(mar = c(5,5.5,4,1)+0.1)
hist(data$maxdep + 20, breaks = breaks, prob = TRUE, bor = "white", 
     xlim = c(0, 1000), ylim = c(0, 0.002), las = 1, ylab = "",
     main = "", xlab = "Max depth (m)", xaxt = "n")
axis(1, at = c(20, 200, 400, 600, 800, 1000), labels = c(20, 200, 400, 600, 800, 1000))

mtext("(e) Unimodal model", side = 3, adj = 0, line = 1, font = 1, cex = 0.8)
mtext("Density", side = 2, line = 4.2, cex = 0.7)  # move y-axis label further from axis

idx = which(mod_uni$xseq <= 1000)
for(i in 1:3) lines(mod_uni$xseq[idx], mod_uni$delta[i] * mod_uni$dens[idx,i], 
                    lwd = lwd, col = color[i], lty = ltys[i])
rect(580, -0.0001/5, 1020, 0.0001)

### unimodal model - zoomed in
par(mar = c(5,1,4,1)+0.1)
hist(data$maxdep + 20, breaks = breaks, prob = TRUE, bor = "white", 
     xlim = c(600, 1000), ylim = c(0, 0.00007), yaxt = "n",
     main = "", ylab = "", xlab = "Max depth (m)")
mtext("(f) Unimodal model - tail", side = 3, adj = 0, line = 1, font = 1, cex = 0.8)
for(i in 1:3) lines(mod_uni$xseq[idx], mod_uni$delta[i] * mod_uni$dens[idx,i], 
                    lwd = lwd, col = color[i], lty = ltys[i])
box()

# dev.off()



## plotting the decoded time series
# pdf("narwal_case_study/figures/narwhal_states_new_10m.pdf", width = 6, height = 4)

# Define a 4-row layout: top row for legend, then the 3 plots
layout(matrix(1:4, nrow = 4), heights = c(0.3, 1, 1, 1.3))  # top row is small

idx = 2800 + 1:200
ylab = "Max depth (m)"
pchs = c(20,17,15)

# define color and labels (ensure colors match your actual state coloring)
state_labels <- c("shallow", "mid", "deep")

# Legend panel (empty plot with legend only)
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = state_labels, lwd = 1, col = color,
       horiz = TRUE, bty = "n", xpd = NA, cex = 0.9, , pch = pchs,
       text.width = strwidth(state_labels))

# panel (a): parametric
par(mar = c(1.2, 4, 0.5, 2) + 0.1)
plot(idx, data$maxdep[idx]+20, type = "h", col = color[mod_par$states[idx]],
     ylab = ylab, bty = "n", las = 1, xaxt = "n", yaxt = "n")
points(idx, data$maxdep[idx]+20, col = color[mod_par$states[idx]], 
       cex = 0.5, pch = pchs[mod_par$states[idx]])
axis(2, at = c(20, 200, 400, 600, 800), labels = c(20, 200, 400, 600, 800))
mtext("(a) Parametric", side = 3, adj = 0.02, line = -1, cex = 0.7)

# panel (b): unconstrained
par(mar = c(1.2, 4, 0.5, 2) + 0.1)
plot(idx, data$maxdep[idx], type = "h", col = color[mod_np$states[idx]],
     ylab = ylab, bty = "n", las = 1, xaxt = "n")
points(idx, data$maxdep[idx]+20, col = color[mod_np$states[idx]], 
       cex = 0.5, pch = pchs[mod_np$states[idx]])
mtext("(b) Nonparametric", side = 3, adj = 0.02, line = -1, cex = 0.7)

# panel (c): constrained
par(mar = c(4, 4, 0.5, 2) + 0.1)
plot(idx, data$maxdep[idx], type = "h", col = color[mod_uni$states[idx]],
     ylab = ylab, bty = "n", las = 1, xlab = "Observation index")
points(idx, data$maxdep[idx]+20, col = color[mod_uni$states[idx]], 
       cex = 0.5, pch = pchs[mod_uni$states[idx]])
mtext("(c) Unimodal", side = 3, adj = 0.02, line = -1, cex = 0.7)

# dev.off()


## persistence of the different models
### number of state switches
length(rle(mod_par$states)$values)
length(rle(mod_np$states)$values) # drastically higher persistence
length(rle(mod_uni$states)$values)

### estimated tpms
round(mod_par$Gamma, 3)
round(mod_np$Gamma, 3)
round(mod_uni$Gamma, 3)

