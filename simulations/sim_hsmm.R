## simulating an HSMM and fitting an HMM

library(LaMa)
library(parallel)
library(scales)
library(RTMBdist)

source("./utils.R")

TapeConfig(matmul = "plain")


# first state: mixture of two Poisson distributions
# second state: regular geometric dwell times

# idea: first state yields small observations, mostly longer stays but sometimes very short ones
# hence, the spline density model might assign observations from state 1 to state 2 because of the state-process missmatch

x = 0:25
dmean1 = c(0.1, 15)
d1 = 0.7 * dpois(x, dmean1[1]) + 0.3 * dpois(x, dmean1[2])
p = 1/10
d2 = dgeom(x, p)

# expected dwell time in state 1
sum((x+1) * d1)

# pdf("./simulations/figures/hsmm_sim_dwelltime.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
plot(x+1, d1, type = "h", lwd = 2, bty = "n", col = color[1], ylim = c(0, 0.6),
     xlab = "dwell time", ylab = "probabilities", main = "bimodal")
plot(x+1, d2, type = "h", lwd = 2, bty = "n", col = color[2], ylim = c(0, 0.1),
     xlab = "dwell time", ylab = "probabilities", main = "geometric")

# dev.off()

# function that simulates the gamma HSMM
sim_data = function(nObs = 1e3, 
                    dwellmean1 = c(0.1, 15),
                    dwellmean2 = 10,
                    prob = c(0.7, 0.3), 
                    mu = c(1, 15), 
                    sigma = c(1, 4),
                    nu = c(1,5)){
  c = rep(1:2, nObs/5)
  s = c()
  #for(i in 1:length(c)){
  i = 1
  while(length(s) <= nObs){
    if(c[i] == 1){
      mix = sample(1:2, 1, prob = prob)
      dwelltime = rpois(1, dwellmean1[mix]) + 1 # shifted Poisson
    } else{
      dwelltime = rgeom(1, 1/dwellmean2) + 1 # shifted geometric with p*(1-p)^{x-1}
    }
    s = c(s, rep(c[i], dwelltime))
    i = i + 1
  }
  s = s[1:nObs]
  # x = rgamma2(nObs, mu[s], sigma[s])
  x = rgengamma(nObs, mu[s], sigma[s], nu[s])
  data.frame(x = x, s = s)
}

# gamma-HMM likelihood function
nll <- function(par) {
  getAll(par)
  Gamma <- tpm(eta)
  delta <- stationary(Gamma)
  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  allprobs <- matrix(1, length(data$x), 2)
  for(j in 1:2) {
    allprobs[,j] <- dgamma2(data$x, mu[j], sigma[j])
  }
  -forward(delta, Gamma, allprobs)
}

# penalised likelihood function
pnll = function(par){
  getAll(par,dat)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  
  alpha = exp(cbind(beta, 0))
  alpha = alpha / rowSums(alpha); REPORT(alpha)
  
  ind = which(!is.na(Z[,1]))
  allprobs = matrix(1, nrow(Z), nrow(alpha))
  allprobs[ind, ] = Z[ind,] %*% t(alpha)
  
  res = -forward(delta, Gamma, allprobs) + penalty(beta, S, lambda)
  if(kappa > 0) {
    res = res + penalty_uni(beta + logweights, m, kappa, rho = 40)
  }
  
  res
}

# fit both models for one random initial parameter set
one_fit = function(dummy, data, mu, sigma){
  
  # random state-dependent mean and sd initials
  thismu = c(runif(1, 0, 5), runif(1, 6, 16))
  thissigma = exp(log(sigma) + rnorm(2, 0, 0.5))
  
  
  # fit correct parametric model
  # environment(nll) <- environment()
  par0 <- list(
    eta = rep(-2, 2) + rnorm(2, 0, 0.5),
    log_mu = log(thismu),
    log_sigma = log(thissigma)
  )
  obj_par <- MakeADFun(nll, par0, silent = TRUE)
  opt_par <- nlminb(obj_par$par, obj_par$fn, obj_par$gr)
  mod_par <- report(obj_par)
  
  # setting up the basis
  k = 30 # number of basis functions
  sDens = smooth_dens_construct(
    data = data["x"],
    par = list(x = list(mean = thismu, sd = thissigma)),
    type = "positive",
    k = k,
    knots = list(x = seq(0, max(data$x) * 1.05, length.out = k - 2))
  )
  
  Z = sDens$Z$x
  S = sDens$S$x
  Z_p = sDens$Z_predict$x
  xseq = sDens$xseq$x
  w = sDens$basis$x$w
  
  # fitting the unconstrained one
  par = list(
    eta = rep(-2, 2) + rnorm(2, 0, 0.5),
    beta = sDens$coef$x
  )
  
  dat = list(
    Z = Z, S = S, lambda = rep(200, 2), kappa = 0
  )
  
  mod = qreml(pnll, par, dat, 
              random = "beta", silent = 1, maxiter = 50)
  
  # computing the weighted state-dependent densities
  mod$dens = t(mod$delta * t(Z_p %*% t(mod$alpha))) # all (weighted by delta) at once
  mod$xseq = xseq
  
  ## fitting the constrained one
  m1 = 1
  m2 = 15:21
  m_grid = expand.grid(m1 = m1, m2 = m2)
  
  llks = rep(NA, nrow(m_grid))
  mods_c = list()
  for(i in 1:nrow(m_grid)){
    dat = list(
      Z = Z, S = S, lambda = rep(100, 2),
      m = m_grid[i, ],
      kappa = 250,
      logweights = t(log(w[-length(w)]) - log(w[length(w)]))[c(1,1), ]
    )
    
    mods_c[[i]] = qreml(pnll, par, dat, 
                        random = "beta", silent = 1, maxiter = 50)
    llks[i] = mods_c[[i]]$llk
  }
  
  mod_c = mods_c[[which.max(llks)]]
  
  # computing the weighted state-dependent densities
  mod_c$dens = t(mod_c$delta * t(Z_p %*% t(mod_c$alpha))) # all (weighted by delta) at once
  mod_c$xseq = xseq
  
  excl = c("Hessian_conditional", "obj_joint", "outer_gr", "relist_par")
  
  ret = list(
    mod_par = mod_par,
    mod = mod[!names(mod) %in% excl],
    mod_c = mod_c[!names(mod_c) %in% excl]
  )
  
  # cleaning up
  rm(mod_par)
  rm(mod)
  rm(mod_c)
  rm(data)
  rm(sDens)
  gc()
  
  ret
}

# safe wrapper to avoid failing
one_fit_safe = function(dummy, data, mu, sigma){
  tryCatch(
    one_fit(dummy, data, mu, sigma),
    error = function(e) {
      message("Error in one_fit: ", conditionMessage(e))
      return(NULL)
    }
  )
}


# Actual simulation -------------------------------------------------------


# setting true parameter values
nObs = 1e3
dwellmean1 = c(0.1, 15)
dwellmean2 = 10
prob = c(0.7, 0.3)
mu = c(2, 6)
sigma = c(1, 0.25)
nu = c(3, 8)

std <- sigma / mu

par(mfrow = c(1,1))
curve(dgengamma(x, mu[1], sigma[1], nu[1]), xlim = c(0.1, 12), bty = "n", lwd = 2)
curve(dgengamma(x, mu[2], sigma[2], nu[2]), add = TRUE, lwd = 2)

n_sim = 200
n_init = 10
results = list()

set.seed(123)

for(i in 1:n_sim){
  message(paste("Simulation run", i, "of", n_sim))
  
  data <- sim_data(nObs, dwellmean1, dwellmean2, prob, mu, sigma, nu)
  
  # parallelise over initial values
  thisres <- mclapply(1:n_init, one_fit_safe,
                      data = data,
                      mu = mu,
                      sigma = c(2, 3), 
                      mc.cores = 4)

  # extract likelihoods safely
  llk_par = lapply(thisres, function(x) x$mod_par$ll)
  llk_par = lapply(llk_par, function(x) if(is.null(x)) NA else x)
  llk_par = unlist(llk_par)
  llk = lapply(thisres, function(x) x$mod$llk)
  llk = lapply(llk, function(x) if(is.null(x)) NA else x)
  llk = unlist(llk)
  llk_c = lapply(thisres, function(x) x$mod_c$llk)
  llk_c = lapply(llk_c, function(x) if(is.null(x)) NA else x)
  llk_c = unlist(llk_c)

  mod_par = thisres[[which.max(llk_par)]]$mod_par
  mod = thisres[[which.max(llk)]]$mod
  mod_c = thisres[[which.max(llk_c)]]$mod_c
  
  ## compute state probabilities for all models
  mod_par$stateprobs = stateprobs(mod = mod_par)
  mod$stateprobs = stateprobs(mod = mod)
  mod_c$stateprobs = stateprobs(mod = mod_c)
  
  ## compute AUC for all models
  mod_par$auc = auc(data$s, mod_par$stateprobs)
  mod$auc = auc(data$s, mod$stateprobs)
  mod_c$auc = auc(data$s, mod_c$stateprobs)
  
  # number of misclassified observations
  states_par = ifelse(mod_par$stateprobs[,1] > 0.5, 1, 2)
  mod_par$misc = sum(states_par != data$s)
  states = ifelse(mod$stateprobs[,1] > 0.5, 1, 2)
  mod$misc = sum(states != data$s)
  states_c = ifelse(mod_c$stateprobs[,1] > 0.5, 1, 2)
  mod_c$misc = sum(states_c != data$s)
  
  results[[i]] = list(mod_par = mod_par, mod = mod, mod_c = mod_c)
  
  # cleaning up
  rm(data)
  rm(thisres)
  rm(mod_par)
  rm(mod)
  rm(mod_c)
  rm(llk_par)
  rm(llk)
  rm(llk_c)
  gc()
}


## Save results
# saveRDS(results, file = "./simulations/results/hsmm_sim_gengamma3.rds")
# results = readRDS("./simulations/results/hsmm_sim_gengamma3.rds")

# simulate data to get empirical approximation of true state distribution
set.seed(123)
bigdata = sim_data(1e5, dwellmean1, dwellmean2, prob, mu, sigma)
delta_hat = prop.table(table(bigdata$s)) # calculate state proportions


## plotting results

# pdf("./simulations/figures/hsmm_sim_2.pdf", width = 7, height = 4.5)

m = matrix(c(1,2,3, 4, 4, 4), nrow = 3, ncol = 2)
layout(mat = m, widths = c(1, 1), heights = c(1,1,1.2))

# parametric
par(mar = c(3,4,2,2))
plot(NA, xlim = c(0, 12), ylim = c(0,0.3), las = 1,
     xlab = "", ylab = "Density", bty = "n", main = "")
mtext("(a) Parametric", side = 3, adj = 0.5, line = 0.3, cex = 0.9)
for(i in 1:length(results)){
  thisres = results[[i]]
  
  this_mu <- thisres$mod_par$mu
  ord <- order(this_mu)
  this_sigma <- thisres$mod_par$sigma
  this_delta <- thisres$mod_par$delta
  
  for(j in 1:2){
    curve(this_delta[ord[j]] * dgamma2(x, this_mu[ord[j]], this_sigma[ord[j]]),
          add = TRUE, col = alpha(color[j], 0.5), lwd = 0.5, n = 500)
  }
}
for(j in 1:2){
  curve(delta_hat[j] * dgengamma(x, mu[j], sigma[j], nu[j]), 
        add = TRUE, lwd = 2, n = 500, col = "#000000", from = 0.01, to = 12)
}
legend("topright", bty = "n",
       legend = c("Estimated", "True"), lwd = c(0.5, 2))

# unconstrained
par(mar = c(3,4,2,2))
plot(NA, xlim = c(0, 12), ylim = c(0,0.3), las = 1,
     xlab = "", ylab = "Density", bty = "n", main = "")
mtext("(b) Nonparametric", side = 3, adj = 0.5, line = 0.3, cex = 0.9)

for(i in 1:length(results)){
  thisres = results[[i]]
  
  dens = thisres$mod$dens
  xseq = thisres$mod$xseq
  ord <- order(colSums(dens))
  
  for(j in 1:2){
    lines(xseq, dens[, ord[j]], col = alpha(color[j], 0.5), lwd = 0.5)
  }
}
for(j in 1:2){
  curve(delta_hat[j] * dgengamma(x, mu[j], sigma[j], nu[j]), 
        add = TRUE, lwd = 2, n = 500, col = "#000000", from = 0.01, to = 12)
}
# legend("topright", bty = "n",
#        legend = c("Estimated", "True"), lwd = c(0.5, 2))

par(mar = c(4.5,4,2,2))
plot(NA, xlim = c(0, 12), ylim = c(0,0.3), las = 1,
     xlab = "Observations", ylab = "Density", bty = "n", main = "")
mtext("(c) Unimodal", side = 3, adj = 0.5, line = 0.3, cex = 0.9)

for(i in 1:length(results)){
  thisres = results[[i]]
  
  dens_c = thisres$mod_c$dens
  xseq = thisres$mod_c$xseq
  
  for(j in 1:2){
    lines(xseq, dens_c[,j], col = alpha(color[j], 0.5), lwd = 0.5)
  }
}
for(j in 1:2){
  curve(delta_hat[j] * dgengamma(x, mu[j], sigma[j], nu[j]), 
        add = TRUE, lwd = 2, n = 500, col = "#000000", from = 0.01, to = 12)
}


## extract AUCs

auc_par = sapply(results, function(x) x$mod_par$auc)
auc = sapply(results, function(x) x$mod$auc)
auc[auc < 0.5] <- 1-auc[auc < 0.5]
auc_c = sapply(results, function(x) x$mod_c$auc)

df = data.frame(AUC = c(auc_par, auc, auc_c), 
                model = c(rep("Parametric", length(auc_par)),
                          rep("Nonparametric", length(auc)), 
                          rep("Unimodal", length(auc_c))))
df$model = factor(df$model, levels = c("Parametric", "Nonparametric", "Unimodal"))

boxplot(AUC ~ model, data = df, main = "", xlab = "Model", yaxt = "n", ylim = c(0.95, 1),
        pch = 20, col = "gray95", lwd = 0.5, outcol = "#00000030", frame = FALSE)
ci_par <- median_CI(df$AUC[df$model == "Parametric"])
ci_np <- median_CI(df$AUC[df$model == "Nonparametric"])
ci_uni <- median_CI(df$AUC[df$model == "Unimodal"])

# add semi-transparent light gray rectangles for the median CIs
rect(xleft = 0.6, xright = 1.4, ybottom = ci_par[1], ytop = ci_par[2],
     col = adjustcolor("black", alpha.f = 0.2), border = NA)
rect(xleft = 1.6, xright = 2.4, ybottom = ci_np[1], ytop = ci_np[2],
     col = adjustcolor("black", alpha.f = 0.2), border = NA)
rect(xleft = 2.6, xright = 3.4, ybottom = ci_uni[1], ytop = ci_uni[2],
     col = adjustcolor("black", alpha.f = 0.2), border = NA)

mtext("(d) AUC", side = 3, adj = 0.5, line = 0.2, cex = 0.9)
axis(2, at = seq(0.95, 1, by = 0.01), labels = seq(0.95, 1, by = 0.01), las = 1)

# dev.off()


## extract estimated tpm entries

par(mfrow = c(1,2))
# parametric
gamma22_par = sapply(results, function(x) x$mod_par$Gamma[2,2])
round(mean(gamma22_par), 4)
round(sd(gamma22_par), 4)
# dwell time
round(mean(1 / (1 - gamma22_par)), 4)
round(sd(1 / (1 - gamma22_par)), 4)

hist(gamma22_par, xlim = c(0.85, 1), bor = "white")
abline(v = 0.9)

par(mfrow = c(1,2))
# unconstrained
gamma22 = sapply(results, function(x) x$mod$Gamma[2,2])
round(mean(gamma22), 4)
round(sd(gamma22), 4)
# dwell time
round(mean(1 / (1 - gamma22)), 4)
round(sd(1 / (1 - gamma22)), 4)

hist(gamma22, xlim = c(0.85, 1), bor = "white")
abline(v = 0.9)

# constrained
gamma22_c = unlist(sapply(results, function(x) x$mod_c$Gamma[2,2]))
round(mean(gamma22_c), 4)
round(sd(gamma22_c), 4)
# dwell time
round(mean(1 / (1 - gamma22_c)), 4)
round(sd(1 / (1 - gamma22_c)), 4)

hist(gamma22_c, xlim = c(0.85, 1), bor = "white")
abline(v = 0.9)



# percentage of erroneously decoded states
perc_wrong_par <- sapply(results, function(r) r$mod_par$misc / nObs)
perc_wrong_par[perc_wrong_par > 0.5] <- 1 - perc_wrong_par[perc_wrong_par > 0.5] # label switching
round(mean(perc_wrong_par), 3)
round(sd(perc_wrong_par), 3)

perc_wrong <- sapply(results, function(r) r$mod$misc / nObs)
perc_wrong[perc_wrong > 0.5] <- 1 - perc_wrong[perc_wrong > 0.5] # label switching
round(mean(perc_wrong), 3)
round(sd(perc_wrong), 3)

perc_wrong_c <- sapply(results, function(r) r$mod_c$misc / nObs)
round(mean(perc_wrong_c), 3)
round(sd(perc_wrong_c), 3)

                              