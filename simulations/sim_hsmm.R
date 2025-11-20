## simulating an HSMM and fitting an HMM

library(LaMa)
library(parallel)
library(scales)

source("./utils.R")


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
                    sigma = c(1, 4)){
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
  x = rgamma2(nObs, mu[s], sigma[s])
  data.frame(x = x, s = s)
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
  # cpen = 0 # initialise with zero
  # for(i in 1:2) {
  #   cpen = cpen - sum(min0_smooth(C[[i]] %*% (beta[i,] + logw), rho = 20))
  # }
  # p + kappa * cpen
  if(kappa > 0) res = res + penalty_uni(beta + logweights, m, kappa)
  
  res
}

# fit both models for one random initial parameter set
one_fit = function(dummy, data, mu, sigma){

  # random state-dependent mean and sd initials
  thismu = c(runif(1, 0, 5), runif(1, 10, 20))
  thissigma = exp(log(sigma) + rnorm(2, 0, 0.3))
  
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
    Z = Z, S = S, lambda = rep(30, 2), kappa = 0
  )
  
  mod = qreml(pnll, par, dat, 
              random = "beta", silent = 1, maxiter = 50)
  
  # computing the weighted state-dependent densities
  mod$dens = t(mod$delta * t(Z_p %*% t(mod$alpha))) # all (weighted by delta) at once
  mod$xseq = xseq
  
  ## fitting the constrained one
  m1 = 1
  m2 = 11:14
  m_grid = expand.grid(m1 = m1, m2 = m2)
  
  llks = rep(NA, nrow(m_grid))
  mods_c = list()
  for(i in 1:nrow(m_grid)){
    dat = list(
      Z = Z, S = S, lambda = rep(40, 2),
      # C = construct_C(m_grid[i,], 2, k),
      m = m_grid[i, ],
      kappa = 1e2,
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
    mod = mod[!names(mod) %in% excl],
    mod_c = mod_c[!names(mod_c) %in% excl]
  )
  
  # cleaning up
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
mu = c(1, 15)
sigma = c(1, 4)

n_sim = 100
n_init = 10
results = list()

for(i in 1:n_sim){
  cat("Simulation run", i, "of", n_sim, "\n")
  
  data = sim_data(nObs, dwellmean1, dwellmean2, prob, mu, sigma)
  
  # parallelise over initial values
  thisres = mclapply(1:n_init, one_fit_safe, 
                     data = data, 
                     mu = mu,
                     sigma = sigma, mc.cores = 4)
  
  # extract likelihoods safely
  llk = lapply(thisres, function(x) x$mod$llk)
  llk = lapply(llk, function(x) if(is.null(x)) NA else x)
  llk = unlist(llk)
  llk_c = lapply(thisres, function(x) x$mod_c$llk)
  llk_c = lapply(llk_c, function(x) if(is.null(x)) NA else x)
  llk_c = unlist(llk_c)
  
  mod = thisres[[which.max(llk)]]$mod
  mod_c = thisres[[which.max(llk_c)]]$mod_c
  
  ## compute state probabilities for all models
  mod$stateprobs = stateprobs(mod = mod)
  mod_c$stateprobs = stateprobs(mod = mod_c)
  
  ## compute AUC for all models
  mod$auc = auc(data$s, mod$stateprobs)
  mod_c$auc = auc(data$s, mod_c$stateprobs)
  
  results[[i]] = list(mod = mod, mod_c = mod_c)
  
  # cleaning up
  rm(data)
  rm(thisres)
  rm(mod)
  rm(mod_c)
  rm(llk)
  rm(llk_c)
  gc()
}

## Save results
saveRDS(results, file = "./simulations/results/hsmm_sim.rds")
results = readRDS("./simulations/results/hsmm_sim.rds")

# simulate data to get empirical approximation of true state distribution
set.seed(123)
bigdata = sim_data(1e5, dwellmean1, dwellmean2, prob, mu, sigma)
delta_hat = prop.table(table(bigdata$s))


## plotting results

# pdf("./simulations/figures/hsmm_sim_w_boxplot2.pdf", width = 7, height = 4.5)

m = matrix(c(1,2,3,3), nrow = 2, ncol = 2)
layout(mat = m, widths = c(1.5, 1.2))

# unconstrained
par(mar = c(4.5,4,2,2))
plot(NA, xlim = c(0, 25), ylim = c(0,0.15), las = 1,
     xlab = "Observations", ylab = "Density", bty = "n", main = "")
mtext("(a) Nonparametric", side = 3, adj = 0.5, line = 0.3, cex = 0.9)

for(i in 1:length(results)){
  thisres = results[[i]]

  dens = thisres$mod$dens
  xseq = thisres$mod$xseq
  
  for(j in 1:2){
    lines(xseq, dens[,j], col = alpha(color[j], 0.5), lwd = 0.5)
  }
}
for(j in 1:2){
  curve(delta_hat[j] * dgamma2(x, mu[j], sigma[j]), 
        add = TRUE, lwd = 2, n = 500, col = "#000000", to = 26)
}
legend("topright", bty = "n",
       legend = c("Estimated", "True"), lwd = c(0.5, 2))

par(mar = c(4.5,4,2,2))
plot(NA, xlim = c(0, 25), ylim = c(0,0.15), las = 1,
     xlab = "Observations", ylab = "Density", bty = "n", main = "")
mtext("(b) Unimodal", side = 3, adj = 0.5, line = 0.3, cex = 0.9)

for(i in 1:length(results)){
  thisres = results[[i]]
  
  dens_c = thisres$mod_c$dens
  xseq = thisres$mod_c$xseq
  
  for(j in 1:2){
    lines(xseq, dens_c[,j], col = alpha(color[j], 0.5), lwd = 0.5)
  }
}
for(j in 1:2){
  curve(delta_hat[j] * dgamma2(x, mu[j], sigma[j]), 
        add = TRUE, lwd = 2, n = 500, col = "#000000", to = 26)
}


## extract AUCs

auc = sapply(results, function(x) x$mod$auc)
auc_c = sapply(results, function(x) x$mod_c$auc)

df = data.frame(AUC = c(auc, auc_c), 
                model = c(rep("Nonparametric", length(auc)), 
                          rep("Unimodal", length(auc_c))))
df$model = factor(df$model, levels = c("Nonparametric", "Unimodal"))

boxplot(AUC ~ model, data = df, main = "", xlab = "Model", yaxt = "n", ylim = c(0.997, 1),
        pch = 20, col = "gray95", lwd = 0.5, outcol = "#00000030", frame = FALSE)
mtext("(c) AUC", side = 3, adj = 0.5, line = 0.2, cex = 0.9)
axis(2, at = c(0.997, 0.998, 0.999, 1), labels = c(0.997, 0.998, 0.999, 1), las = 1)

# dev.off()


## extract estimated tpm entries

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


