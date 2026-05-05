## simulating an HSMM and fitting an HMM

library(LaMa)
library(parallel)
library(scales)

TapeConfig(matmul = "plain")

source("./utils.R")

# first state: mixture of two Poisson distributions
# second state: regular geometric dwell times

# idea: first state yields small observations, mostly longer stays but sometimes very short ones
# hence, the spline density model might assign observations from state 1 to state 2 because of the state-process missmatch

# Simulate data
# 2-states: skewnormal and t
# simulating data function
sim_data <- function(nsamples, par, Gamma, N = 2){
  # state process
  delta <- stationary(Gamma)
  s <- rep(NA, nsamples)
  s[1] <- sample(1:N, 1, prob = delta)
  for(t in 2:nsamples){
    s[t] <- sample(1:N, 1, prob = Gamma[s[t-1],])
  }
  # observed process
  x <- rep(NA, nsamples)
  x[s == 1] <- rskewnorm(sum(s == 1), par$xi, par$omega, par$alpha)
  x[s == 2] <- rt(sum(s == 2), par$df) * par$sigma + par$mu
  return(data.frame(x = x, s = s))
}

## parametric likelihood functions

# parametric likelihood skewnormal/ t
nll_st <- function(par){
  getAll(par, dat)
  Gamma = tpm(beta0)
  delta = stationary(Gamma)
  omega = exp(logomega); REPORT(omega)
  REPORT(xi); REPORT(alpha)
  sd = exp(logsd); REPORT(sd)
  df = exp(logdf); REPORT(df)
  REPORT(mu)
  allprobs = matrix(1, length(x), 2)
  ind = which(!is.na(x))
  allprobs[ind,1] = dskewnorm(x[ind], xi, omega, alpha)
  allprobs[ind,2] = dt((x[ind]-mu)/sd, df) / sd
  -forward(delta, Gamma, allprobs)
}




# fit both models for one random initial parameter set
one_fit = function(dummy, data, truepar){
  
  library(LaMa)
  
  ## likelihood functions
  # parametric likelihood normal
  nll_norm <- function(par){
    getAll(par, dat)
    Gamma = tpm(beta0)
    delta = stationary(Gamma)
    sigma = exp(logsigma)
    REPORT(mu); REPORT(sigma)
    allprobs = matrix(1, length(x), 2)
    ind = which(!is.na(x))
    allprobs[ind,] = cbind(dnorm(x[ind], mu[1], sigma[1]),
                           dnorm(x[ind], mu[2], sigma[2]))
    -forward(delta, Gamma, allprobs)
  }
  environment(nll_norm) <- environment()
  
  # nonparametric likelihood function
  pnll = function(par){
    getAll(par, dat)
    # regular stationary HMM stuff
    Gamma = tpm(beta0)
    delta = stationary(Gamma)
    # smooth state-dependent densities
    alpha = exp(cbind(beta, 0))
    alpha = alpha / rowSums(alpha) # multinomial logit link
    REPORT(alpha)
    allprobs = matrix(1, nrow(Z), N)
    ind = which(!is.na(Z[,1])) # only for non-NA obs.
    allprobs[ind,] = Z[ind,] %*% t(alpha)
    # forward algorithm + P-spline penalty
    res = -forward(delta, Gamma, allprobs) + penalty(beta, S, lambda)
    # unimodality constraints formulated as penalty
    if(kappa > 0) res = res + penalty_uni(beta + logweights, m, kappa)
    res
  }
  
  ## Fit wrong parametric model
  # random initial values
  par = list(
    beta0 = runif(2, -4, -1),
    mu = truepar$mean + rnorm(2, 0, 0.5),
    logsigma = log(truepar$sd) + rnorm(2, 0, 0.5)
  )
  # observations
  dat = list(x = data$x)
  
  # model fitting and reporting
  obj_par_norm = MakeADFun(nll_norm, par)
  opt_par_norm = tryCatch(
    optim(obj_par_norm$par, obj_par_norm$fn, obj_par_norm$gr, method = "BFGS"),
    error = function(e) NULL)
  mod_par_norm = obj_par_norm$report()
  mod_par_norm$llk = -opt_par_norm$value
  
  
  # exlude this from saved models
  # excl = c("allprobs", "Hessian_conditional", "obj_joint", "outer_gr", "relist_par")
  excl = c("Hessian_conditional", "obj_joint", "outer_gr", "relist_par")
  
  par = list(
    beta0 = runif(2, -3, -1),
    mu = rnorm(2, c(truepar$xi, truepar$mu), 0.5),
    logsigma = rnorm(2, log(c(truepar$omega, truepar$sd)), 0.5)
  )
  
  ## Fit spline model
  k = 40
  # random initial values based on the ones for the parametric model
  par0_smooth = list(
    mean = par$mu,
    sd = exp(par$logsigma)
  )
  
  sDens = smooth_dens_construct(data["x"], k = k,
                                par = list(x = par0_smooth))
  Z = sDens$Z$x
  S = sDens$S$x
  Z_p = sDens$Z_predict$x
  xseq = sDens$xseq$x
  coef = sDens$coef$x
  w = sDens$basis$x$w
  
  par = list(beta0 = runif(2, -3, -1), 
             beta = coef)
  dat = list(Z = Z, S = S, lambda = rep(10, 2), N = 2, kappa = 0,
             # C = construct_C(c(1,1), 2, k), # doesn't matter here bc kappa = 0
             logweights = t(log(w[-length(w)]) - log(w[length(w)]))[c(1,1), ]) 
  
  ## unconstrained fit
  mod_np = qreml(pnll, par, dat, random = "beta", 
                 silent = 1, maxiter = 50)
  
  # computing the weighted state-dependent densities
  mod_np$dens = t(mod_np$delta * t(Z_p %*% t(mod_np$alpha)))
  mod_np$xseq = xseq
  mod_np = mod_np[!names(mod_np) %in% excl] # excluding big objects
  
  
  ## constrained fit
  # use the same initial parameter values as for the unconstrained fit
  
  ms = set_m_grid(mod_np$par$beta) # m_grid based on previous optimum
  # add indices I want to check for sure
  m1 = 12:13
  m2 = 16:17
  ms = rbind(ms, expand.grid(m1 = m1, m2 = m2))
  # drop duplicates
  ms = unique(ms)
  
  dat$lambda = rep(50, 2)
  dat$kappa = 1e3 * 5
  llks = rep(NA, nrow(ms))
  mods_c = list()
  for(i in 1:nrow(ms)){
    # thism = ms[i,]
    # dat$C = construct_C(thism, 2, k)
    dat$m = ms[i,]
    mods_c[[i]] = qreml(pnll, par, dat, random = "beta", 
                        silent = 1, maxiter = 30, alpha = 0.2)
    llks[i] =  mods_c[[i]]$llk
  }
  # optimal one
  mod_c =  mods_c[[which.max(llks)]]
  rm(mods_c)
  dens_c = t(mod_c$delta * t(Z_p %*% t(mod_c$alpha)))
  mod_c$dens = dens_c
  mod_c$xseq = xseq
  mod_c = mod_c[!names(mod_c) %in% excl]
  
  ## returning
  out = list(
    mod_par_norm = mod_par_norm,
    mod_np = mod_np,
    mod_c = mod_c
  )
  
  ## cleanup for memory
  rm(mod_np); rm(mod_c); rm(mod_par_norm); rm(obj_par_norm)
  rm(data)
  rm(sDens)
  rm(dat)
  rm(w)
  rm(Z_p); rm(Z); rm(S); rm(xseq); rm(coef)
  
  ## garbage collector
  gc()
  
  ## returning
  return(out)
}

# safe wrapper to avoid failing
one_fit_safe = function(dummy, data, truepar){
  tryCatch(
    one_fit(dummy, data, truepar),
    error = function(e) {
      message("Error in one_fit: ", conditionMessage(e))
      return(NULL)
    }
  )
}


# Actual simulation -------------------------------------------------------

set.seed(123)

# setting true parameter values for skew normal and t distributions
truepar = list(xi = 0, omega = 1, alpha = 6,
               mu = 3, sigma = 1, df = 3)

# adding true implied means and standard deviations
d = truepar$alpha / sqrt(1 + truepar$alpha^2)
truepar$mean = c(truepar$xi + truepar$omega * d * sqrt(2 / pi), # mean of skew normal
                 truepar$mu)
truepar$sd = c(sqrt(truepar$omega^2 * (1 - (2 * d^2 / pi))), # sd of skew normal
               truepar$sigma * sqrt((truepar$df - 2) / truepar$df)) # sd of t

# true transition probability matrix used for simulating the chains
Gamma = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)

n_sim = 100
n_init = 15
nObs = 500
results = list()

t1 = Sys.time()
i = 1
while(i <= n_sim){
  message(paste("Simulation run", i, "of", n_sim))
  
  data = sim_data(nObs, truepar, Gamma)
  
  ## fit parametric models -- one set of initial values: true parameters
  ## Fit correct parametric model
  # random initial values
  par = list(
    beta0 = qlogis(rep(0.1, 2)),
    xi = truepar$xi,
    logomega = log(truepar$omega),
    alpha = truepar$alpha,
    mu = truepar$mu,
    logsd = log(truepar$sigma),
    logdf = log(truepar$df)
  )
  # observations
  dat = list(x = data$x)
  
  # model fitting and reporting
  obj_par = MakeADFun(nll_st, par, silent = TRUE)
  opt_par = tryCatch(
    optim(obj_par$par, obj_par$fn, obj_par$gr, method = "BFGS"),
    error = function(e) NULL)
  if(!is.list(opt_par)){
    next
  }
  mod_par = obj_par$report()
  mod_par$llk = -opt_par$value
  
  
  ## fit nonparametric models
  # parallelise over initial values
  thisres = mclapply(1:n_init, one_fit_safe, 
                     data = data, 
                     truepar = truepar,
                     mc.cores = 4)
  
  # extract likelihoods to find best model
  llk_norm = lapply(thisres, function(x) x$mod_par_norm$llk)
  llk_norm = lapply(llk_norm, function(x) if(is.null(x)) NA else x)
  llk_norm = unlist(llk_norm)
  
  llk_np = lapply(thisres, function(x) x$mod_np$llk)
  llk_np = lapply(llk_np, function(x) if(is.null(x)) NA else x)
  llk_np = unlist(llk_np)
  
  llk_c = lapply(thisres, function(x) x$mod_c$llk)
  llk_c = lapply(llk_c, function(x) if(is.null(x)) NA else x)
  llk_c = unlist(llk_c)
  
  if(all(is.na(llk_norm))){
    mod_par_norm = NULL
  } else{
    mod_par_norm = thisres[[which.max(llk_norm)]]$mod_par_norm
  }
  
  if(all(is.na(llk_np))){
    mod_np = NULL
  } else{
    mod_np = thisres[[which.max(llk_np)]]$mod_np
  }
  
  if(all(is.na(llk_c))){
    mod_c = NULL
  } else{
    mod_c = thisres[[which.max(llk_c)]]$mod_c
  }
  
  # if any fit did not work at all, skip this iteration
  if(is.null(mod_par) || is.null(mod_par_norm) || is.null(mod_np) || is.null(mod_c)){
    next
  }
  
  ## compute state probabilities for all models
  mod_par$stateprobs = stateprobs(mod = mod_par)
  mod_par_norm$stateprobs = stateprobs(mod = mod_par_norm)
  mod_np$stateprobs = stateprobs(mod = mod_np)
  mod_c$stateprobs = stateprobs(mod = mod_c)
  
  ## compute viterbi sequence for all models
  mod_par$states = viterbi(mod = mod_par)
  mod_par_norm$states = viterbi(mod = mod_par_norm)
  mod_np$states = viterbi(mod = mod_np)
  mod_c$states = viterbi(mod = mod_c)
  
  ## compute AUC for all models
  mod_par$auc = auc(data$s, mod_par$stateprobs)
  mod_par_norm$auc = auc(data$s, mod_par_norm$stateprobs)
  mod_np$auc = auc(data$s, mod_np$stateprobs)
  mod_c$auc = auc(data$s, mod_c$stateprobs)
  
  results[[i]] = list(
    mod_par = mod_par,
    mod_par_norm = mod_par_norm,
    mod_np = mod_np,
    mod_c = mod_c
  )
  
  # cleanup
  rm(data); rm(dat); rm(par); rm(thisres)
  rm(mod_par); rm(mod_par_norm); rm(mod_np); rm(mod_c)
  rm(llk_np); rm(llk_c)
  
  # garbage collector
  gc()
  
  # increase counter
  i = i + 1
}
Sys.time() - t1
gc()


## Save results
# saveRDS(results, file = "./simulations/results/simple_sim100_new.rds")
# results = readRDS("./simulations/results/simple_sim100_new.rds")


## Filter results
results = results[sapply(results, length) > 1]
nresults = length(results)


## plotting results

# pdf("./simulations/figures/simple_sim_new.pdf", width = 8, height = 5)

# plotting the results
par(mfrow = c(2,2), mar = c(5,4,2,1) + 0.1)

# parametric true
plot(NA, xlim = c(-2, 8), ylim = c(0, 0.5), bty = "n", las = 1,
     ylab = "Density", xlab = "Observations")
mtext("(a) Skew normal - t", side = 3, adj = 0.5, line = 0.2, cex = 0.9)
for(i in 1:nresults){
  thismod = results[[i]]$mod_par
  delta = thismod$alpha / sqrt(1 + thismod$alpha^2)
  thism = thismod$xi + thismod$omega * delta * sqrt(2 / pi)
  thiscol = color[order(c(thism, thismod$mu))]
  
  curve(thismod$delta[1] * dskewnorm(x, thismod$xi, thismod$omega, thismod$alpha),
        lwd = 0.5, col = alpha(thiscol[1], 0.2), add = TRUE, n = 500)
  curve(thismod$delta[2] * dt((x-thismod$mu)/thismod$sd, thismod$df) / thismod$sd,
        add = TRUE, lwd = 0.5, col = alpha(thiscol[2], 0.2), n = 500)
}
curve(0.5 * dskewnorm(x, truepar$xi, truepar$omega, truepar$alpha),
      lwd = 2, col = "#000000", n = 500, add = TRUE)
curve(0.5 * dt((x-truepar$mu) / truepar$sigma, truepar$df) / truepar$sigma,
      lwd = 2, col = "#000000", n = 500, add = TRUE)


# parametric wrong
plot(NA, xlim = c(-2, 8), ylim = c(0, 0.5), bty = "n", las = 1,
     ylab = "Density", xlab = "Observations")
mtext("(b) Normal", side = 3, adj = 0.5, line = 0.2, cex = 0.9)
for(i in 1:nresults){
  thismod = results[[i]]$mod_par_norm
  thiscol = color[order(thismod$mu)]
  
  for(j in 1:2){
    curve(thismod$delta[j] * dnorm(x, thismod$mu[j], thismod$sigma[j]),
          lwd = 0.5, col = alpha(thiscol[j], 0.2), add = TRUE, n = 500)
  }
}
curve(0.5 * dskewnorm(x, truepar$xi, truepar$omega, truepar$alpha),
      lwd = 2, col = "#000000", n = 500, add = TRUE)
curve(0.5 * dt((x-truepar$mu) / truepar$sigma, truepar$df) / truepar$sigma,
      lwd = 2, col = "#000000", n = 500, add = TRUE)


# nonparametric unconstrained
label_switch = rep(FALSE, nresults)
plot(NA, xlim = c(-2, 8), ylim = c(0, 0.5), bty = "n", las = 1,
     ylab = "Density", xlab = "Observations")
mtext("(c) Nonparametric", side = 3, adj = 0.5, line = 0.2, cex = 0.9)
for(i in 1:nresults){
  thismod = results[[i]]$mod_np
  thisxseq = thismod$xseq
  thisdens = thismod$dens
  thisdens2 = t(thisdens) / rowSums(t(thisdens)) # turn into discrete pmf for xseq
  thismean = sapply(1:2, function(j) sum(thisdens2[j,] * thisxseq))
  label_switch[i] = thismean[1] > thismean[2]
  thiscol = color[order(thismean)]
  
  for(j in 1:2){
    lines(thisxseq, thisdens[,j], lwd = 0.5, col = alpha(thiscol[j], 0.2))
  }
}
curve(0.5 * dskewnorm(x, truepar$xi, truepar$omega, truepar$alpha),
      lwd = 2, col = "#000000", n = 500, add = TRUE)
curve(0.5 * dt((x-truepar$mu) / truepar$sigma, truepar$df) / truepar$sigma,
      lwd = 2, col = "#000000", n = 500, add = TRUE)

# nonparametric constrained
label_switch_c = rep(FALSE, nresults)
plot(NA, xlim = c(-2, 8), ylim = c(0, 0.5), bty = "n", las = 1,
     ylab = "Density", xlab = "Observations")
mtext("(d) Unimodal", side = 3, adj = 0.5, line = 0.2, cex = 0.9)
for(i in 1:nresults){
  thismod = results[[i]]$mod_c
  thisxseq = thismod$xseq
  thisdens = thismod$dens
  thisdens2 = t(thisdens) / rowSums(t(thisdens))
  thismean = sapply(1:2, function(j) sum(thisdens2[j,] * thisxseq))
  label_switch_c[i] = thismean[1] > thismean[2]
  thiscol = color[order(thismean)]
  
  for(j in 1:2){
    lines(thisxseq, thisdens[,j], lwd = 0.5, col = alpha(thiscol[j], 0.2))
  }
}
curve(0.5 * dskewnorm(x, truepar$xi, truepar$omega, truepar$alpha),
      lwd = 2, col = "#000000", n = 500, add = TRUE)
curve(0.5 * dt((x-truepar$mu) / truepar$sigma, truepar$df) / truepar$sigma,
      lwd = 2, col = "#000000", n = 500, add = TRUE)

# dev.off()


# compute mean AUC for all model classes

auc_par = sapply(results, function(x) x$mod_par$auc)
auc_par_norm = sapply(results, function(x) x$mod_par_norm$auc)
auc_par_norm[auc_par_norm < 0.5] <- 1-auc_par_norm[auc_par_norm < 0.5]
auc_np = sapply(results, function(x) x$mod_np$auc)
auc_c = sapply(results, function(x) x$mod_c$auc)
# correct for label switching
auc_c[label_switch_c] = 1 - auc_c[label_switch_c]
auc_c[auc_c < 0.5] = 1 - auc_c[auc_c < 0.5]

names = c("Skew normal - t", "Normal", "Nonparametric", "Unimodal")
AUC = data.frame(
  auc = c(auc_par, auc_par_norm, auc_np, auc_c),
  model = rep(names, each = nresults)
)
AUC$model = factor(AUC$mod, levels = names)

# compute CIs for median
ci_snt <- median_CI(AUC$auc[AUC$model == "Skew normal - t"])
ci_no <- median_CI(AUC$auc[AUC$model == "Normal"])
ci_np <- median_CI(AUC$auc[AUC$model == "Nonparametric"])
ci_uni <- median_CI(AUC$auc[AUC$model == "Unimodal"])


# pdf("./simulations/figures/simple_sim_boxplot_CI.pdf", width = 5, height = 3.5)

# par(mfrow = c(1,2))
par(mfrow = c(1,1), mar = c(5,7,3,2)+0.1)
# boxplot(auc ~ model, data = AUC, ylab = "AUC", xlab = "Model", main = "(a)",
#         pch = 20, col = "gray95", lwd = 0.5, outcol = "#00000030", frame = FALSE)
boxplot(auc ~ model, data = AUC, ylim = c(0.97, 1), ylab = "", xlab = "AUC", main = "",
        horizontal = TRUE,
        pch = 20, col = "gray95", lwd = 0.5, outcol = "#00000030", frame = FALSE, las = 1)
abline(v = median(AUC$auc[which(AUC$model == "Skew normal - t")]), col = color[2], lty = 2)
# mtext("AUC", side = 2, line = 4, cex = 1)  # move y-axis label further from axis

# add semi-transparent light gray rectangles for the median CIs
rect(xleft = ci_snt[1], xright = ci_snt[2], ybottom = 0.6, ytop = 1.4,
     col = adjustcolor("black", alpha.f = 0.2), border = NA)
rect(xleft = ci_no[1], xright = ci_no[2], ybottom = 1.6, ytop = 2.4,
     col = adjustcolor("black", alpha.f = 0.2), border = NA)
rect(xleft = ci_np[1], xright = ci_np[2], ybottom = 2.6, ytop = 3.4,
     col = adjustcolor("black", alpha.f = 0.2), border = NA)
rect(xleft = ci_uni[1], xright = ci_uni[2], ybottom = 3.6, ytop = 4.4,
     col = adjustcolor("black", alpha.f = 0.2), border = NA)

# dev.off()

