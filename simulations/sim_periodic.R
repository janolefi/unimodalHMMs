## simulating an HSMM and fitting an HMM

library(LaMa)
library(parallel)
library(scales)

source("./utils.R")

beta_tpm <- matrix(c(-3.5, -2, 3, 4, -2, -1, 0, 0, -0.2, -0.2), nrow = 2)
Gamma <- tpm_p(beta = beta_tpm, degree = 2)
d <- ddwell(1:24, Gamma)
plot(d$`state 1`, type = "h")
plot(d1[1:24], type = "h")


x <- 1:24
d1 <- d$`state 1`
d2 <- d$`state 2`

# expected dwell time in state 1
sum((x+1) * d1)


par(mfrow = c(1,2))
plot(x, d1, type = "h", lwd = 2, bty = "n", col = color[1], ylim = c(0, 0.6),
     xlab = "dwell time", ylab = "probabilities", main = "State 1")
plot(x, d2, type = "h", lwd = 2, bty = "n", col = color[2], ylim = c(0, 0.1),
     xlab = "dwell time", ylab = "probabilities", main = "State 2")


# simulate data
mu <- c(1, 15)
sigma <- c(1, 4)

sim_data <- function(nObs) {
  s <- rep(NA, nObs)
  s[1] <- sample(1:2, 1)
  tod <- rep(1:24, ceiling(nObs / 24))[1:nObs]
  for(t in 2:nObs) {
    s[t] <- sample(1:2, 1, prob = Gamma[s[t-1],,tod[t]])
  }
  x <- rgamma2(nObs, mu[s], sigma[s])
  data.frame(x = x, s = s, tod = tod)
}

# penalised likelihood function
pnll = function(par){
  getAll(par,dat)
  Gamma <- tpm(eta)
  delta <- stationary(Gamma)
  
  alpha <- exp(cbind(beta, 0))
  alpha <- alpha / rowSums(alpha); REPORT(alpha)
  
  ind <- which(!is.na(Z[,1]))
  allprobs <- matrix(1, nrow(Z), nrow(alpha))
  allprobs[ind,] <- Z[ind,] %*% t(alpha)
  
  -forward(delta, Gamma, allprobs) + penalty(beta, S, lambda)
}
pnll2 = function(par){
  getAll(par,dat)
  Gamma <- tpm_g(X, beta_tpm); REPORT(beta_tpm)
  delta <- stationary_p(Gamma[,,1:24], t=1)
  
  alpha <- exp(cbind(beta, 0))
  alpha <- alpha / rowSums(alpha); REPORT(alpha)
  
  ind <- which(!is.na(Z[,1]))
  allprobs <- matrix(1, nrow(Z), nrow(alpha))
  allprobs[ind,] <- Z[ind,] %*% t(alpha)
  
  -forward_g(delta, Gamma[,,tod], allprobs) + penalty(beta, S, lambda)
}


k <- 25 # number of basis functions
nSim <- 100
nObs <- 1000
densities1 <- densities2 <- list()

X <- cosinor(1:24, period = c(24, 12))

set.seed(1234)
for(i in 1:nSim) {
  cat("Simulation run", i, "of", nSim, "\n")
  data <- sim_data(nObs)
  
  sDens <- smooth_dens_construct(
    data = data["x"],
    par = list(x = list(mean = mu, sd = sigma)),
    type = "positive",
    k = k,
    knots = list(x = seq(0, max(data$x) * 1.05, length.out = k - 2))
  )
  
  Z <- sDens$Z$x
  S <- sDens$S$x
  Z_p <- sDens$Z_predict$x
  xseq <- sDens$xseq$x

  # fitting unconstrained model only
  par <- list(
    eta = rep(-2, 2),
    beta = sDens$coef$x
  )
  dat <- list(
    Z = Z, S = S, lambda = rep(30, 2), kappa = 0
  )
  mod1 <- qreml(pnll, par, dat, 
               random = "beta", silent = 1)
  par <- list(
    beta_tpm = beta_tpm, 
    beta = sDens$coef$x
  )
  dat$tod <- data$tod
  dat$X <- X
  mod2 <- qreml(pnll2, par, dat, 
                random = "beta", silent = 1)
  
  gc()
  
  # computing the weighted state-dependent densities
  densities1[[i]] <- list(xseq = xseq, 
                          dens = t(mod1$delta * t(Z_p %*% t(mod1$alpha))))
  densities2[[i]] <- list(xseq = xseq, 
                          dens = t(mod2$delta * t(Z_p %*% t(mod2$alpha))))
  cat("\n\n")
}

densities <- densities1


# number of bimodal densities
sum(sapply(densities, function(d) d$dens[which.min(abs(d$xseq - 15)),1]) > 0.002) / nSim # 8 percent second mode in state 1
sum(sapply(densities, function(d) d$dens[1,2]) > 0.02) / nSim # 85 percent second mode in state 2

# get Monte Carlo estimate of state distribution
sim <- sim_data(1e5)
delta <- prop.table(table(sim$s))

# pdf("./simulations/figures/periodic_sim.pdf", width = 7, height = 6)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

# plot dwell-time distributions
plot(x, d1, type = "h", lwd = 2, bty = "n", col = color[1], ylim = c(0, 0.6),
     xlab = "dwell time", ylab = "probabilities", main = "Dwell-time State 1")
plot(x, d2, type = "h", lwd = 2, bty = "n", col = color[2], ylim = c(0, 0.1),
     xlab = "dwell time", ylab = "probabilities", main = "Dwell-time State 2")

par(mar = c(5,4,2.5,2))
# plot estimated and true state-dependent densities
plot(NA, xlim = c(0, 30), ylim = c(0, 0.15), 
     ylab = "Density", xlab = "Observations", bty = "n", main = "Weighted state-dependent densities")
# estimated
for(i in 1:nSim) {
  for(j in 1:2){
    lines(densities[[i]]$xseq, 
          densities[[i]]$dens[,j], col = alpha(color[j], 0.5), lwd = 1, lty = c(1,3)[j])
  }
}
# true
for(j in 1:2) {
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, lwd = 3, lty = c(1,3)[j])
}
legend("topright", bty = "n",
       lwd = c(1,1,3,3), lty = c(1,3,1,3),
       legend = c("State 1 estimated ", "State 2 estimated",
                  "State 1 true", "State 2 true"))

# dev.off()

