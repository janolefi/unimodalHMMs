## defining color vector
color = c("#E69F00", "#56B4E9", "#009E73", "plum")

set_m_grid <- function(beta) {
  N <- nrow(beta)
  k <- ncol(beta)

  m <- vector("list", N)  # Store mode indices for each row

  for (i in 1:N) {
    first_diff = diff(beta[i,])
    first_diff_sign <- sign(first_diff)
    rle_result <- rle(first_diff_sign)
    change_points <- cumsum(rle_result$lengths)
    mode_indices <- change_points[which(rle_result$values[-length(rle_result$values)] == 1)] + 1

    # exclude stupidly small ones
    threshold = min(beta[i,] + 0.2 * diff(range(beta[i,])))
    mode_indices <- mode_indices[which(beta[i, mode_indices] > threshold)]

    if (beta[i, 1] > beta[i, 2]) mode_indices <- c(1, mode_indices)
    if (beta[i, k] > beta[i, k - 1]) mode_indices <- c(mode_indices, k)

    if (length(mode_indices) == 1) {
      m[[i]] <- max(1, mode_indices - 2):min(k, mode_indices + 2)
    } else if (length(mode_indices) > 1) {
      m[[i]] <- seq(max(1, min(mode_indices) - 2), min(k, max(mode_indices) + 2))
    } else {
      m[[i]] <- integer(0)  # No modes found
    }
  }

  # Generate all possible combinations (Cartesian product)
  grid <- expand.grid(m, KEEP.OUT.ATTRS = FALSE)
  colnames(grid) <- paste0("m", 1:N)  # Name columns dynamically

  return(grid)
}

mlogit <- function(x) {
  s <- 1 - rowSums(x)
  return(log(t(t(x) / s)))
}

tpr_fpr = function(states, stateprobs, tau){
  # states: true simulated state sequence
  # stateprobs: state probabilities
  # tau: threshold for state probabilities
  
  # true positive rate
  tpr = sum(states == 1 & stateprobs[,1] > tau) / sum(states == 1)
  
  # false positive rate
  fpr = sum(states == 2 & stateprobs[,1] > tau) / sum(states == 2)
  
  return(c(tpr, fpr))
}

auc = function(states, stateprobs, npoints = 1000){
  tau_grid = seq(0, 1, length.out = npoints)
  h = diff(tau_grid)[1]
  tpr = rep(NA, npoints)
  fpr = rep(NA, npoints)
  for(i in 1:npoints){
    tpr[i] = tpr_fpr(states, stateprobs, tau_grid[i])[1]
    fpr[i] = tpr_fpr(states, stateprobs, tau_grid[i])[2]
  }
  # AUC = area under the curve
  val <- sum(-diff(fpr) * (tpr[-1] + tpr[-npoints]) / 2)
  val[val < 0.5] <- 1 - val[val < 0.5] # label switch
  return(val)
}

median_CI <- function(x, conf = 0.95) {
  x <- sort(x)
  n <- length(x)
  
  # Compute k1 and k2 from Binomial(n, 0.5)
  alpha <- 1 - conf
  k1 <- qbinom(alpha/2, n, 0.5)
  k2 <- qbinom(1 - alpha/2, n, 0.5)
  
  # Adjust for R 1-indexing
  k1 <- max(1, k1)        # can't be less than 1
  k2 <- min(n, k2 + 1)    # upper bound inclusive
  
  # Return the order-statistic based CI
  c(lower = x[k1], upper = x[k2])
}

