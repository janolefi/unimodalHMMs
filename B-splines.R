library(splines2)
library(LaMa) # only for dskewnorm()

cbPalette <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#0072B2", # blue
  "red3",    # red
  "yellow4", # dark yellow/olive
  "plum",    # purple
  "#F0E442", # bright yellow
  "#CC79A7", # reddish purple / pink
  "#D55E00"  # vermillion / strong reddish-orange
)
cbPalette <- rep(cbPalette, 3)


to <- 10
x <- seq(0, to, length = 1000)
d <- 3

k <- 10 # number of basis functions
knots <- seq(0, to, length = k-2)
knots <- knots[2:(length(knots)-1)]
# get B-spline design matrix
B <- bSpline(x, knots = knots, Boundary.knots = c(0,to), degree = d, intercept = TRUE)

# compute integrals of basis functions
w <- colSums(t(t(B))) * diff(x)[1]
# compute normalised B-splines basis
B_norm <- t(t(B) / w)



# plotting

pdf("B-splines_new.pdf", width = 7, height = 5)

par(mfrow = c(3,2),
    mar = c(4.5,4,2.7,1) + 0.1)

lwd = 1.5

# original B-splines basis
plot(x, B[,1], type = "l", bty = "n", ylab = "Basis function", yaxt = "n",
     main = "(a) B-spline basis", ylim = c(0,1), col = cbPalette[1], lwd = lwd)
for (i in 2:ncol(B)) {
  # lines(x, B[,i], col = "lightgray", lwd = 2)
  lines(x, B[,i], col = cbPalette[i], lwd = lwd)
}
axis(2, at = seq(0,1,by=0.2), labels = seq(0,1,by=0.2),las = 1)

# normalised B-splines basis
plot(x, B_norm[,1], type = "l", bty = "n", ylim = c(0, 1), ylab = "Basis function",
     yaxt = "n",
     main = "(b) Normalised basis", lwd = lwd, col = cbPalette[1])
# lines(x, B_norm[,1], col = cbPalette[1], lwd = 2)
for (i in 2:ncol(B_norm)) {
  lines(x, B_norm[,i], col = cbPalette[i], lwd = lwd)
}
axis(2, at = seq(0,1,by=0.2), labels = seq(0,1,by=0.2),las = 1)


# plot weighted B-splines basis -> density
# k <- 15
# knots <- seq(0, to, length = k-2)
# knots <- knots[2:(length(knots)-1)]
# B <- bSpline(x, knots = knots, Boundary.knots = c(0,to), degree = d, intercept = TRUE)

# compute integrals of basis functions
w <- colSums(t(t(B))) * diff(x)[1]
B_norm <- t(t(B) / w)
# expected value of each basis function density
basis_pos <- colSums(x * t(t(B_norm))) * diff(x)[1]

# coefficients
b <- dskewnorm(basis_pos, 2.1, 2.2, 4, log = TRUE)
b <- b - min(b)
a <- exp(b)
a <- a / sum(a)
a[2] <- 0.01
a[9] <- 0.001
a[10] <- 0
a <- a / sum(a)

dspline <- t(t(B_norm) * a)
plot(x, dspline[,1] , type = "l", bty = "n", ylim = c(0, 0.3),
     main = "(c) Unimodal Density", ylab = "Density",
     col = cbPalette[1], lwd = lwd, las = 1)
for (i in 2:ncol(B_norm)) {
  lines(x, dspline[,i], col = cbPalette[i], lwd = lwd)
}
lines(x, rowSums(dspline), lwd = 2.5)


## bottom panel with alphas

plot(basis_pos, a, type = "h", bty = "n", ylim = c(0,0.6), lwd = 2, 
     col = cbPalette[4],
     main = "(d) Coefficients", xlab = "x", ylab = expression(alpha), las = 1)
points(basis_pos[4], a[4], type = "h", lwd = 2.5, col = cbPalette[7])
lines(x, rowSums(dspline)*2, lty = 2)
# axis(4, at = seq(0, 0.6, by = 0.1), labels = seq(0, 0.3, by = 0.05))
# mtext("Density", side = 4, line = 2, cex = 0.7)



# bimodal density
# coefficients
b1 <- dnorm(basis_pos, 2.5, 1, log = TRUE)
b2 <- dnorm(basis_pos, 7.5, 1, log = TRUE)
b1 <- b1 - min(c(b1,b2))
b2 <- b2 - min(c(b1,b2))
a1 <- exp(b1)
a2 <- exp(b2)
a_bm <- a1 + a2
a_bm <- a_bm / sum(a_bm)
a_bm[8] <- 0.1
a_bm[2] <- 0.01
a_bm[9] <- 0.01
a_bm[10] <- 0
a_bm[1] <- 0
a_bm <- a_bm / sum(a_bm)

dspline_bm <- t(t(B_norm) * a_bm)
plot(x, dspline_bm[,1] , type = "l", bty = "n", ylim = c(0, 0.2),
     main = "(e) Bimodal density", ylab = "Density", las = 1,
     col = cbPalette[1], lwd = lwd)
for (i in 2:ncol(B_norm)) {
  lines(x, dspline_bm[,i], col = cbPalette[i], lwd = lwd)
}
lines(x, rowSums(dspline_bm), lwd = 2.5)


## bottom panel with alphas
plot(basis_pos, a_bm, type = "h", bty = "n", ylim = c(0, 0.4), lwd = 2, 
     col = cbPalette[4], las = 1,
     main = "(f) Bimodal coefficients", xlab = "x", ylab = expression(alpha))
# for (i in 2:ncol(B_norm)) {
#   lines(x, dspline_bm[,i], col = cbPalette[i], lwd = lwd)
# }
points(basis_pos[c(4,7)], a_bm[c(4,7)], type = "h", lwd = 2.5, col = cbPalette[7])

lines(x, rowSums(dspline_bm)*2, lty = 2)
# axis(4, at = seq(0, 0.6, by = 0.1), labels = seq(0, 0.3, by = 0.05))
# mtext("Density", side = 4, line = 2, cex = 0.7)

dev.off()




# Initialisation strategy explained ---------------------------------------

par(mfrow = c(1,1))

# Below, we briefly explain how the spline coefficients are initialised using a 
# parametric reference density. This happens automatically when calling
# smooth_dens_construct()

# coefficients
sd = 2 # standard deviation of Gaussian
b <- dnorm(basis_pos, 5, sd, log = TRUE)
b <- b - min(b)
a <- exp(b)

# with all a = 1, we get a straight line for regular basis
# hence for regular basis, if a's are shaped like a Gaussian, this is inherited by the spline function
plot(x, rowSums(B), type = "l", ylim = c(0, 5)) # rowSums(B) = B %*% rep(1, ncol(B))
# but the normalisation changes this, so we need to account for that
lines(x, rowSums(B_norm))
# approximate value of B_norm %*% at basis positions
scaling = sapply(1:k, function(i) rowSums(B_norm)[which.min(abs(basis_pos[i] - x))])
lines(x, B_norm %*% (1 / scaling)) # much better 
# shift does not matter because a will be normalised

# rescale a
a = a / scaling
# sum to 1
a <- a / sum(a)


dspline <- t(t(B_norm) * a)

plot(x, rowSums(dspline), type = "l", bty = "n",
     main = "(c) Sum of basis functions", ylab = "Density", lwd = 2, 
     ylim = c(0, dnorm(5, 5, sd)))
curve(dnorm(x, 5, sd), add = TRUE, col = "blue", n = 500)

