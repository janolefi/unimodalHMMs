library(splines2)
library(LaMa) # only for dskewnorm()

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "red3", "yellow4" ,"plum", "#F0E442")
cbPalette <- rep(cbPalette, 3)


to <- 10
x <- seq(0, to, length = 1000)
d <- 3

k <- 8 # number of basis functions
knots <- seq(0, to, length = k-2)
knots <- knots[2:(length(knots)-1)]
# get B-spline design matrix
B <- bSpline(x, knots = knots, Boundary.knots = c(0,to), degree = d, intercept = TRUE)

# compute integrals of basis functions
w <- colSums(t(t(B))) * diff(x)[1]
# compute normalised B-splines basis
B_norm <- t(t(B) / w)



# plotting

# pdf("B-splines.pdf", width = 7, height = 3)

par(mfrow = c(1,3),
    mar = c(5,4,4,0.5))

lwd = 1.5

# original B-splines basis
plot(x, B[,1], type = "l", bty = "n", ylab = "Basis function",
     main = "(a) B-spline basis", ylim = c(0,2), col = cbPalette[1], lwd = lwd)
for (i in 2:ncol(B)) {
  # lines(x, B[,i], col = "lightgray", lwd = 2)
  lines(x, B[,i], col = cbPalette[i], lwd = lwd)
}

# normalised B-splines basis
plot(x, B_norm[,1], type = "l", bty = "n", ylim = c(0, 2), ylab = "Basis function",
     main = "(b) Normalised basis", lwd = lwd, col = cbPalette[1])
# lines(x, B_norm[,1], col = cbPalette[1], lwd = 2)
for (i in 2:ncol(B_norm)) {
  lines(x, B_norm[,i], col = cbPalette[i], lwd = lwd)
}


# plot weighted B-splines basis -> density
k <- 15
knots <- seq(0, to, length = k-2)
knots <- knots[2:(length(knots)-1)]
B <- bSpline(x, knots = knots, Boundary.knots = c(0,to), degree = d, intercept = TRUE)

# compute integrals of basis functions
w <- colSums(t(t(B))) * diff(x)[1]
B_norm <- t(t(B) / w)
# expected value of each basis function density
basis_pos <- colSums(x * t(t(B_norm))) * diff(x)[1]

# coefficients
b <- dskewnorm(basis_pos, 2.2, 2.1, 4, log = TRUE)
b <- b - min(b)
a <- exp(b)
a <- a / sum(a)

dspline <- t(t(B_norm) * a)
plot(x, dspline[,1] , type = "l", bty = "n", ylim = c(0, 0.35),
     main = "(c) Sum of basis functions", ylab = "Density",
     col = cbPalette[1], lwd = 1)
for (i in 2:ncol(B_norm)) {
  lines(x, dspline[,i], col = cbPalette[i], lwd = 1)
}
lines(x, rowSums(dspline), lwd = lwd)

# dev.off()




# Initialisation strategy explained ---------------------------------------

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

