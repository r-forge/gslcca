##### Check results against CCA
library(gslcca)
data(clonidine)

result <- gslcca(spectra, "Critical Exponential",
    time = Time, treatment = Treatment, subject.smooth = TRUE,
    data = clonidine, subset = Rat == "42")

## result return Y and x after 
## a) smoothing Y
## b) partialling out covariates specified by partial
Y <- result$y[[1]] #power spectra
X <- result$x[[1]] #nonlinear design matrix

## perform regular CCA
cc <- cancor(X, Y, xcenter = FALSE, ycenter = FALSE)

## check results

## - same correlation?
all.equal(cc$cor[1], result$cor[1])

## - same coefficients of X? (modulo sign)
B <- result$xcoef[,1]
all.equal(abs(cc$xcoef[,1]), abs(B), check.attributes = FALSE)

## - generalized coefficients of reduced rank Y equal coefficients of Y?
library(MASS)
xscores <- X %*% cc$xcoef[,1]
ycoef <- ginv(t(Y) %*% Y) %*% t(Y) %*% xscores
ycoef <- c(ycoef/sqrt(sum((Y %*% ycoef)^2)))
all.equal(abs(ycoef), abs(unlist(result$ycoef[,1])), check.attributes = FALSE)

## check standardisation of coefficients
A <- result$ycoef[,1]
all.equal(sum((Y %*% A)^2), 1)
all.equal(sum((X %*% B)^2), 1)

## check relation between "correlation" and optimised value
opt <- 1- exp(result$opt[[1]]$value) # opt = log(1 - RSS)
all.equal(result$cor, sqrt(opt))
all.equal(result$cor, lm(Y%*%A ~ 0 + X %*% B)$coef, check.attributes = FALSE)

## - coefficients of Y from CCA same as those obtained from reduced Y?
r <- length(cc$ycoef[,1])
A <- lm(X %*% B ~ I(Y[,1:r]/result$cor) - 1)$coef
all.equal(A/sqrt(sum((Y[,1:r] %*% A)^2)), cc$ycoef[,1], check.attributes = FALSE)

#### Use gslcca for CCA

## dummy function so that optimisation over par does nothing
fn <- function(time, par, X){
    X
}

pop <- as.matrix(LifeCycleSavings[, 2:3])
oec <- as.matrix(LifeCycleSavings[, -(2:3)])
res1 <- cancor(pop, oec)

res2 <- gslcca(oec, formula = ~fn(time, par, pop), time = rep(1, nrow(oec)), 
               partial = ~1, subject.smooth = FALSE, start = list(par = 1))
               
all.equal(abs(res1$xcoef[,1]), abs(c(res2$xcoef)), check.attributes = FALSE)
all.equal(abs(res1$ycoef[,1]), abs(c(res2$ycoef)), check.attributes = FALSE)
all.equal(res1$cor[1], res2$cor)
