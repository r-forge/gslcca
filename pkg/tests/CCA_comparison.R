## Check results against CCA
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

## - generalized coefficients of reduced Y equal coefficients of Y?
library(MASS)
xscores <- X %*% cc$xcoef[,1]
ycoef <- ginv(t(Y) %*% Y) %*% t(Y) %*% xscores
ycoef <- c(ycoef/sqrt(sum((Y %*% ycoef)^2)))
all.equal(abs(ycoef), abs(unlist(result$ycoef[,1])), check.attributes = FALSE)

## check standardisation of coefficients
A <- result$ycoef[,1]
all.equal(sum((Y %*% A)^2), 1)
all.equal(sum((X %*% B)^2), 1)

