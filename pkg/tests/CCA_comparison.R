## Check results against CCA

data(clonidine)

result <- gslcca(clonidine$spectra, "Critical Exponential",
    time = From, treatment = Treatment, subject.smooth = TRUE,
    data = clonidine$design, subset = Rat == "42")

Y <- result$y[[1]] #power spectra
X <- result$x[[1]] #nonlinear design matrix

## split X into F and G
## specific to partialing out intercept only
partial <- 1:ncol(X) == ncol(X)
F <- X[, !partial]
G <- X[, partial]

## partial out G
Y2 <- residuals(lm(result$y[[1]] ~ 0 + G))
F2 <- residuals(lm(F ~ 0 + G))

## perform regular CCA
cc <- cancor(F2, Y2, xcenter = FALSE, ycenter = FALSE)

## check results

## - same correlation?
all.equal(cc$cor[1], result$cor[1])

## - coefficients of F2 equal coefficients of F? (modulo sign)
B <- result$xcoef[!partial,1]
all.equal(abs(cc$xcoef[,1]), abs(B), check.attributes = FALSE)

## - generalized coefficients of Y2 equal coefficients of Y?
library(MASS)
xscores <- F2 %*% cc$xcoef[,1]
ycoef <- ginv(t(Y2) %*% Y2) %*% t(Y2) %*% xscores
ycoef <- c(ycoef/sqrt(sum((Y2 %*% ycoef)^2)))
all.equal(abs(ycoef), abs(unlist(result$ycoef[,1])), check.attributes = FALSE)

## So, CCA gives obs = Y2 %*% A; fit = X2 %*% B; model fit = p * obs
## Re-express in terms of original Y, i.e. want YA = FB +GC, 
## use linear regression to find C:
A <- result$ycoef[,1]
C <- lm(Y%*%A - F%*%B ~ 0 + G)$coef
all.equal(C, result$xcoef[partial,1], check.attributes = FALSE)

## check standardisation of coefficients
all.equal(sum((Y2 %*% A)^2), 1)
all.equal(sum((F2 %*% B)^2), 1)

