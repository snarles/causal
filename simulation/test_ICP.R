####
##  Robustness test for ICP
####
#  Model:
#   Y = Xbeta + N(0, sigma^2(1 + Xbeta)^2)
##

library(pracma)
library(bnlearn)
library(InvariantCausalPrediction)
library(parallel)

## generate variables with hetero parameter
generate_vars <- function(coefmatrix, noises, hetero = 0) {
  p <- dim(coefmatrix)[1]
  stopifnot(p == dim(noises)[2])
  if (is.null(dim(noises))) {
    noises <- t(noises)
  } else {
    n <- dim(noises)[1]
  }
  ans <- matrix(0, n, p)
  ans[, 1] <- noises[, 1]
  for (i in 2:p) {
    ans[, i] <- ans[, 1:(i-1), drop = FALSE] %*% coefmatrix[1:(i-1), i, drop = FALSE] + 
      noises[, i]
    ans[, i] <- Ynew + (1 + hetero * Ynew) * noises[, i]
  }
  colnames(ans) <- colnames(coefmatrix)
  ans
}

## Setup graph
nms <- c("X1", "X2", "Y", "Z1", "Z2")
p <- length(nms)
adjmat <- rbind(c(0,0,1,1,0),
                c(0,0,1,0,1),
                c(0,0,0,1,0),
                c(0,0,0,0,1))
colnames(adjmat) <- nms
rownames(adjmat) <- colnames(adjmat)
dag <- empty.graph(colnames(adjmat))
amat(dag) <- adjmat
graphviz.plot(dag)
graphviz.plot(cpdag(dag))


## Setup SEM
coefmatrix <- adjmat ## coefficients all 1
mus <- rep(0, p)
vars <- rep(1, p)

## Setup manipulation
del.mus <- rbind(0, eye(5)[c(1, 2, 4, 5)])
