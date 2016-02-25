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
library(MASS)

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
    Ynew <- ans[, 1:(i-1), drop = FALSE] %*% coefmatrix[1:(i-1), i, drop = FALSE]
    ans[, i] <- Ynew + (1 + hetero * Ynew) * noises[, i]
  }
  colnames(ans) <- colnames(coefmatrix)
  ans
}

generate_manipulations <- function(coefmatrix, mus, vars, del.mus, del.vars, ExpInd, hetero = 0) {
  n <- length(ExpInd)
  ki <- max(ExpInd)
  pieces <- lapply(1:ki, function(ii) {
    mvrnorm(n = sum(ExpInd==ii), 
            mu = mus + del.mus[ii, ], 
            Sigma = diag(vars + del.vars[ii, ]))
  })
  noises <- do.call(rbind, pieces)
  ans <- generate_vars(coefmatrix, noises, hetero = hetero)
  ans
}

## Setup graph
nms <- c("X1", "X2", "Y", "Z1", "Z2")
Yind <- 3
p <- length(nms)
adjmat <- rbind(c(0,0,1,1,0),
                c(0,0,1,0,1),
                c(0,0,0,1,1),
                c(0,0,0,0,0),
                c(0,0,0,0,0))
colnames(adjmat) <- nms
rownames(adjmat) <- colnames(adjmat)
dag <- empty.graph(colnames(adjmat))
amat(dag) <- adjmat

## Setup SEM
coefmatrix <- adjmat ## coefficients all 1
mus <- rep(0, p)
vars <- rep(1, p)

## Setup manipulation
del.mus <- rbind(0, eye(5)[c(1, 2, 4, 5), ])
del.vars <- 0 * del.mus

# ## Demo run
# ExpInd <- rep(1:5, each = 1000)
# dat <- generate_manipulations(coefmatrix, mus, vars, del.mus, del.vars, ExpInd, hetero = 0)
# dat <- generate_manipulations(coefmatrix, mus, vars, del.mus, del.vars, ExpInd, hetero = 1)
# Y <- dat[, Yind]
# X <- dat[, -Yind]
# (resICP <- ICP(X, Y, ExpInd))
# (resH <- hiddenICP(X, Y, ExpInd))

####
##  Simulation code
####

run_stuff <- function(procs = c("ICP", "hiddenICP"), neach, hetero = 0, ...) {
  ExpInd <- rep(1:5, each = neach)
  dat <- generate_manipulations(coefmatrix, mus, vars, del.mus, del.vars, ExpInd, 
                                hetero = hetero)
  Y <- dat[, Yind]
  X <- dat[, -Yind]
  ans <- list()
  if ("ICP" %in% procs) {
    ans[["ICP"]] <- ICP(X, Y, ExpInd, ...)
  }
  if ("hiddenICP" %in% procs) {
    ans[["hiddenICP"]] <- hiddenICP(X, Y, ExpInd, ...)
  }
  ans
}

run_stuff("ICP", neach = 200, hetero = 0)$ICP
