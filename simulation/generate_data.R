####
##  Simulate data from a linear (or non-linear) SEM model
####

###
#  Build a DAG using random edges and enforcing monotonicity
###

p <- 5
library(pracma)
adjmat <- (rand(p, p) < 0.5)
adjmat[row(adjmat) >= col(adjmat)] <- 0
colnames(adjmat) <- paste0(LETTERS[1:p], letters[1:p], letters[1:p])
rownames(adjmat) <- colnames(adjmat)
dag <- empty.graph(colnames(adjmat))
amat(dag) <- adjmat
graphviz.plot(dag)
graphviz.plot(cpdag(dag))

###
#  Generate SEM coefficients consistent with the DAG
###

coefmatrix <- randn(p, p) * adjmat
n <- 20
noises <- randn(n, p) - .5

###
#  Function for generating variables from coefmatrix and noise terms
###

## noise is n x p matrix
generate_vars <- function(coefmatrix, noises) {
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
  }
  colnames(ans) <- colnames(coefmatrix)
  ans
}

###
#  Result
###
layout(t(1:2))
n <- 100000
dat <- data.frame(generate_vars(coefmatrix, noises))
res <- gs(dat)
res <- iamb(dat)
res <- hc(dat, restart = 10)
graphviz.plot(res)
graphviz.plot(dag)
