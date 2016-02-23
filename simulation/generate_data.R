####
##  Simulate data from a linear (or non-linear) SEM model
####

library(pracma)
library(bnlearn)
library(InvariantCausalPrediction)

###
#  Build a DAG using random edges and enforcing monotonicity
###

p <- 7
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
mu_noise <- rnorm(p)

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
#  Result of BNLEARN
###
layout(matrix(1:4, 2, 2))
n <- 2000
noises <- t(t(randn(n, p)) + mu_noise)
dat <- data.frame(generate_vars(coefmatrix, noises))
# res <- gs(dat)
# res <- iamb(dat)
res <- hc(dat, restart = 10)
graphviz.plot(dag, main = "truth")
graphviz.plot(res, main = "est")
graphviz.plot(cpdag(dag), main = "truth")
graphviz.plot(cpdag(res), main = "est")



###
#  To add interventions, just change up the noise matrix
###

ind_targ <- floor(p/2) # targeted variable
nints <- 4
neach <- floor(n/nints)
ExpInd <- rep(1:nints, neach)
targeted_vars <- sample((1:p)[-ind_targ], max(ExpInd), replace = FALSE)
n <- length(ExpInd)

noises <- t(t(randn(n, p)) + mu_noise)
for (i in 1:max(ExpInd)) {
  noises[ExpInd==i, targeted_vars[i]] <- noises[ExpInd==i, targeted_vars[i]] + rnorm(1)
}
dat <- generate_vars(coefmatrix, noises)
X <- dat[, -ind_targ]
Y <- dat[, ind_targ]
ICP(X, Y, ExpInd)
coefmatrix[, ind_targ]
