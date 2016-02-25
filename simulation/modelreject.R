library(InvariantCausalPrediction)
library(pracma)
neach <- 200
Xs <- randn(neach * 2, 2)
ExpInd <- rep(c(0, 1), each = neach)
Xs[, 1] <- Xs[, 1] + ExpInd
Y <- Xs[, 1] + rnorm(neach * 2)
ICP(Xs, Y, ExpInd)
