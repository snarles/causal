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

# ## hiddenICP seems broken, skipping hereafter
# run_stuff(neach = 1e5, hetero = 0)
# run_stuff(neach = 1e5, hetero = 1)
# 
# 
# ## some interesting cases
# 
# run_stuff("ICP", neach = 1000, hetero = 0)
# run_stuff("ICP", neach = 1000, hetero = 1)
# run_stuff("ICP", neach = 100, hetero = 1)
# run_stuff("ICP", neach = 50, hetero = 1)
# run_stuff("ICP", neach = 20, hetero = 1)
# run_stuff("ICP", neach = 20, hetero = 1, alpha = 0.5)

## automatically extract outcomes
ex_outcomes <- function(res) {
  mr <- NA
  if ("modelReject" %in% names(res)) {
    mr <-   res$modelReject
  }
  cfs <- res$ConfInt
  if (!is.na(mr) && mr) {
    pvs <-   res$pvalues + NA
  } else {
    pvs <- res$pvalues
    pvs[cfs[1, ]==cfs[2, ]] <- 1
  }
  names(pvs) <- colnames(X)
  c(modelReject = mr, pvs)
}

## run experiments

run_exps <- function(mc.reps, mc.cores, procs, neach, hetero, ...) {
  pres <- mclapply(1:mc.reps, function(i) {
    set.seed(i)
    run_stuff(procs, neach, hetero, ...)
  }, mc.cores = mc.cores)
  ans <- list()
  for (st in procs) {
    temp <- lapply(pres, function(v) ex_outcomes(v[[st]]))
    ans[[st]] <- do.call(rbind, temp)
  }
  ans
}

mcc <- 40
mc.reps <- 120
neachs <- c(50, 100, 150, 200)#), 250)#, 300, 350)
heteros <- c(0, 0.1, 0.3, 0.5)#, 1)

res <- array(list(NULL), c(length(neachs), length(heteros)))

for (i in 1:length(neachs)) {
  for (j in 1:length(heteros)) {
    res[[i, j]] <- run_exps(mc.reps, mcc, c("ICP", "hiddenICP"), 
                          neachs[i], 
                          heteros[j], alpha = 0.01)
  }
}

saveRDS(res, file = "aws_results01.rds")

####
##  Interpret results
####

neachs <- c(50, 100, 150, 200)#), 250)#, 300, 350)
heteros <- c(0, 0.1, 0.3, 0.5)#, 1)

res <- readRDS("aws_results01.rds")
res <- res[1:4, 1:4]
res2 <- list()

alph <- 0.05
for (i in 1:dim(res)[1]) {
  for (j in 1:dim(res)[2]) {
    icr <- data.frame(res[[i, j]]$ICP[1:100, ])
    success <- (icr$modelReject == 0) & 
      (icr$X1 < alph) & (icr$X2 < alph) & 
      (icr$Z1 > alph) & (icr$Z2 > alph)
    failure <- (icr$modelReject == 0) & 
      (pmin(icr$Z1, icr$Z2) < alph)
    reject <- icr$modelReject==1
    res2 <- c(res2, 
              list(c(n = neachs[i] * 5,
                     h = heteros[j],
                success = sum(success), 
                failure = sum(failure), 
                reject = sum(reject),
                nullset = sum(!success & !failure & !reject))))
  }
}
res2 <- do.call(rbind, res2)
View(res2)

###
#  Draw the thing
###
cc <- c("green", "red", "black", "white")

par(mar = c(0, 0, 0, 0))
lmat <- matrix(0, 5, 5)
lmat[1, ] <- 1:5
lmat[, 1] <- c(1, 6:9)
lmat[-1, -1] <- 10:25
layout(lmat)
plot(NA, NA, axes = FALSE, ann = FALSE, xlim = c(-1, 1), ylim = c(-1, 1))
for (i in 1:4) {
  plot(NA, NA, axes = FALSE, ann = FALSE, xlim = c(-1, 1), ylim = c(-1, 1))
  text(0, 0, paste0("h=", heteros[i]), cex = 3)
}
for (i in 1:4) {
  plot(NA, NA, axes = FALSE, ann = FALSE, xlim = c(-1, 1), ylim = c(-1, 1))
  text(0, 0, paste0("n=", neachs[i]*5), cex = 3)
}

for (i in 1:16) {
  pie(res2[i, 3:6], labels = rep("", 4), col = cc)
}
