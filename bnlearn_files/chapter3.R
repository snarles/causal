### R code from vignette source '/home/fizban/book-pratiqueR/english/chap3/chapter3.rnw'

###################################################
### code chunk number 1: chapter3.rnw:40-41
###################################################
options(digits = 3)


###################################################
### code chunk number 2: rjags1
###################################################
library(rjags)
sp <- c(0.5, 0.5)
mu <- c(6.1, 6.25)
sigma <- 0.05
jags.data <- list(sp = sp, mu = mu, sigma = sigma,
               cdiam = 6.20)
model1 <- jags.model(file = "inclu.sc.jam", data = jags.data)


###################################################
### code chunk number 3: chapter3.rnw:180-181
###################################################
set.seed(43210)


###################################################
### code chunk number 4: rjags2
###################################################
update(model1, n.iter = 10000)
simu1 <- coda.samples(model = model1, variable.names = "csup",
                      n.iter = 20000, thin = 20)
sim1 <- simu1[[1]]


###################################################
### code chunk number 5: chapter3.rnw:202-203
###################################################
sum(sim1 == 1) / length(sim1)


###################################################
### code chunk number 6: chapter3.rnw:208-209
###################################################
options(digits=4)


###################################################
### code chunk number 7: chapter3.rnw:211-214
###################################################
d.s1 <- dnorm(6.2, mean = mu[1], sd = sigma)
d.s2 <- dnorm(6.2, mean = mu[2], sd = sigma)
d.s1 / (d.s1 + d.s2)


###################################################
### code chunk number 8: rjags3
###################################################
limits <- c(6.16, 6.19)
dsd <- matrix(c(diff(c(0, pnorm(limits, mu[1], sigma), 1)),
                diff(c(0, pnorm(limits, mu[2], sigma), 1))),
              3, 2)
dimnames(dsd) <- list(D = c("thin", "average", "thick"),
                      S = c("s1", "s2"))
dsd


###################################################
### code chunk number 9: rjags3bis
###################################################
jointd <- dsd * sp


###################################################
### code chunk number 10: rjags3ter
###################################################
dds <- t(jointd / rowSums(jointd))
dds


###################################################
### code chunk number 11: pestdag
###################################################
library(bnlearn)
pest.dag <- model2network("[PR][CL][G1|PR:CL][G2|G1][TR|G1][LO|G2:TR]")


###################################################
### code chunk number 12: chapter3.rnw:416-422
###################################################
dat0 <- list(p.PR = c(0.7, 0.2, 0.1),
             a.CL = 3, b.CL = 1,
             g.G1 = c(1, 3, 10),
             k.G2 = 10,
             m.TR = 5, s.TR = 2.5,
             r.LO = 1/3, d.LO = 1)


###################################################
### code chunk number 13: chapter3.rnw:474-475
###################################################
set.seed(123)


###################################################
### code chunk number 14: pest.PR123
###################################################
exp.loss  <- rep(NA, 3)
names(exp.loss) <- paste("PR=", 1:3, sep = "")
qua.loss <- exp.loss
for (PR in 1:3) {
  dat1 <- dat0
  dat1$PR <- PR
  mopest <- jags.model(file = "inclu.pest.jam", data = dat1,
              quiet = TRUE)
  update(mopest, 3000)
  sipest <- 
    coda.samples(model = mopest, variable.names = "LO",
                 n.iter  =  50000)
  summa <- summary(sipest)
  exp.loss[PR] <- summa$statistics["Mean"]
  qua.loss[PR] <- summa$quantiles["75%"]
}#FOR
mean3 <- mean(sipest[[1]][, "LO"])
round(c(exp.loss, MEAN = mean(exp.loss)), 1)


###################################################
### code chunk number 15: pest.PR1
###################################################
set.seed(567)
nbs <- 50000
PR <- 1
g <- function(pr) c(1, 3, 10)[pr] 
CL <- rbeta(nbs, 3, 1)
G1 <- rpois(nbs, CL * g(PR))
G2 <- rpois(nbs, G1 * 10)
il <- function(x) { 
  exp((x - 5) / 2.5)/(1 + exp((x - 5) / 2.5))
}#IL
TR <- rbinom(nbs, 1, il(G1))
x.lo <- G2 * (1 - (1-1/3)*TR)
LO <- rchisq(nbs, 1, ncp = x.lo)
round(mean(LO), 1)


###################################################
### code chunk number 16: pest.PRquantile
###################################################
round(qua.loss)


