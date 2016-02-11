### R code from vignette source '/home/fizban/book-pratiqueR/english/chap2/chapter2.rnw'

###################################################
### code chunk number 1: dag.bnlearn
###################################################
library(bnlearn)
dag.bnlearn <- model2network("[G][E][V|G:E][N|V][W|V][C|N:W]")
dag.bnlearn


###################################################
### code chunk number 2: cropdag
###################################################
graphviz.plot(dag.bnlearn)


###################################################
### code chunk number 3: inc.dag.bnlearn
###################################################
library(bnlearn)
dag.bnlearn <- model2network("[G][E][V|G:E][N|V][W|V][C|N:W]")
dag.bnlearn


###################################################
### code chunk number 4: tot.indep
###################################################
nano <- nodes(dag.bnlearn)
for (n1 in nano) {
  for (n2 in nano) {
    if (dsep(dag.bnlearn, n1, n2))
      cat(n1, "and", n2, "are independent.\n")
  }#FOR
}#FOR


###################################################
### code chunk number 5: dsep1
###################################################
dsep(dag.bnlearn, "V", "V")


###################################################
### code chunk number 6: con.indep
###################################################
for (n1 in nano[nano != "V"]) { 
  for (n2 in nano[nano != "V"]) {
    if (n1 < n2) {
      if (dsep(dag.bnlearn, n1, n2, "V")) 
        cat(n1, "and", n2, "are independent given V.\n")
    }#THEN
  }#FOR
}#FOR


###################################################
### code chunk number 7: dsep2
###################################################
dsep(dag.bnlearn, "E", "V", "V")


###################################################
### code chunk number 8: path
###################################################
bnlearn::path(dag.bnlearn, from = "E", to = "C")


###################################################
### code chunk number 9: proba1
###################################################
disE <- list(coef = c("(Intercept)" = 50), sd = 10)
disG <- list(coef = c("(Intercept)" = 50), sd = 10)
disV <- list(coef = c("(Intercept)" = -10.35534, 
               E = 0.70711, G = 0.5), sd = 5)
disN <- list(coef = c("(Intercept)" = 45, V = 0.1),
               sd = 9.949874)


###################################################
### code chunk number 10: proba1bis
###################################################
disW <- list(coef = c("(Intercept)" = 15, V = 0.7),
               sd = 7.141428)
disC <- list(coef = c("(Intercept)" = 0, N = 0.3, W = 0.7), 
             sd = 6.25)
dis.list = list(E = disE, G = disG, V = disV, N = disN, 
                 W = disW, C = disC)


###################################################
### code chunk number 11: proba2
###################################################
gbn.bnlearn <- custom.fit(dag.bnlearn, dist = dis.list)


###################################################
### code chunk number 12: proba3
###################################################
gbn.bnlearn$G


###################################################
### code chunk number 13: chapter2.rnw:309-310
###################################################
gbn.bnlearn$C


###################################################
### code chunk number 14: rbmn1
###################################################
library(rbmn)
gbn.rbmn <- bnfit2nbn(gbn.bnlearn)


###################################################
### code chunk number 15: rbmn2
###################################################
gema.rbmn <- nbn2gema(gbn.rbmn)
mn.rbmn <- gema2mn(gema.rbmn)
print8mn(mn.rbmn)


###################################################
### code chunk number 16: rbmn3
###################################################
str(mn.rbmn);


###################################################
### code chunk number 17: data1
###################################################
set.seed(4567)
cropdata1 <- rbn(gbn.bnlearn, n = 200)
set.seed(1234)
cropdata2 <- rbn(gbn.bnlearn, n = 20000)


###################################################
### code chunk number 18: data1d
###################################################
dim(cropdata1)
round(head(cropdata1), 2)


###################################################
### code chunk number 19: estima1
###################################################
est.para <- bn.fit(dag.bnlearn, data = cropdata1)


###################################################
### code chunk number 20: chapter2.rnw:436-437
###################################################
options(digits=3)


###################################################
### code chunk number 21: chapter2.rnw:443-444 (eval = FALSE)
###################################################
## est.para$C <- lm(C ~ N + W, data = cropdata1)


###################################################
### code chunk number 22: chapter2.rnw:505-508 (eval = FALSE)
###################################################
## library(penalized)
## est.para$C <- penalized(C ~ N + W, lambda1 = 0, lambda2 = 1.5, 
##                 data = cropdata1)


###################################################
### code chunk number 23: estima2
###################################################
est.para$E


###################################################
### code chunk number 24: estima3
###################################################
est.para$C


###################################################
### code chunk number 25: estima3bis
###################################################
est.para$C <- lm(C ~ N + W - 1, data = cropdata1)
est.para$C


###################################################
### code chunk number 26: estima4
###################################################
lmC <- lm(C ~ N + W, data = cropdata1[, c("N", "W", "C")])
coef(lmC)


###################################################
### code chunk number 27: estima5
###################################################
confint(lmC)


###################################################
### code chunk number 28: chapter2.rnw:612-613
###################################################
cormat <- cor(cropdata1[, c("C", "W", "N")])


###################################################
### code chunk number 29: chapter2.rnw:621-624
###################################################
library(corpcor)
invcor <- cor2pcor(cormat)
dimnames(invcor) <- dimnames(cormat)


###################################################
### code chunk number 30: chapter2.rnw:626-627
###################################################
invcor


###################################################
### code chunk number 31: chapter2.rnw:639-640
###################################################
ci.test("C", "W", "N", test = "cor", data = cropdata1)


###################################################
### code chunk number 32: learning1
###################################################
stru1 <- iamb(cropdata1, test = "cor")


###################################################
### code chunk number 33: learning4
###################################################
wl <- matrix(c("V", "N"), ncol = 2)
wl
stru2 <- iamb(cropdata1, test = "cor", whitelist = wl)
all.equal(dag.bnlearn, stru2)


###################################################
### code chunk number 34: learning6
###################################################
dim(cropdata2)
stru3 <- iamb(cropdata2, test = "cor")
all.equal(dag.bnlearn, stru3)


###################################################
### code chunk number 35: chapter2.rnw:718-719
###################################################
options(digits=8)


###################################################
### code chunk number 36: chapter2.rnw:755-757
###################################################
score(dag.bnlearn, data = cropdata2, type = "bic-g")
score(dag.bnlearn, data = cropdata2, type = "bge")


###################################################
### code chunk number 37: chapter2.rnw:771-772
###################################################
options(digits = 3)


###################################################
### code chunk number 38: rbmn1
###################################################
print8nbn(gbn.rbmn)


###################################################
### code chunk number 39: rbmn2
###################################################
print8gema(gema.rbmn)


###################################################
### code chunk number 40: rbmn3
###################################################
print8mn(condi4joint(mn.rbmn, par = "C", pour = "V", x2 = 80))
print8mn(condi4joint(mn.rbmn, par = "V", pour = "C", x2 = 80))


###################################################
### code chunk number 41: rbmn4
###################################################
unlist(condi4joint(mn.rbmn, par = "C", pour = "V", x2 = NULL))


###################################################
### code chunk number 42: chapter2.rnw:892-893
###################################################
set.seed(1234)


###################################################
### code chunk number 43: simu1
###################################################
nbs <- 4
VG <- rnorm(nbs, mean = 50, sd = 10)
VE <- rnorm(nbs, mean = 50, sd = 10)
VV <- rnorm(nbs, mean = -10.355 + 0.5 * VG + 0.707 * VE, 
        sd = 5)
VN <- rnorm(nbs, mean = 45 + 0.1 * VV, sd = 9.95)
cbind(VV, VN)


###################################################
### code chunk number 44: chapter2.rnw:906-907
###################################################
set.seed(1234)


###################################################
### code chunk number 45: simu2
###################################################
sim <- rbn(gbn.bnlearn, n = 4)
sim[, c("V", "N")]


###################################################
### code chunk number 46: data2
###################################################
set.seed(4567)
cropdata1 <- rbn(gbn.bnlearn, n = 200)
set.seed(1234)
cropdata2 <- rbn(gbn.bnlearn, n = 20000)


###################################################
### code chunk number 47: simu3
###################################################
head(cpdist(gbn.bnlearn, nodes = c("C", "N", "W"),
       evidence = (C > 80)))


###################################################
### code chunk number 48: simu4
###################################################
head(cpdist(gbn.bnlearn, nodes = c("V"), 
       evidence = list(G = 10, E = 90), method = "lw"))


###################################################
### code chunk number 49: simu5
###################################################
cpquery(gbn.bnlearn, event = (V > 70),
  evidence = list(G = 10, E = 90), method = "lw")


###################################################
### code chunk number 50: igraph1
###################################################
library(igraph)
igraph.options(print.full = TRUE)
dag0.igraph <- graph.formula(G-+V, E-+V, V-+N, V-+W, 
                               N-+C, W-+C)
dag0.igraph


###################################################
### code chunk number 51: igraph2
###################################################
dag.igraph <- igraph.from.graphNEL(as.graphNEL(dag.bnlearn))


###################################################
### code chunk number 52: igraph3
###################################################
V(dag.igraph)
E(dag.igraph)


###################################################
### code chunk number 53: igraph4
###################################################
par(mfrow = c(2, 2), mar = rep(3, 4), cex.main = 2)
plot(dag.igraph, main = "\n1: defaults")
dag2 <- dag.igraph
V(dag2)$label <- V(dag2)$name
plot(dag2, main = "\n2: with labels")
ly <- matrix(c(2, 3, 1, 1, 2, 3,
               1, 4, 4, 2, 3, 2), 6)
plot(dag2, layout = ly, main = "\n3: positioning")
colo <- c("black", "darkgrey", "darkgrey", rep(NA, 3))
lcolo <- c(rep("white", 3), rep(NA, 3))
par(mar = rep(0, 4), lwd = 1.5)
plot(dag2, layout = ly, frame = TRUE,
     main = "\n4: final",
     vertex.color = colo, vertex.label.color = lcolo,
     vertex.label.cex = 3, vertex.size = 50,
     edge.arrow.size = 0.8, edge.color = "black")


###################################################
### code chunk number 54: gausplot1
###################################################
gbn.fit <- bn.fit(dag.bnlearn, cropdata2)
bn.fit.qqplot(gbn.fit)


###################################################
### code chunk number 55: gausplot2
###################################################
bn.fit.qqplot(gbn.fit$V)


###################################################
### code chunk number 56: chapter2.rnw:1132-1137
###################################################
condi4joint <- function(...) {
  res <- rbmn::condi4joint(...)
  res$rho = zapsmall(res$rho)
  return(res)
}


###################################################
### code chunk number 57: rbmn5
###################################################
C.EV <- condi4joint(mn.rbmn, par = "C", pour = c("E", "V"),
          x2 = NULL)
C.EV$rho


###################################################
### code chunk number 58: dsep1
###################################################
dsep(gbn.bnlearn, "E", "C", "V")


###################################################
### code chunk number 59: graphe1
###################################################
set.seed(5678)
cropdata3 <- cpdist(gbn.bnlearn, nodes = c("E", "V", "C"),
                    evidence = TRUE, n = 1000)
plot(cropdata3$V, cropdata3$C, type = "n",
     main = "C | V, E; E is the point size")
cexlim <- c(0.1, 2.4)
cexE <- cexlim[1] + diff(cexlim) / diff(range(cropdata3$E)) *
                    (cropdata3$E - min(cropdata3$E))
points(cropdata3$V, cropdata3$C, cex = cexE)
cqa <- quantile(cropdata3$C, seq(0, 1, 0.1))
abline(h = cqa, lty = 3)


