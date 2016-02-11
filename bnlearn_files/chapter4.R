### R code from vignette source '/home/fizban/book-pratiqueR/english/chap4/chapter4.rnw'

###################################################
### code chunk number 1: chapter4.rnw:4-6
###################################################
library(bnlearn)
options(width = 65)


###################################################
### code chunk number 2: chapter4.rnw:233-236
###################################################
X <- paste("[X1][X3][X5][X6|X8][X2|X1][X7|X5][X4|X1:X2]",
      "[X8|X3:X7][X9|X2:X7][X10|X1:X9]", sep = "")
dag <- model2network(X)


###################################################
### code chunk number 3: chapter4.rnw:239-241
###################################################
skel <- skeleton(dag)
vstructs(dag)


###################################################
### code chunk number 4: chapter4.rnw:247-248
###################################################
cp1 <- cpdag(dag)


###################################################
### code chunk number 5: chapter4.rnw:283-290
###################################################
dag2 <- dag
dag2 <- set.arc(dag2, "X7", "X5")
dag2 <- set.arc(dag2, "X4", "X2")
dag2 <- set.arc(dag2, "X1", "X2")
dag2 <- set.arc(dag2, "X1", "X4")
cp2 <- cpdag(dag2)
all.equal(cp1, cp2)


###################################################
### code chunk number 6: chapter4.rnw:358-359
###################################################
dsep(dag, x = "X9", y = "X5", z = c("X2", "X7", "X10"))


###################################################
### code chunk number 7: chapter4.rnw:419-421
###################################################
mb(dag, node = "X9")
mb(dag, node = "X7")


###################################################
### code chunk number 8: chapter4.rnw:427-429
###################################################
par.X9 <- parents(dag, node = "X9")
ch.X9 <- children(dag, node = "X9")


###################################################
### code chunk number 9: chapter4.rnw:433-434
###################################################
sp.X9 <- sapply(ch.X9, parents, x = dag)


###################################################
### code chunk number 10: chapter4.rnw:439-441
###################################################
sp.X9 <- sp.X9[sp.X9 != "X9"]
unique(c(par.X9, ch.X9, sp.X9))


###################################################
### code chunk number 11: chapter4.rnw:449-452
###################################################
V <- setdiff(nodes(dag), "X9")
S <- mb(dag, "X9")
sapply(setdiff(V, S), dsep, bn = dag, y = "X9", z = S)


###################################################
### code chunk number 12: chapter4.rnw:455-458
###################################################
V <- setdiff(nodes(dag), "X7")
S <- mb(dag, "X7")
sapply(setdiff(V, S), dsep, bn = dag, y = "X7", z = S)


###################################################
### code chunk number 13: chapter4.rnw:472-476
###################################################
belongs <- logical(0)
for (node in S)
  belongs[node] <- "X7" %in% mb(dag, node)
belongs


###################################################
### code chunk number 14: chapter4.rnw:517-518
###################################################
mg1 <- moral(dag)


###################################################
### code chunk number 15: chapter4.rnw:523-525
###################################################
all.equal(moral(dag),
          moral(set.arc(dag, from = "X7", to = "X3")))


###################################################
### code chunk number 16: chapter4.rnw:531-536
###################################################
mg2 <- dag
vs <- vstructs(dag)
for (i in seq(nrow(vs)))
  mg2 <- set.edge(mg2, from = vs[i, "X"], to = vs[i, "Y"],
           check.cycles = FALSE)


###################################################
### code chunk number 17: chapter4.rnw:541-543
###################################################
mg2 <- skeleton(mg2)
all.equal(mg1, mg2)


###################################################
### code chunk number 18: chapter4.rnw:1041-1056
###################################################
set.seed(4567);
dag.bnlearn <- model2network("[G][E][V|G:E][N|V][W|V][C|N:W]")
disE <- list(coef = c("(Intercept)" = 50), sd = 10)
disG <- list(coef = c("(Intercept)" = 50), sd = 10)
disV <- list(coef = c("(Intercept)" = -10.35534,
             E = 0.70711, G = 0.5), sd = 5)
disN <- list(coef = c("(Intercept)" = 45, V = 0.1), sd = 9.949874)
disW <- list(coef = c("(Intercept)" = 15, V = 0.7), sd = 7.141428)
disC <- list(coef = c("(Intercept)" = 0, N = 0.3, W = 0.7),
             sd = 6.25);
dis.liste = list(E = disE, G = disG, V = disV, N = disN,
                 W = disW, C = disC)
gbn.bnlearn <- custom.fit(dag.bnlearn, dist = dis.liste)
cropdata1 <- cpdist(gbn.bnlearn, nodes = nodes(gbn.bnlearn),
                    evidence = TRUE, n = 200)


###################################################
### code chunk number 19: chapter4.rnw:1067-1069
###################################################
bn.cor <- gs(cropdata1, test = "cor", alpha = 0.05)
modelstring(bn.cor)


###################################################
### code chunk number 20: chapter4.rnw:1077-1081
###################################################
bn.zf <- gs(cropdata1, test = "zf", alpha = 0.05)
bn.mc <- gs(cropdata1, test = "mc-cor", B = 1000)
all.equal(bn.cor, bn.zf)
all.equal(bn.cor, bn.mc)


###################################################
### code chunk number 21: chapter4.rnw:1086-1088
###################################################
bn.iamb <- iamb(cropdata1, test = "cor", alpha = 0.05)
all.equal(bn.cor, bn.iamb)


###################################################
### code chunk number 22: chapter4.rnw:1125-1127
###################################################
ci.test("N", "V", test = "cor", data = cropdata1)
ci.test("N", "V", "C", test = "cor", data = cropdata1)


###################################################
### code chunk number 23: chapter4.rnw:1132-1135
###################################################
bn.cor <- gs(cropdata1, test = "cor", alpha = 0.05, 
            whitelist = c("V", "N"))
all.equal(bn.cor, dag.bnlearn)


###################################################
### code chunk number 24: chapter4.rnw:1208-1209
###################################################
survey <- read.table("../chap1/survey.txt", header = TRUE)


###################################################
### code chunk number 25: chapter4.rnw:1220-1223
###################################################
learned <- hc(survey, score = "bic")
modelstring(learned)
score(learned, data = survey, type = "bic")


###################################################
### code chunk number 26: chapter4.rnw:1295-1299
###################################################
survey.dag <- model2network("[A][S][E|A:S][O|E][R|E][T|O:R]")
learned.start <- hc(survey, score = "bic", start = survey.dag)
modelstring(learned.start)
all.equal(cpdag(learned), cpdag(learned.start))


###################################################
### code chunk number 27: chapter4.rnw:1303-1304
###################################################
set.seed(1234)


###################################################
### code chunk number 28: chapter4.rnw:1306-1307
###################################################
hc(survey, score = "bic", start = random.graph(names(survey)))


###################################################
### code chunk number 29: chapter4.rnw:1375-1376
###################################################
mmhc(survey)


###################################################
### code chunk number 30: chapter4.rnw:1381-1382
###################################################
rsmax2(survey, restrict = "mmpc", maximize = "hc")


###################################################
### code chunk number 31: chapter4.rnw:1391-1394
###################################################
rsmax2(survey, restrict = "si.hiton.pc", test = "x2",
  maximize = "tabu", score = "bde", 
  maximize.args = list(iss = 5))


###################################################
### code chunk number 32: chapter4.rnw:1702-1737
###################################################
res <- model2network("[A][S][E|A:S][O|E][R|E][T|O:R]")

A.lv <- c("young", "adult", "old")
S.lv <- c("M", "F")
E.lv <- c("high", "uni")
O.lv <- c("emp", "self")
R.lv <- c("small", "big")
T.lv <- c("car", "train", "other")

A.prob <- matrix(c(0.30, 0.50, 0.20), ncol = 3,
            dimnames = list(c(""), A = A.lv))

S.prob <- matrix(c(0.60, 0.40), ncol = 2,
            dimnames = list(c(""), S = S.lv))

E.prob <- array(c(0.75, 0.25, 0.72, 0.28, 0.88, 0.12, 0.64, 0.36,
            0.70, 0.30, 0.90, 0.10), dim = c(2, 3, 2),
            dimnames = list(E = E.lv, A = A.lv, S = S.lv))

O.prob <- matrix(c(0.96, 0.04, 0.92, 0.08), nrow = 2,
            dimnames = list(O = O.lv, E = E.lv))

R.prob <- matrix(c(0.25, 0.75, 0.20, 0.80), nrow = 2,
            dimnames = list(R = R.lv, E = E.lv))

T.prob <- array(c(0.48, 0.42, 0.10, 0.56, 0.36, 0.08, 0.58, 0.24,
            0.18, 0.70, 0.21, 0.09), dim = c(3, 2, 2),
            dimnames = list(T = T.lv, O = O.lv, R = R.lv))

res3 <- model2network("[A][S][E|A:S][O|E][R|E][T|O:R]")

tdp <- list(A = A.prob, S = S.prob, E = E.prob, O = O.prob, R = R.prob, T = T.prob)
bn <- custom.fit(res3, tdp)
set.seed(123)
options(digits = 5)


###################################################
### code chunk number 33: chapter4.rnw:1746-1748
###################################################
cpquery(bn, event = (S == "M") & (T == "car"),
          evidence = (E == "high"), n = 10^6)


###################################################
### code chunk number 34: chapter4.rnw:1754-1756
###################################################
particles <- rbn(bn, 10^6)
head(particles, n = 5)


###################################################
### code chunk number 35: chapter4.rnw:1761-1763
###################################################
partE <- particles[(particles[, "E"] == "high"), ]
nE <- nrow(partE)


###################################################
### code chunk number 36: chapter4.rnw:1768-1771
###################################################
partEq <- 
  partE[(partE[, "S"] == "M") & (partE[, "T"] == "car"), ]
nEq <- nrow(partEq)


###################################################
### code chunk number 37: chapter4.rnw:1775-1776
###################################################
nEq/nE


###################################################
### code chunk number 38: chapter4.rnw:1820-1822
###################################################
mutbn <- mutilated(bn, list(E = "high"))
mutbn$E


###################################################
### code chunk number 39: chapter4.rnw:1828-1829
###################################################
set.seed(123)


###################################################
### code chunk number 40: chapter4.rnw:1831-1836
###################################################
particles <- rbn(bn, 10^6)
partQ <- particles[(particles[, "S"] == "M") & 
                  (particles[, "T"] == "car"), ]
nQ <- nrow(partQ)
nQ/10^6


###################################################
### code chunk number 41: chapter4.rnw:1843-1844
###################################################
w <- logLik(bn, particles, nodes = "E", by.sample = TRUE)


###################################################
### code chunk number 42: chapter4.rnw:1848-1852
###################################################
wEq <- sum(exp(w[(particles[, "S"] == "M") & 
                (particles[, "T"] == "car")]))
wE <- sum(exp(w))
wEq/wE


###################################################
### code chunk number 43: chapter4.rnw:1862-1863
###################################################
set.seed(678)


###################################################
### code chunk number 44: chapter4.rnw:1865-1867
###################################################
cpquery(bn, event = (S == "M") & (T == "car"),
          evidence = list(E = "high"), method = "lw")


###################################################
### code chunk number 45: chapter4.rnw:1957-1959
###################################################
data(marks)
head(marks)


###################################################
### code chunk number 46: chapter4.rnw:1979-1984
###################################################
latent <- factor(c(rep("A", 44), "B", 
                  rep("A", 7), rep("B", 36)))
modelstring(hc(marks[latent == "A", ]))
modelstring(hc(marks[latent == "B", ]))
modelstring(hc(marks))


###################################################
### code chunk number 47: chapter4.rnw:1996-1998
###################################################
dmarks <- discretize(marks, breaks = 2, method = "interval")
modelstring(hc(cbind(dmarks, LAT = latent)))


