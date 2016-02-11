### R code from vignette source '/home/fizban/book-pratiqueR/english/chap6/chapter6.rnw'

###################################################
### code chunk number 1: chapter6.rnw:90-93
###################################################
library(bnlearn)
sachs <- read.table("sachs.data.txt", header = TRUE)
head(sachs)


###################################################
### code chunk number 2: chapter6.rnw:104-107
###################################################
dag.iamb <- inter.iamb(sachs, test = "cor")
narcs(dag.iamb)
directed.arcs(dag.iamb)


###################################################
### code chunk number 3: chapter6.rnw:115-121
###################################################
sachs.modelstring <- 
  paste("[PKC][PKA|PKC][Raf|PKC:PKA][Mek|PKC:PKA:Raf]",
        "[Erk|Mek:PKA][Akt|Erk:PKA][P38|PKC:PKA]",
        "[Jnk|PKC:PKA][Plcg][PIP3|Plcg][PIP2|Plcg:PIP3]")
dag.sachs <- model2network(sachs.modelstring)
unlist(compare(dag.sachs, dag.iamb))


###################################################
### code chunk number 4: chapter6.rnw:126-127
###################################################
unlist(compare(skeleton(dag.sachs), skeleton(dag.iamb)))


###################################################
### code chunk number 5: chapter6.rnw:246-248 (eval = FALSE)
###################################################
## dsachs <- discretize(sachs, method = "hartemink",
##             breaks = 3, ibreaks = 60, idisc = "quantile")


###################################################
### code chunk number 6: chapter6.rnw:250-254
###################################################
dsachs <- read.table("sachs.discretised.txt", header = TRUE)
options(width = 60)
set.seed(123)
options(digits = 4)


###################################################
### code chunk number 7: chapter6.rnw:274-276
###################################################
boot <- boot.strength(dsachs, R = 500, algorithm = "hc",
          algorithm.args = list(score = "bde", iss = 10))


###################################################
### code chunk number 8: chapter6.rnw:284-285
###################################################
options(digits = 6)


###################################################
### code chunk number 9: chapter6.rnw:294-295
###################################################
boot[(boot$strength > 0.85) & (boot$direction >= 0.5), ]


###################################################
### code chunk number 10: chapter6.rnw:308-309
###################################################
avg.boot <- averaged.network(boot, threshold = 0.85)


###################################################
### code chunk number 11: chapter6.rnw:318-319
###################################################
avg.boot <- skeleton(avg.boot)


###################################################
### code chunk number 12: chapter6.rnw:341-347 (eval = FALSE)
###################################################
## nodes <- names(dsachs)
## start <- random.graph(nodes = nodes, method = "ic-dag", 
##            num = 500, every = 50)
## netlist <- lapply(start, function(net) {
##   hc(dsachs, score = "bde", iss = 10, start = net)
## })


###################################################
### code chunk number 13: chapter6.rnw:354-357 (eval = FALSE)
###################################################
## rnd <- custom.strength(netlist, nodes = nodes)
## rnd[(rnd$strength > 0.85) & (rnd$direction >= 0.5), ]
## avg.start <- averaged.network(rnd, threshold = 0.85)


###################################################
### code chunk number 14: chapter6.rnw:359-362
###################################################
load("sachs.avg.start.rda")
rnd[(rnd$strength > 0.85) & (rnd$direction >= 0.5), ]
avg.start <- averaged.network(rnd, threshold = 0.85)


###################################################
### code chunk number 15: chapter6.rnw:370-371
###################################################
all.equal(cpdag(avg.boot), cpdag(avg.start))


###################################################
### code chunk number 16: chapter6.rnw:378-380
###################################################
score(cextend(cpdag(avg.start)), dsachs, type = "bde",
  iss = 10)


###################################################
### code chunk number 17: chapter6.rnw:398-409 (eval = FALSE)
###################################################
## library(catnet)
## netlist <- vector(500, mode = "list")
## ndata <- nrow(dsachs)
## nodes <- names(dsachs)
## netlist <- lapply(netlist, function(net) {
##   boot <- dsachs[sample(ndata, replace = TRUE), ]
##   top.ord <- cnSearchOrder(boot)
##   best <- cnFindBIC(top.ord, ndata)
##   cnMatEdges(best)
## })
## sann <- custom.strength(netlist, nodes = nodes)


###################################################
### code chunk number 18: chapter6.rnw:412-413
###################################################
load("sachs.avg.catnet.rda")


###################################################
### code chunk number 19: chapter6.rnw:436-438
###################################################
sann[(sann$strength > 0.85) & (sann$direction >= 0.5), ]
avg.catnet <- averaged.network(sann, threshold = 0.85)


###################################################
### code chunk number 20: chapter6.rnw:449-451
###################################################
narcs(avg.catnet)
narcs(avg.start)


###################################################
### code chunk number 21: chapter6.rnw:457-461
###################################################
score(cextend(cpdag(avg.catnet)), dsachs, type = "bde",
  iss = 10)
score(cextend(cpdag(avg.start)), dsachs, type = "bde",
  iss = 10)


###################################################
### code chunk number 22: chapter6.rnw:475-477
###################################################
all.equal(averaged.network(boot, threshold = 0.50),
          averaged.network(boot, threshold = 0.70))


###################################################
### code chunk number 23: chapter6.rnw:487-488
###################################################
averaged.network(boot)


###################################################
### code chunk number 24: chapter6.rnw:555-558
###################################################
unlist(compare(cpdag(dag.sachs), cpdag(avg.boot)))
unlist(compare(cpdag(dag.sachs), 
               cpdag(averaged.network(boot))))


###################################################
### code chunk number 25: chapter6.rnw:601-603
###################################################
isachs <- read.table("sachs.interventional.txt", 
            header = TRUE, colClasses = "factor")


###################################################
### code chunk number 26: chapter6.rnw:610-613
###################################################
wh <- matrix(c(rep("INT", 11), names(isachs)[1:11]), ncol = 2)
dag.wh <- tabu(isachs, whitelist = wh, score = "bde",
           iss = 10, tabu = 50)


###################################################
### code chunk number 27: chapter6.rnw:628-632
###################################################
tiers <- list("INT", names(isachs)[1:11])
bl <- tiers2blacklist(nodes = tiers)
dag.tiers <- tabu(isachs, blacklist = bl,
              score = "bde", iss = 1, tabu = 50)


###################################################
### code chunk number 28: chapter6.rnw:663-667
###################################################
INT <- sapply(1:11, function(x) {
                          which(isachs$INT == x) })
nodes <- names(isachs)[1:11]
names(INT) <- nodes


###################################################
### code chunk number 29: chapter6.rnw:677-685 (eval = FALSE)
###################################################
## start <- random.graph(nodes = nodes, method = "melancon",
##            num = 500, burn.in = 10^5, every = 100)
## netlist <- lapply(start, function(net) {
##   tabu(isachs[, 1:11], score = "mbde", exp = INT,
##      iss = 1, start = net, tabu = 50)
## })
## intscore <- custom.strength(netlist, nodes = nodes,
##               cpdag = FALSE)


###################################################
### code chunk number 30: chapter6.rnw:687-688
###################################################
load("sachs.intscore.rda")


###################################################
### code chunk number 31: chapter6.rnw:717-719
###################################################
dag.mbde <- averaged.network(intscore)
unlist(compare(dag.sachs, dag.mbde))


###################################################
### code chunk number 32: chapter6.rnw:754-758
###################################################
isachs <- isachs[, 1:11]
for (i in names(isachs))
  levels(isachs[, i]) = c("LOW", "AVG", "HIGH")
fitted <- bn.fit(dag.sachs, isachs, method = "bayes")


###################################################
### code chunk number 33: chapter6.rnw:765-766
###################################################
options(digits=3)


###################################################
### code chunk number 34: chapter6.rnw:775-777
###################################################
library(gRain)
jtree <- compile(as.grain(fitted))


###################################################
### code chunk number 35: chapter6.rnw:782-783
###################################################
jlow <- setEvidence(jtree, nodes = "Erk", states  = "LOW")


###################################################
### code chunk number 36: chapter6.rnw:788-790
###################################################
querygrain(jtree, nodes = "Akt")$Akt
querygrain(jlow, nodes = "Akt")$Akt


###################################################
### code chunk number 37: chapter6.rnw:813-818
###################################################
causal.sachs <- drop.arc(dag.sachs, "PKA", "Erk")
causal.sachs <- drop.arc(causal.sachs, "Mek", "Erk")
cfitted <- bn.fit(causal.sachs, isachs, method = "bayes")
cjtree <- compile(as.grain(cfitted))
cjlow <- setEvidence(cjtree, nodes = "Erk", states  = "LOW")


###################################################
### code chunk number 38: chapter6.rnw:823-825
###################################################
querygrain(cjtree, nodes = "PKA")$PKA
querygrain(cjlow, nodes = "PKA")$PKA


###################################################
### code chunk number 39: chapter6.rnw:833-835
###################################################
querygrain(jtree, nodes = "PKA")$PKA
querygrain(jlow, nodes = "PKA")$PKA


###################################################
### code chunk number 40: chapter6.rnw:840-841
###################################################
names(which.max(querygrain(jlow, nodes = c("PKA"))$PKA))


###################################################
### code chunk number 41: bc01
###################################################
library(rbmn)
data(boco)
round(head(boco), 1)
boco$B <- boco$W / boco$H^2 * 10^4
dim(boco)
n <- nrow(boco)
vr <- colnames(boco)[5:13]
co <- c("A", "H", "W", "C", "B")


###################################################
### code chunk number 42: chapter6.rnw:998-999
###################################################
set.seed(1234)


###################################################
### code chunk number 43: bc02
###################################################
str <- sort(sample(n, round(n/2)))
dtr <- boco[str, ]
dva <- boco[-str, ]


###################################################
### code chunk number 44: bc02b
###################################################
satua <- lm(cbind(TF, LF, AF, TL, LL, AL, TB, LB, AB) ~ 
              A + H + W + C + B, data = dtr)
r.dof <- anova(satua)["Residuals", "Df"]
satup <- predict(satua, newdata = dva)
satabias <- abs(dva[, vr] - satup)
satstdev <- outer(rep(1, nrow(dtr)),
              sqrt(colSums(residuals(satua)^2)/r.dof), "*")
satsep <- sqrt(satabias^2 + satstdev^2)
satgsco <- cbind(colMeans(satabias), colMeans(satstdev),
                 colMeans(satsep))
colnames(satgsco) <- c("|Bias|", "Sd.Dev", "SEP")


###################################################
### code chunk number 45: chapter6.rnw:1048-1051
###################################################
round(satgsco, 2)
satsupe <- colSums(satgsco)
round(satsupe, 2)


###################################################
### code chunk number 46: bc03
###################################################
library(bnlearn)
dag1 <- hc(dtr)
paste(substr(modelstring(dag1), 1, 40), "...", sep = "")


###################################################
### code chunk number 47: bc04
###################################################
wl1 <- cbind(from = rep(co, each = 9), to = rep(vr, 5))
dag2 <- hc(dtr, whitelist = wl1)
paste(substr(modelstring(dag2), 1, 40), "...", sep = "")


###################################################
### code chunk number 48: bc05
###################################################
bl1 <- wl1[, 2:1]
dag3 <- hc(dtr, blacklist = bl1)
paste(substr(modelstring(dag3), 1, 40), "...", sep = "")
all.equal(dag2, dag3)


###################################################
### code chunk number 49: bc06
###################################################
iwl <- 1:15
wl2 <- wl1[iwl, ]
bl2 <- bl1[-iwl, ]
dag4 <- hc(dtr, whitelist = wl2, blacklist = bl2)
paste(substr(modelstring(dag4), 1, 40), "...", sep = "")


###################################################
### code chunk number 50: bc07
###################################################
bn2 <- bn.fit(dag2, data = dtr)


###################################################
### code chunk number 51: bc08
###################################################
library(rbmn)
mn2 <- gema2mn(nbn2gema(bnfit2nbn(bn2)))
bias <- stde <- dva[, vr]
for (ind in 1:nrow(dva)) {
  mni <- condi4joint(mn2, par = vr, pour = co, 
                       unlist(dva[ind, co]))
  bias[ind, vr] <- dva[ind, vr] - mni$mu[vr]
  stde[ind, vr] <- sqrt(diag(mni$gamma)[vr])
}#FOR
sep <- sqrt(bias^2 + stde^2)


###################################################
### code chunk number 52: bc09
###################################################
gscores <- cbind(colMeans(abs(bias)), colMeans(stde),
                 colMeans(sep))
colnames(gscores) <- c("|Bias|", "Sd.Dev", "SEP")
round(gscores, 2)
superf <- colSums(gscores)
round(superf, 2)


###################################################
### code chunk number 53: bc10
###################################################
library(graph)
library(igraph)
load("bc.poco.rda")
cbind(posi, colo)[c(1:3, 10:11), ]
idag2 <- igraph.from.graphNEL(as.graphNEL(dag2))
nad <- V(idag2)$label <- V(idag2)$name
edcol <- rep("lightgrey", nrow(arcs(dag2)))
aa <- which((arcs(dag2)[, 1] %in% vr) & 
            (arcs(dag2)[, 2] %in% vr))
va <- as.numeric(E(idag2, P = t(arcs(dag2)[aa, ])))
edcol[va] <- "black"
plot(idag2, layout = posi[nad, ], main = "DAG 2",
       edge.color = edcol, vertex.color = colo[nad])


###################################################
### code chunk number 54: bc11
###################################################
av.tl <- anova(lm(TL ~ H + W, data = dva))
1 - av.tl["Residuals", "Sum Sq"] / sum(av.tl[, "Sum Sq"])


