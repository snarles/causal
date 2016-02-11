### R code from vignette source '/home/fizban/book-pratiqueR/english/chap5/chapter5.rnw'

###################################################
### code chunk number 1: chapter5.rnw:146-151
###################################################
library(bnlearn)
data(marks)
latent <- factor(c(rep("A", 44), "B",
                  rep("A", 7), rep("B", 36)))
marks$LAT <- latent


###################################################
### code chunk number 2: chapter5.rnw:156-158
###################################################
library(deal)
net <- network(marks)


###################################################
### code chunk number 3: chapter5.rnw:161-162
###################################################
net


###################################################
### code chunk number 4: chapter5.rnw:166-167
###################################################
prior <- jointprior(net, N = 5)


###################################################
### code chunk number 5: chapter5.rnw:171-173
###################################################
net <- learn(net, marks, prior)$nw
.z <- capture.output({best <- autosearch(net, marks, prior)})


###################################################
### code chunk number 6: chapter5.rnw:175-177 (eval = FALSE)
###################################################
## net <- learn(net, marks, prior)$nw
## best <- autosearch(net, marks, prior)


###################################################
### code chunk number 7: chapter5.rnw:196-197
###################################################
mstring <- deal::modelstring(best$nw)


###################################################
### code chunk number 8: chapter5.rnw:226-230
###################################################
dag.bnlearn <- model2network(
    "[ANL][MECH][LAT|ANL:MECH][VECT|LAT][ALG|LAT][STAT|LAT]")
dag.deal <- model2network(mstring)
unlist(bnlearn::compare(cpdag(dag.deal), cpdag(dag.bnlearn)))


###################################################
### code chunk number 9: chapter5.rnw:269-273
###################################################
library(catnet)
dmarks <- discretize(marks, breaks = 2, method = "interval")
ord <- cnSearchSA(dmarks, maxParentSet = 2)
ord


###################################################
### code chunk number 10: chapter5.rnw:279-281
###################################################
nets <- ord@nets
nets[[1]]


###################################################
### code chunk number 11: chapter5.rnw:286-288
###################################################
best <- cnFindBIC(ord, nrow(dmarks))
best


###################################################
### code chunk number 12: chapter5.rnw:296-297
###################################################
cnSamples(best, numsamples = 4)


###################################################
### code chunk number 13: chapter5.rnw:301-303
###################################################
em <- empty.graph(names(dmarks))
arcs(em) <- cnMatEdges(best)


###################################################
### code chunk number 14: chapter5.rnw:338-341
###################################################
library(pcalg)
marks <- marks[, colnames(marks) != "LAT"]
suffStat <- list(C = cor(marks), n = nrow(marks))


###################################################
### code chunk number 15: chapter5.rnw:345-347
###################################################
pc.fit <- pc(suffStat, indepTest = gaussCItest,
            p = ncol(marks), alpha = 0.05)


###################################################
### code chunk number 16: chapter5.rnw:352-353
###################################################
pc.fit@graph


###################################################
### code chunk number 17: chapter5.rnw:363-365 (eval = FALSE)
###################################################
## fci.fit <- fci(suffStat, indepTest = gaussCItest,
##              p = ncol(marks), alpha = 0.05)


