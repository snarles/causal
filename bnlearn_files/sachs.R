# analyzing the Sachs data with bnlearn.
library(bnlearn)

sachs = read.table("bnlearn_files/sachs.data.txt", header = TRUE)
dsachs = discretize(sachs, method = "hartemink", breaks = 3, ibreaks = 60)
boot = boot.strength(data = dsachs, R = 500, algorithm = "hc",
         algorithm.args = list(score = "bde", iss = 10))
boot[(boot$strength > 0.85) & (boot$direction >= 0.5), ]
avg.boot = averaged.network(boot, threshold = 0.85)

nodes = names(dsachs)
start = random.graph(nodes = nodes, method = "melancon", num = 500)
netlist = lapply(start, function(net) {
  hc(dsachs, score = "bde", iss = 10, start = net) })
rnd = custom.strength(netlist, nodes = nodes)
rnd[(rnd$strength > 0.85) & (rnd$direction >= 0.5), ]
avg.start = averaged.network(rnd, threshold = 0.85)

all.equal(cpdag(avg.boot), cpdag(avg.start))
all.equal(moral(avg.boot), moral(avg.start))
score(avg.boot, dsachs, type = "bde", iss = 10)
score(avg.start, dsachs, type = "bde", iss = 10)

# analyzing the Sachs data with catnet.
library(catnet)

netlist = vector(500, mode = "list")
ndata = nrow(sachs)
netlist = lapply(netlist, function(net) {
  boot = dsachs[sample(ndata, replace = TRUE), ]
  nets = cnSearchOrder(boot)
  best = cnFindBIC(nets, ndata)
  cnMatEdges(best)
})
sa = custom.strength(netlist, nodes = nodes)
sa[(sa$strength > 0.85) & (sa$direction >= 0.5), ]
avg.catnet = averaged.network(sa, threshold = 0.85)
avg.catnet

# estimated threshold.
all.equal(avg.boot, averaged.network(boot, threshold = 0.60))
all.equal(avg.start, averaged.network(rnd, threshold = 0.60))
all.equal(avg.catnet, averaged.network(sa, threshold = 0.60))

averaged.network(boot)

# analysis with the auto-tuned threshold from the AIME paper.
isachs = read.table("bnlearn_files/sachs.interventional.txt", header = TRUE, colClasses = "factor")

par(mfrow = c(1, 2))

# treatment modeled with all edges in a whitelist.
wh = matrix(c(rep("INT", 11), names(isachs)[1:11]), ncol = 2)
bn.wh = tabu(isachs, whitelist = wh, score = "bde", iss = 10, tabu = 50)

# treatment modeled with the edges chosen by the model.
tiers = list("INT", names(isachs)[1:11])
bl = tiers2blacklist(tiers)
bn.tiers = tabu(isachs, blacklist = bl, score = "bde", iss = 10, tabu = 50)


highlight = list(arcs = wh, nodes = "INT", col = "grey")
graphviz.plot(bn.wh, highlight = highlight)
graphviz.plot(bn.tiers, highlight = highlight)


boot = boot.strength(data = isachs, R = 500, algorithm = "tabu",
         algorithm.args = list(score = "bde", iss = 10, whitelist = wh,  tabu = 50))
bn.wh = averaged.network(boot)

boot = boot.strength(data = isachs, R = 500, algorithm = "tabu",
         algorithm.args = list(score = "bde", iss = 10, blacklist = bl,  tabu = 50))
bn.tiers = averaged.network(boot)

highlight = list(arcs = wh, nodes = "INT", col = "grey")
graphviz.plot(bn.wh, highlight = highlight)
graphviz.plot(bn.tiers, highlight = highlight)

set.seed(123)
nodes = names(isachs)[1:11]
start = random.graph(nodes = nodes, method = "melancon", num = 500, burn.in = 10^5, every = 100)
perturbed = sapply(1:11, function(x) which(isachs$INT == x))
names(perturbed) = nodes
netlist = lapply(start, function(net) {
  tabu(isachs[, 1:11], score = "mbde", exp = perturbed, iss = 10, start = net, tabu = 50) })
arcs = custom.strength(netlist, nodes = nodes)

graphviz.plot(averaged.network(arcs, threshold = 0.85))
shd(true, averaged.network(arcs, threshold = 0.85))
unlist(compare(true, averaged.network(arcs, threshold = 0.85)))


