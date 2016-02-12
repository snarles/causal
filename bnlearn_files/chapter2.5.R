####
##  2.5.1 Model Averaging
####

library(bnlearn)
sachs = read.table("bnlearn_files/sachs.data.txt", header = TRUE)
# View(sachs)
dsachs = discretize(sachs, method = "hartemink", breaks = 3, 
                    ibreaks = 60, idisc = "quantile")
boot = boot.strength(data = dsachs, R = 500, algorithm = "hc",
                     algorithm.args = list(score = "bde",
                                           iss = 10))
boot[(boot$strength > 0.85) & (boot$direction >= 0.5), ]
avg.boot = averaged.network(boot, threshold = 0.85)
avg.boot
graphviz.plot(avg.boot, layout = "fdp")
graphviz.plot(cpdag(avg.boot), layout = "fdp")


nodes = names(dsachs)
start = random.graph(nodes = nodes, method = "ic-dag", num = 500)
netlist = lapply(start, function(net) {
  hc(dsachs, score = "bde", iss = 10, start = net)
})
rnd = custom.strength(netlist, nodes = nodes)
rnd[(rnd$strength > 0.85) & (rnd$direction >= 0.5), ]
avg.start = averaged.network(rnd, threshold = 0.85)
all.equal(cpdag(avg.boot), cpdag(avg.start))

score(cextend(cpdag(avg.start)), dsachs, type = "bde", iss = 10)


library(catnet)
netlist = vector(500, mode = "list")
ndata = nrow(dsachs)
netlist = lapply(netlist, function(net) {
  boot = dsachs[sample(ndata, replace = TRUE), ]
  nets = cnSearchOrder(boot)
  best = cnFindBIC(nets, ndata)
  cnMatEdges(best)
})
sa = custom.strength(netlist, nodes = nodes)
sa[(sa$strength > 0.85) & (sa$direction >= 0.5), ]
avg.catnet = averaged.network(sa, threshold = 0.85)

####
##  2.5.2 Choosing the significance threshold
####

all.equal(averaged.network(boot, threshold = 0.50),
          averaged.network(boot, threshold = 0.70))
averaged.network(boot)

####
##  2.5.3 Handling interventional data
####

isachs = read.table("bnlearn_files/sachs.interventional.txt",
                    header = TRUE, colClasses = "factor")
View(isachs)
wh = matrix(c(rep("INT", 11), names(isachs)[1:11]), ncol = 2)
bn.wh = tabu(isachs, whitelist = wh, score = "bde", iss = 10, tabu = 50)
tiers = list("INT", names(isachs)[1:11])
bl = tiers2blacklist(nodes = tiers)
bn.tiers = tabu(isachs, blacklist = bl, score = "bde", iss = 10, tabu = 50)

INT = sapply(1:11, function(x) {which(isachs$INT == x)})
nodes = names(isachs)[1:11]
names(INT) = nodes
start = random.graph(nodes = nodes, method = "melancon",
                     num = 500, burn.in = 10^5, every = 100)
netlist = lapply(start, function(net) {
  tabu(isachs[, 1:11], score = "mbde", exp = INT,
       iss = 10, start = net, tabu = 50)  
})
arcs = custom.strength(netlist, nodes = nodes)
bn.mbde = averaged.network(arcs, threshold = 0.85)

graphviz.plot(bn.wh, layout = "dot")
graphviz.plot(bn.tiers, layout = "dot")
graphviz.plot(bn.mbde, layout = "dot")

####
##  My own extra stuff
####

isachs = read.table("bnlearn_files/sachs.interventional.txt",
                    header = TRUE, colClasses = "factor")

