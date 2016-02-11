#### Copied from the 2003 book

library(bnlearn)
data(marks)
str(marks)
ug = empty.graph(names(marks))
## "The following generates an error"
#arcs(ug, ignore.cyles = TRUE) <- matrix(
#   c("MECH", "VECT", "MECH", "ALG", "VECT", "MECH",
#     "VECT", "ALG", "ALG", "MECH", "ALG", "VECT",
#     "ALG", "ANL", "ALG", "STAT", "ANL", "ALG",
#     "ANL", "STAT", "STAT", "ALG", "STAT", "ANL"),
#     ncol=2, byrow = TRUE,
#     dimnames = list(c(), c("from", "to")))

dag = empty.graph(names(marks))
arcs(dag) <- matrix(
  c("VECT", "MECH", "ALG", "MECH", "ALG", "VECT",
    "ANL", "ALG", "STAT", "ALG", "STAT", "ANL"),
  ncol=2, byrow = TRUE,
  dimnames = list(c(), c("from", "to")))
dag
plot(dag)
node.ordering(dag)
score(dag, data = marks, type = "loglik-g")
dag.eq = reverse.arc(dag, "STAT", "ANL")
score(dag.eq, data = marks, type = "loglik-g")

vstructs(dag, moral=FALSE)
vstructs(dag.eq, moral=FALSE)

all.equal(cpdag(dag), cpdag(dag.eq))
all.equal(moral(dag), moral(dag.eq))

dag2 = drop.arc(dag, from = "STAT", to = "ANL")
dag3 = drop.arc(dag, from = "ALG", to = "VECT")
vstructs(dag2, moral = FALSE)
