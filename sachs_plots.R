# PLOTS FOR PRESENTATION
# See sachs_analysis.R for motivation 

## SETUP

source("source.R")
source("analysis_source.R")
library(readxl)
sachs0 = read.table("bnlearn_files/sachs.data.txt", header = TRUE)
sachs_ints <- list()
sachs_ints[[1]] <- readxl::read_excel("sachs_data/1. cd3cd28.xls")
sachs_ints[[2]] <- readxl::read_excel("sachs_data/2. cd3cd28icam2.xls")
sachs_ints[[3]] <- readxl::read_excel("sachs_data/3. cd3cd28+aktinhib.xls")
sachs_ints[[4]] <- readxl::read_excel("sachs_data/4. cd3cd28+g0076.xls")
sachs_ints[[5]] <- readxl::read_excel("sachs_data/5. cd3cd28+psitect.xls")
sachs_ints[[6]] <- readxl::read_excel("sachs_data/6. cd3cd28+u0126.xls")
sachs_ints[[7]] <- readxl::read_excel("sachs_data/7. cd3cd28+ly.xls")
sachs_ints[[8]] <- readxl::read_excel("sachs_data/8. pma.xls")
sachs_ints[[9]] <- readxl::read_excel("sachs_data/9. b2camp.xls")
nms <- c("Raf", "Mek", "PLCg", "PIP2", "PIP3", "Erk", "Akt", "PKA", "PKC", "p38", "Jnk")
colnames(sachs0) <- nms
for (i in 1:9) colnames(sachs_ints[[i]]) <- nms
sachs <- c(list(sachs0), sachs_ints)
sachs_all <- do.call(rbind, sachs)

####
##  DIAGRAMS GROUP 1: Effect of interventions
####

###
# PIP3 vs PIP2
# C1 == C3 != C5
# Ok!
###

pdfcmd <- function(fname) {
  pdf(paste0("images/", fname), width = 4.5, height = 4.5)
}

pngcmd <- function(fname) {
  png(paste0("images/", fname), width = 900, height = 500)
}

yl <- c(0, 8); xl <- c(0, 8)
par(bg = "white")
fmla <- log(PIP2) ~ log(PIP3)
pdfcmd("plot01_01.pdf")
plot1(fmla, sachs_ints[[1]], main = "Control")
dev.off()
pdfcmd("plot01_02.pdf")
plot1(fmla, sachs_ints[[3]], main = "Akt-")
dev.off()
pdfcmd("plot01_03.pdf")
plot1(fmla, sachs_ints[[5]], main = "PIP2-")
dev.off()
pdfcmd("plot01_04.pdf")
triplot(fmla, sachs_ints[c(1,3,5)])
dev.off()

###
#  PKC vs Akt
#  C0 =? C8 =! C9
#  OK!
###

yl <- c(0, 9); xl <- c(0, 7)
par(bg = "white")
fmla <- log(Akt) ~ log(PKC)
pdfcmd("plot02_01.pdf")
plot1(fmla, sachs[[1+0]], main = "Control")
dev.off()
pdfcmd("plot02_02.pdf")
plot1(fmla, sachs[[1+8]], main = "PKC+")
dev.off()
pdfcmd("plot02_03.pdf")
plot1(fmla, sachs[[1+9]], main = "PKA+")
dev.off()
pdfcmd("plot02_04.pdf")
triplot(fmla, sachs[1 + c(0,8,9)])
dev.off()

###
#  p38 Jnk | PKA, PKC
###

pdfcmd2 <- function(fname) {
  pdf(paste0("images/", fname), width = 6, height = 3)
}

zattach(sachs_all)
par(bg = "white")
fm <- log(p38) ~ log(Jnk); control <- c("PKA", "PKC"); cf <- log(PKA) ~ log(PKC)
layout(1); par(mar = c(5.1, 4.1, 4.1, 2.1))

layout(1)
pdfcmd("plot03_01.pdf")
marginal_plot()
dev.off()




## getting w1s
# w1 <- getwind()
# condition_w(w1)
# condition_w(w1, marginal = TRUE)
# printw(w1)

w1 <- list(x=c(1.44902468836522,2.07741413494026),y=c(6.85303764838221,7.22635573516386)) # 77
pdfcmd2("plot03_02.pdf")
par(mar = c(5, 4, 2, 2))
layout(matrix(1:3, 1, 3))
plotwind(w1)
condition_w(w1)
condition_w(w1, marginal = TRUE)
dev.off()

w1 <- list(x=c(2.52935390244464,3.22162907609093),y=c(6.28237064068292,6.54201476660439)) # 362
pdfcmd2("plot03_03.pdf")
par(mar = c(5, 4, 2, 2))
layout(matrix(1:3, 1, 3))
plotwind(w1)
condition_w(w1)
condition_w(w1, marginal = TRUE)
dev.off()

w1 <- list(x=c(4.12691199547453,4.87243910555514),y=c(2.24346201523775,2.50310614115923)) # 49
pdfcmd2("plot03_04.pdf")
par(mar = c(5, 4, 2, 2))
layout(matrix(1:3, 1, 3))
plotwind(w1)
condition_w(w1)
condition_w(w1, marginal = TRUE)
dev.off()



####
##  DIAGRAMS: INVARIANCE
####

pdfcmd3 <- function(fname) {
  pdf(paste0("images/", fname), width = 6, height = 3)
}

plot2 <- function(...) {plot(pch = 19, cex = 0.4, ...)}

###
# Akt | PKA, Erk invariant under perturbation C1 vs C2 vs C3
# Akt | PKA, Erk NOT invariant under perturbation C1 vs C6
###

layout(1)
fmla <- log(Akt) ~ log(PKA) + log(Erk)
dat1 <- subsample(sachs_ints[[1]], 1)
dat2 <- subsample(sachs_ints[[9]], 1)
# dat1 <- subsample(sachs_ints[[1]], 0.5)
# dat2 <- subsample(sachs_ints[[3]], 0.5)
yX1 <- model.frame(fmla, data = dat1)
yX2 <- model.frame(fmla, data = dat2)
y1 <- yX1[, 1]
y2 <- yX2[, 1]
X1 <- yX1[, -1]
X2 <- yX2[, -1]
# Xa1 <- data.frame(cbind(y = y1, quadd(X1)))
# Xa2 <- data.frame(cbind(y = y2, quadd(X2)))
Xa1 <- data.frame(cbind(y = y1, X1))
Xa2 <- data.frame(cbind(y = y2, X2))
res1 <- lm(y ~ ., data = Xa1)
res2 <- lm(y ~ ., data = Xa2)
yh11 <- predict(res1, Xa1)
yh21 <- predict(res1, Xa2)
yh12 <- predict(res2, Xa1)
yh22 <- predict(res2, Xa2)

r11 <- y1 - yh11; r12 <- y1 - yh12
r21 <- y2 - yh21; r22 <- y2 - yh22

pdfcmd3("plot04_01.pdf")
layout(matrix(1:2, 1, 2))
plot2(yh11, yh12, main = "Control"); abline(0, 1, col = "red")
plot2(yh21, yh22, main = "PKA+"); abline(0, 1, col = "red")
dev.off()

pdfcmd3("plot04_02.pdf")
layout(matrix(1:2, 1, 2))
plot2(r11, r12, main = "Control"); abline(0, 1, col = "red")
plot2(r21, r22, main = "PKA+"); abline(0, 1, col = "red")
dev.off()

for (i in 1:3) {
  str1 <- paste0("plot04_01r", i, ".pdf")
  str2 <- paste0("plot04_02r", i, ".pdf")
  layout(1)
  dat1 <- subsample(sachs_ints[[1]], 0.5)
  dat2 <- subsample(sachs_ints[[9]], 0.5)
  yX1 <- model.frame(fmla, data = dat1)
  yX2 <- model.frame(fmla, data = dat2)
  y1 <- yX1[, 1]
  y2 <- yX2[, 1]
  X1 <- yX1[, -1]
  X2 <- yX2[, -1]
  # Xa1 <- data.frame(cbind(y = y1, quadd(X1)))
  # Xa2 <- data.frame(cbind(y = y2, quadd(X2)))
  Xa1 <- data.frame(cbind(y = y1, X1))
  Xa2 <- data.frame(cbind(y = y2, X2))
  res1 <- lm(y ~ ., data = Xa1)
  res2 <- lm(y ~ ., data = Xa2)
  yh11 <- predict(res1, Xa1)
  yh21 <- predict(res1, Xa2)
  yh12 <- predict(res2, Xa1)
  yh22 <- predict(res2, Xa2)
  r11 <- y1 - yh11; r12 <- y1 - yh12
  r21 <- y2 - yh21; r22 <- y2 - yh22
  pdfcmd3(str1)
  layout(matrix(1:2, 1, 2))
  plot2(yh11, yh12); abline(0, 1, col = "red")
  plot2(yh21, yh22); abline(0, 1, col = "red")
  dev.off()
  pdfcmd3(str2)
  layout(matrix(1:2, 1, 2))
  plot2(r11, r12); abline(0, 1, col = "red")
  plot2(r21, r22); abline(0, 1, col = "red")
  dev.off()
}




###
# Mek | Raf, Erk, PKA NOT invariant under perturbation C2, C4
###

layout(1)
fmla <- log(Raf) ~ log(Mek) + log(Erk) + log(PKA)
dat1 <- subsample(sachs_ints[[2]], 1)
dat2 <- subsample(sachs_ints[[4]], 1)
yX1 <- model.frame(fmla, data = dat1)
yX2 <- model.frame(fmla, data = dat2)
y1 <- yX1[, 1]
y2 <- yX2[, 1]
X1 <- yX1[, -1]
X2 <- yX2[, -1]
# Xa1 <- data.frame(cbind(y = y1, quadd(X1)))
# Xa2 <- data.frame(cbind(y = y2, quadd(X2)))
Xa1 <- data.frame(cbind(y = y1, X1))
Xa2 <- data.frame(cbind(y = y2, X2))
res1 <- lm(y ~ ., data = Xa1)
res2 <- lm(y ~ ., data = Xa2)
yh11 <- predict(res1, Xa1)
yh21 <- predict(res1, Xa2)
yh12 <- predict(res2, Xa1)
yh22 <- predict(res2, Xa2)

r11 <- y1 - yh11; r12 <- y1 - yh12
r21 <- y2 - yh21; r22 <- y2 - yh22

pdfcmd3("plot05_01.pdf")
layout(matrix(1:2, 1, 2))
plot2(yh11, yh12, main = "Control"); abline(0, 1, col = "red")
plot2(yh21, yh22, main = "PKC-"); abline(0, 1, col = "red")
dev.off()

pdfcmd3("plot05_02.pdf")
layout(matrix(1:2, 1, 2))
plot2(r11, r12, main = "Control"); abline(0, 1, col = "red")
plot2(r21, r22, main = "PKC-"); abline(0, 1, col = "red")
dev.off()

for (i in 1:3) {
  str1 <- paste0("plot05_01r", i, ".pdf")
  str2 <- paste0("plot05_02r", i, ".pdf")
  layout(1)
  dat1 <- subsample(sachs_ints[[2]], .5)
  dat2 <- subsample(sachs_ints[[4]], .5)
  yX1 <- model.frame(fmla, data = dat1)
  yX2 <- model.frame(fmla, data = dat2)
  y1 <- yX1[, 1]
  y2 <- yX2[, 1]
  X1 <- yX1[, -1]
  X2 <- yX2[, -1]
  # Xa1 <- data.frame(cbind(y = y1, quadd(X1)))
  # Xa2 <- data.frame(cbind(y = y2, quadd(X2)))
  Xa1 <- data.frame(cbind(y = y1, X1))
  Xa2 <- data.frame(cbind(y = y2, X2))
  res1 <- lm(y ~ ., data = Xa1)
  res2 <- lm(y ~ ., data = Xa2)
  yh11 <- predict(res1, Xa1)
  yh21 <- predict(res1, Xa2)
  yh12 <- predict(res2, Xa1)
  yh22 <- predict(res2, Xa2)
  r11 <- y1 - yh11; r12 <- y1 - yh12
  r21 <- y2 - yh21; r22 <- y2 - yh22
  pdfcmd3(str1)
  layout(matrix(1:2, 1, 2))
  plot2(yh11, yh12); abline(0, 1, col = "red")
  plot2(yh21, yh22); abline(0, 1, col = "red")
  dev.off()
  pdfcmd3(str2)
  layout(matrix(1:2, 1, 2))
  plot2(r11, r12); abline(0, 1, col = "red")
  plot2(r21, r22); abline(0, 1, col = "red")
  dev.off()
}
