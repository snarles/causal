# Data files are stored as Excel Spreadsheets with columns as indicated in the headers.  The stimulations used are as indicated in the Materials and Methods.
# 
# 
# For the main result, the following conditions were used:
#   
#   Stimulation:    				        File                        .
# 1. CD3, CD28						          cd3cd28.xls
# 2. CD3, CD28, ICAM-2				      cd3cd28icam2.xls
# 3. CD3, CD28, akt-inhibitor				cd3cd28+aktinhib.xls
# 4. CD3, CD28, G0076				        cd3cd28+g0076.xls
# 5. CD3, CD28, Psitectorigenin			cd3cd28+psitect.xls
# 6. CD3, CD28, U0126				        cd3cd28+u0126.xls
# 7. CD3, CD28, LY294002				    cd3cd28+ly.xls
# 8. PMA							              pma.xls
# 9. beta2camp						          b2camp.xls
# 
# 
# 
# For the simulated western blots, the following were ALSO used, in addition to the 9 above:
#   
#   Stimulation:					                File                        .
# 10. CD3, CD28, ICAM-2, akt-inhib	     	cd3cd28icam2+aktinhib.xls
# 11. CD3, CD28, ICAM-2, G0076			      cd3cd28icam2+g0076.xls
# 12. CD3, CD28, ICAM-2, Psitectorigenin 	cd3cd28icam2+psit.xls
# 13. CD3, CD28, ICAM-2, U0126			      cd3cd28icam2+u0126.xls
# 14. CD3, CD28, ICAM-2, LY294002		      cd3cd28icam2+ly.xls

####
##  Load Sachs data
####

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

head(sachs0)
head(sachs_ints[[1]])

nms <- c("Raf", "Mek", "PLCg", "PIP2", "PIP3", "Erk", "Akt", "PKA", "PKC", "p38", "Jnk")
colnames(sachs0) <- nms
for (i in 1:9) colnames(sachs_ints[[i]]) <- nms

###
#  Marginal histograms
###

sachs <- c(list(sachs0), sachs_ints)
layout(matrix(1:12, 4, 3))
par("mar") # 5.1 4.1 4.1 2.1
par(mar = c(2, 2, 2, 2))
for (i in 1:10) hist(log(sachs[[i]]$Raf), main = paste(i))
layout(matrix(1:12, 4, 3))
for (i in 1:10) hist(log(sachs[[i]]$Mek), main = paste(i))
layout(matrix(1:12, 4, 3))
for (i in 1:10) hist(log(sachs[[i]]$PLCg), main = paste(i))

sachs_all <- do.call(rbind, sachs)


## Pairs plot
png("images/sachs_analysis_pairs2.png", width=1000, height=1000)
pairs(log(sachs_all), pch = ".", col = rgb(0,0,0,0.1))
dev.off()

####
##  DIAGRAMS GROUP 1: Effect of interventions
####


###
# Raf vs Erk
# C1 == C3 != C6
# Not a good example since Raf also changes
###
layout(1)
yl <- c(0, 6); xl <- c(0, 9)
par(bg = "white")
plot1(log(Erk) ~ log(Raf), sachs_ints[[1]], main = "Control")
plot1(log(Erk) ~ log(Raf), sachs_ints[[3]], main = "Akt-")
plot1(log(Erk) ~ log(Raf), sachs_ints[[6]], main = "Mek-")

triplot(log(Erk) ~ log(Raf), sachs_ints[c(1,3,6)])


###
# PIP3 vs PIP2
# C1 == C3 != C5
# Ok!
###

yl <- c(0, 8); xl <- c(0, 8)
par(bg = "white")
fmla <- log(PIP2) ~ log(PIP3)
plot1(fmla, sachs_ints[[1]], main = "Control")
plot1(fmla, sachs_ints[[3]], main = "Akt-")
plot1(fmla, sachs_ints[[5]], main = "PIP2-")

triplot(fmla, sachs_ints[c(1,3,5)])

###
#  PKC vs Akt
#  C0 =? C8 =! C9
#  OK!
###

yl <- c(0, 9); xl <- c(0, 7)
par(bg = "white")
fmla <- log(Akt) ~ log(PKC)
plot1(fmla, sachs[[1+0]], main = "Control")
plot1(fmla, sachs[[1+8]], main = "PKC+")
plot1(fmla, sachs[[1+9]], main = "PKA+")

triplot(fmla, sachs[1 + c(0,8,9)])


####
##  DIAGRAMS: CONDITIONAL INDEPENDENCE PLOTS
####

###
#  p38 Jnk | PKA, PKC
###

zattach(sachs_all)
par(bg = "white")
fm <- log(p38) ~ log(Jnk); control <- c("PKA", "PKC"); cf <- log(PKA) ~ log(PKC)
layout(1); par(mar = c(5.1, 4.1, 4.1, 2.1))

layout(1)
marginal_plot()
## subplots
layout(matrix(1:3, 1, 3))
w1 <- getwind()
condition_w(w1)
condition_w(w1, marginal = TRUE)

####
##  DIAGRAMS: INVARIANCE
####

###
# p38 | PKC, PKA invariant under perturbation C1 vs C2
# Nope, not invariant!
###

layout(1)
fmla <- log(p38) ~ log(PKC) + log(PKA)
dat1 <- sachs_ints[[6]]
dat2 <- sachs_ints[[9]]
yX1 <- model.frame(fmla, data = dat1)
yX2 <- model.frame(fmla, data = dat2)
y1 <- yX1[, 1]
y2 <- yX2[, 1]
X1 <- yX1[, -1]
X2 <- yX2[, -1]
Xa1 <- data.frame(cbind(y = y1, quadd(X1)))
Xa2 <- data.frame(cbind(y = y2, quadd(X2)))
res1 <- lm(y ~ ., data = Xa1)
res2 <- lm(y ~ ., data = Xa2)
yh11 <- predict(res1, Xa1)
yh21 <- predict(res1, Xa2)
yh12 <- predict(res2, Xa1)
yh22 <- predict(res2, Xa2)

layout(matrix(1:2, 1, 2))
plot(yh11, yh12); abline(0, 1, col = "red")
plot(yh21, yh22); abline(0, 1, col = "red")

