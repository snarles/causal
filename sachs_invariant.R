####
##  Apply InvariantCausalPrediction to Sachs data
####

library(InvariantCausalPrediction)


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

nms <- c("Raf", "Mek", "PLCg", "PIP2", "PIP3", "Erk", "Akt", "PKA", "PKC", "p38", "Jnk")
colnames(sachs0) <- nms
for (i in 1:9) colnames(sachs_ints[[i]]) <- nms

for (i in 1:9) sachs_ints[[i]] <- cbind(ExpInd = i, sachs_ints[[i]])
sachs0 <- cbind(ExpInd = 0, sachs0)
sachs <- c(list(sachs0), sachs_ints)
sachs_all <- do.call(rbind, sachs)

ExpInd <- sachs_all[, "ExpInd"]
sachs_data <- sachs_all[, -1]

####
##  Which interventions include which variables??
####

(prots <- colnames(sachs_data))

interventions <- list()
glob <- c("PIP2", "PIP3", "Raf")
interventions[[1]] <- glob
interventions[[2]] <- glob
interventions[[3]] <- glob
interventions[[4]] <- c(glob, "PKC")
interventions[[5]] <- c(glob, "PIP2")
interventions[[6]] <- c(glob, "Erk")
interventions[[7]] <- glob
interventions[[8]] <- "PKC"
interventions[[9]] <- "PKA"

intmat <- matrix(FALSE, 9, 11)
colnames(intmat) <- prots
for (i in 1:11) {
  pp <- prots[i]
  for (j in 1:9) {
    if (pp %in% interventions[[j]]) intmat[j, i] <- TRUE
  }
}

####
##  Apply method
####

apply_ICP <- function(prot_no, hidden = FALSE, ...) {
  pp <- prots[prot_no]
  valid_exp <- c(0, which(!intmat[, pp]))
  exp_filt <- ExpInd %in% valid_exp
  XY <- sachs_data[exp_filt, ]
  EI <- ExpInd[exp_filt]
  ##unique(EI)
  X <- as.matrix(XY[, -prot_no])
  Y <- XY[, prot_no]
  if (hidden) {
    resICP <- hiddenICP(X=X, Y=Y, ExpInd=EI, ...)        
  } else {
    resICP <- ICP(X=X, Y=Y, ExpInd=EI, ...)    
  }
  resICP  
}

## No hidden
apply_ICP(1, FALSE, alpha = 0.1)
apply_ICP(2, FALSE, alpha = 0.1)
apply_ICP(3, FALSE, alpha = 0.1)
apply_ICP(4, FALSE, alpha = 0.1) ## PIP2 ~ PLCg
apply_ICP(5, FALSE, alpha = 0.1)
apply_ICP(6, FALSE, alpha = 0.1)
apply_ICP(7, FALSE, alpha = 0.1)
apply_ICP(8, FALSE, alpha = 0.1)
apply_ICP(9, FALSE, alpha = 0.1)
apply_ICP(10, FALSE, alpha = 0.1)
apply_ICP(11, FALSE, alpha = 0.1)

## Hidden
apply_ICP(1, TRUE, alpha = 0.1)
apply_ICP(2, TRUE, alpha = 0.1)
apply_ICP(3, TRUE, alpha = 0.1)
apply_ICP(4, TRUE, alpha = 0.1) ## PIP2 ~ PLCg
apply_ICP(5, TRUE, alpha = 0.1)
apply_ICP(6, TRUE, alpha = 0.1)
apply_ICP(7, TRUE, alpha = 0.1)
apply_ICP(8, TRUE, alpha = 0.1)
apply_ICP(9, TRUE, alpha = 0.1)
apply_ICP(10, TRUE, alpha = 0.1)
apply_ICP(11, TRUE, alpha = 0.1)
