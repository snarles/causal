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

library(readxl)

sachs0 = read.table("bnlearn_files/sachs.data.txt", header = TRUE)
sachs_ints <- list()
sachs_ints[[1]] <- readxl::read_excel("cd3cd28.xls")
sachs_ints[[2]] <- readxl::read_excel("cd3cd28icam2.xls")
sachs_ints[[3]] <- readxl::read_excel("cd3cd28+aktinhib.xls")
sachs_ints[[4]] <- readxl::read_excel("cd3cd28+g0076.xls")
sachs_ints[[5]] <- readxl::read_excel("cd3cd28+psitect.xls")
sachs_ints[[6]] <- readxl::read_excel("cd3cd28+u0126.xls")
sachs_ints[[7]] <- readxl::read_excel("cd3cd28+ly.xls")
sachs_ints[[8]] <- readxl::read_excel("pma.xls")
sachs_ints[[9]] <- readxl::read_excel("b2camp.xls")

# 1. CD3, CD28  					          cd3cd28.xls
# 2. CD3, CD28, ICAM-2				      cd3cd28icam2.xls
# 3. CD3, CD28, akt-inhibitor				cd3cd28+aktinhib.xls
# 4. CD3, CD28, G0076				        cd3cd28+g0076.xls
# 5. CD3, CD28, Psitectorigenin			cd3cd28+psitect.xls
# 6. CD3, CD28, U0126				        cd3cd28+u0126.xls
# 7. CD3, CD28, LY294002				    cd3cd28+ly.xls
# 8. PMA							              pma.xls
# 9. beta2camp						          b2camp.xls
