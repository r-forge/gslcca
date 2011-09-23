library(gslcca)
clonidine <- readSpectra("Clonidine Light.txt", info = 1:3, 
                    treatment = c("Control", "Low Dose", 
                                  "Middle Dose", "High Dose"))
 
## subsample every 10 minutes
select <- rep(c(TRUE, FALSE), length.out = nrow(clonidine$spectra))
clonidine <- clonidine[select,]

save("clonidine", file = "clonidine.RData")
library(tools)
checkRdaFiles("clonidine.RData")

resaveRdaFiles("clonidine.RData")
checkRdaFiles("clonidine.RData")
