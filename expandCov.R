rm(list = ls(all = T))
library(dplyr)

source("functions.R")
source("cvInitNG.R")
source("expandFunctions.R")

call <- ~ EVI + ALT + POPDEN + LSTN + RFALL + RDIST + WDIST + NLIGHT + LATNUM + LONGNUM

covLab <- attr(terms(call), "term.labels")

# newLdat$cluster <- ldat2$cluster  

newLdat <- model.matrix(call, data = ldat2) %>% 
  .[,-1] %>% 
  as.data.frame() %>%
  addQuadTerm(covLab) %>%
  addInterTerm(covLab) %>%
  addCubiTerm(covLab) %>%
  addCubiInterTerm(covLab) %>%
  scale() %>%
  as.data.frame() %>%
  removeHiCor(.90)
  
newLdat$cluster <- ldat2$cluster

write.csv(newLdat, "Data/NG/NGLocDatExpCubi.csv", row.names = F)


#################
rm(list = ls(all = T))
library(dplyr)

source("functions.R")
source("cvInitGH.R")
source("expandFunctions.R")

call <- ~ EVI + ALT + POPDEN + LSTD + LSTN + RFALL + RDIST + WDIST + NLIGHT + LATNUM + LONGNUM

covLab <- attr(terms(call), "term.labels")

# newLdat$cluster <- ldat2$cluster  

newLdat <- model.matrix(call, data = ldat2) %>% 
  .[,-1] %>% 
  as.data.frame() %>%
  addQuadTerm(covLab) %>%
  addInterTerm(covLab) %>%
  addCubiTerm(covLab) %>%
  addCubiInterTerm(covLab) %>%
  scale() %>%
  as.data.frame() %>%
  removeHiCor(.90)

newLdat$cluster <- ldat2$cluster

write.csv(newLdat, "Data/GH/GHLocDatExpCubi.csv", row.names = F)

#################
rm(list = ls(all = T))
library(dplyr)

source("functions.R")
source("cvInitTG.R")
source("expandFunctions.R")

call <- ~ EVI + ALT + POPDEN + LSTD + LSTN + RFALL + RDIST + WDIST + NLIGHT + LATNUM + LONGNUM

covLab <- attr(terms(call), "term.labels")

# newLdat$cluster <- ldat2$cluster  

newLdat <- model.matrix(call, data = ldat2) %>% 
  .[,-1] %>% 
  as.data.frame() %>%
  addQuadTerm(covLab) %>%
  addInterTerm(covLab) %>%
  addCubiTerm(covLab) %>%
  addCubiInterTerm(covLab) %>%
  scale() %>%
  as.data.frame() %>%
  removeHiCor(.90)

newLdat$cluster <- ldat2$cluster

write.csv(newLdat, "Data/TG/TGLocDatExpCubi.csv", row.names = F)
#################
rm(list = ls(all = T))
library(dplyr)

source("functions.R")
source("cvInitBF.R")
source("expandFunctions.R")

call <- ~ EVI + ALT + POPDEN + LSTD + LSTN + RFALL + RDIST + WDIST + NLIGHT +
  LATNUM + LONGNUM

covLab <- attr(terms(call), "term.labels")

newLdat <- model.matrix(call, data = ldat2) %>% 
  .[,-1] %>% 
  as.data.frame() %>%
  addQuadTerm(covLab) %>%
  addInterTerm(covLab) %>%
  addCubiTerm(covLab) %>%
  addCubiInterTerm(covLab) %>%
  scale() %>%
  as.data.frame() %>%
  removeHiCor(.90)

newLdat$cluster <- ldat2$cluster

write.csv(newLdat, "Data/BF/BFLocDatExpCubi.csv", row.names = F)
