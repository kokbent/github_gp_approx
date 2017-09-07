# file.copy("../Geospatial/DHSGH2014/GHTime.csv", "Data/GH/")

source("cvFunctions.R")

# Initialize for Cross Validations
malaWeal <- read.csv("Data/GH/DHS2014MalaEduWeal.csv") %>% 
  dplyr::select(cluster, hh, micro_mala, wealth, Edu = maxEdu, age)

malaWeal$age <- classifyAge(malaWeal$age)

wdat <- malaWeal
ldat <- read.csv("Data/GH/DHS2014EnvAnalysis.csv") %>% dplyr::select(-LC) %>%
  left_join(read.csv("Data/GH/GHTime.csv"), by = c("DHSCLUST" = "cluster"))

colnames(ldat)[6] <- "POPDEN"
ldat$GRUMP[ldat$GRUMP == 2] <- 0
ldat$GRUMP2 <- ldat$GRUMP
medUrbPop <- aggregate(ldat$POPDEN ~ ldat$GRUMP, FUN = median)[1,2]
ldat$GRUMP2[ldat$POPDEN < medUrbPop & ldat$GRUMP == 0] <- 0.5
ldat$LCType <- as.factor(ldat$LCType)
ldat$GRUMP2 <- as.factor(ldat$GRUMP2)

latlong <- ldat[,1:3]
ldat[,2:3] <- scale(ldat[,2:3])
ldat$meanDate <- scale(ldat$meanDate) %>% as.vector()
ind <- which(names(ldat) == "meanDate")
names(ldat)[ind] <- "date"

mdat <- malaWeal %>% inner_join(ldat, by = c("cluster" = "DHSCLUST")) %>%
  arrange(cluster)

ldat2 <- data.frame(cluster = unique(mdat$cluster)) %>%
  left_join(ldat, by = c("cluster" = "DHSCLUST"))

latlong <- data.frame(cluster = unique(mdat$cluster)) %>%
  left_join(latlong, by = c("cluster" = "DHSCLUST"))

shuff <- read.csv("Data/GH/shuff.csv")[,1]
h <- 40
u <- 0.01
v <- 0.1