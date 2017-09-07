# Initialize for Cross Validations
source("cvFunctions.R")

malaWeal <- read.csv("Data/TG/MalaEduWeal.csv") %>% 
  dplyr::select(cluster, hh, micro_mala, wealth, Edu = maxEdu, age) %>%
  filter(age < 72) %>%
  arrange(cluster)

malaWeal$age <- classifyAge(malaWeal$age)

wdat <- malaWeal
ldat <- read.csv("Data/TG/TGDHS2013Env.csv") %>%
  left_join(read.csv("Data/TG/TGTime.csv"), by = c("DHSCLUST" = "cluster"))

ldat$meanDate <- scale(ldat$meanDate) %>% as.vector()
names(ldat)[17] <- "date"

ldat$GRUMP[ldat$GRUMP == 2] <- 0
ldat$GRUMP2 <- ldat$GRUMP
medUrbPop <- aggregate(ldat$POPDEN ~ ldat$GRUMP, FUN = median)[1,2]
ldat$GRUMP2[ldat$POPDEN < medUrbPop & ldat$GRUMP == 0] <- 0.5
ldat$LCType <- as.factor(ldat$LCType)
ldat$GRUMP2 <- as.factor(ldat$GRUMP2)

latlong <- ldat[,1:3]
ldat[,2:3] <- scale(ldat[,2:3])

mdat <- malaWeal %>% inner_join(ldat, by = c("cluster" = "DHSCLUST")) %>%
  arrange(cluster)

ldat2 <- data.frame(cluster = unique(mdat$cluster)) %>%
  left_join(ldat, by = c("cluster" = "DHSCLUST"))

latlong <- data.frame(cluster = unique(mdat$cluster)) %>%
  left_join(latlong, by = c("cluster" = "DHSCLUST"))

shuff <- read.csv("Data/TG/shuff.csv")[,1]

h <- 32
u <- 1
v <- 0.01
