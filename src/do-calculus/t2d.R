library(tidyverse)
library(dplyr)

otu_table <- read.csv('~/academics/Causal-Discovery-for-Biomedical-Applications/data/t2d/filtered_otu_table_norm.csv')
metadata <- read.csv('~/academics/Causal-Discovery-for-Biomedical-Applications/data/t2d/metadata.csv')
data <- merge(metadata, otu_table, by='SampleID')
data <- select(data, -c(SampleID))
data$IRIS <- as.factor(data$IRIS)
data.otu <- select(data, -c(IRIS, Gender, Ethnicity))

##### MODEL 1 #####
model.1.data <- select(data, IRIS, genus_Butyricimonas, 
                       genus_Clostridium.XlVb,
                       genus_Odoribacter,
                       genus_unclassified_Bacteria,
                       genus_unclassified_Firmicutes)
model1 <- glm(IRIS ~ ., family=binomial, data=model.1.data)
summary(model1)

##### MODEL 2 #####

# Butyricimonas
buty.data <- select(data, IRIS, genus_Butyricimonas, 
                    genus_Clostridium.XlVa,
                    genus_Oscillibacter)
buty <- glm(IRIS ~ ., family=binomial, data=buty.data)
summary(buty)

# Clostridium XIVb
clos.xivb.data <- select(data, IRIS, genus_Clostridium.XlVb,
                         genus_Flavonifractor, 
                         genus_Parasutterella,
                         genus_Blautia)
clos.xivb <- glm(IRIS ~ ., family=binomial, data=clos.xivb.data)
summary(clos.xivb)

# Odoribacter (mediation)
odor.data <- select(data, IRIS, genus_Odoribacter,
                    genus_unclassified_Bacteria,
                    genus_Alistipes,
                    genus_Barnesiella)
u.bact.med <- lm(genus_unclassified_Bacteria ~ genus_Odoribacter, data=data)
summary(u.bact.med)
odor.u.bact <- coefficients(u.bact.med)[2]
odor <- glm(IRIS ~ ., family=binomial, data=odor.data)
summary(odor)
odor.effect <- coefficients(odor)[2] + (odor.u.bact * coefficients(odor)[3])

# unclassified Bacteria
u.bact.data <- select(data, IRIS, genus_unclassified_Bacteria,
                      genus_Odoribacter,
                      genus_unclassified_Firmicutes,
                      genus_Oscillibacter,
                      genus_unclassified_Lachnospiraceae,
                      Ethnicity)
u.bact <- glm(IRIS ~ ., family=binomial, data=u.bact.data)
summary(u.bact)

# unclassified Firmicutes (mediation)
u.firm.data <- select(data, IRIS, genus_unclassified_Firmicutes,
                      genus_unclassified_Bacteria,
                      genus_Clostridium.sensu.stricto)
u.bact.med <- lm(genus_unclassified_Bacteria ~ genus_unclassified_Firmicutes, data=data)
summary(u.bact.med)
u.firm.u.bact <- coefficients(u.bact.med)[2]
u.firm <- glm(IRIS ~ ., family=binomial, data=u.firm.data)
summary(u.firm)
u.firm.effect <- coefficients(u.firm)[2] + (u.firm.u.bact * coefficients(u.firm)[3])

