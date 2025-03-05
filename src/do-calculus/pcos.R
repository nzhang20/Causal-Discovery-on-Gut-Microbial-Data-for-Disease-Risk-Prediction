library(tidyverse)
library(dplyr)

otu_table <- read.csv('~/academics/Causal-Discovery-for-Biomedical-Applications/data/pcos/filtered_otu_table_norm.csv')
metadata <- read.csv('~/academics/Causal-Discovery-for-Biomedical-Applications/data/pcos/metadata.csv')
data <- merge(metadata, otu_table, by='X')
data <- select(data, -c(X))
data$group <- as.factor(data$group)
data$region <- as.factor(data$region)
data$study_site <- as.factor(data$study_site)

##### MODEL 1 #####
model.1.data <- select(data, group, Alistipes, Blautia,
                       Burkholderia, Desulfovibrio, 
                       Holdemanella, Knoellia,
                       Prevotellaceae_NK3B31_group,
                       Ruminococcus, Ruminococcus_gnavus_group)
model1 <- glm(group ~ ., family=binomial, data=model.1.data)
summary(model1)

##### MODEL 2 #####

# Alistipes
alis.data <- select(data, group, Alistipes,
                    Bacteroides, UCG.005, UCG.002, region)
alis <- glm(group ~ ., family=binomial, data=alis.data)
summary(alis)

# Blautia
blau.data <- select(data, group, Blautia,
                    region, Subdoligranulum)
blau <- glm(group ~ ., family=binomial, data=blau.data)
summary(blau)

# Burkholderia
burk.data <- select(data, group, Burkholderia,
                    Bacteroides, Schaedlerella, Kineothrix)
burk <- glm(group ~ ., family=binomial, data=burk.data)
summary(burk)

# Desulfovibrio
desulf.data <- select(data, group, Desulfovibrio,
                      region, Clostridia_vadinBB60_group)
desulf <- glm(group ~ ., family=binomial, data=desulf.data)
summary(desulf)

# Holdemanella
hold.data <- select(data, group, Holdemanella, 
                    region, Christensenellaceae_R.7_group, Subdoligranulum,
                    UCG.010, Bacteroides)
hold <- glm(group ~ ., family=binomial, data=hold.data)
summary(hold)

# Knoellia
knoe.data <- select(data, group, Knoellia, Duganella, Acinetobacter)
knoe <- glm(group ~ ., family=binomial, data=knoe.data)
summary(knoe)

# Prevotellaceae_NK3B31_group
prev.nk3b31.data <- select(data, group, Prevotellaceae_NK3B31_group)
prev.nk3b31 <- glm(group ~ ., family=binomial, data=prev.nk3b31.data)
summary(prev.nk3b31)

# Ruminococcus
rumi.data <- select(data, group, Ruminococcus, Subdoligranulum,
                    Clostridia_UCG.014, Eubacterium)
rumi <- glm(group ~ ., family=binomial, data=rumi.data)
summary(rumi)

# Ruminococcus_gnavus_group
rumi.gnav.data <- select(data, group, Ruminococcus_gnavus_group, 
                         Bacteroides, study_site, Lachnoclostridium)
rumi.gnav <- glm(group ~ ., family=binomial, data=rumi.gnav.data)
summary(rumi.gnav)
