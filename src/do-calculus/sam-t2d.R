library(tidyverse)
library(dplyr)

otu_table <- read.csv('~/academics/Causal-Discovery-for-Biomedical-Applications/data/sam-t2d/filtered_otu_table_norm.csv')
metadata <- read.csv('~/academics/Causal-Discovery-for-Biomedical-Applications/data/sam-t2d/metadata.csv')
colnames(otu_table)[1] <- 'sample.id'
data <- merge(metadata, otu_table, by='sample.id')
data <- select(data, -c(sample.id))
data$t2d <- as.factor(data$t2d)
data.otu <- select(data, -c(t2d, study, region_num))

##### MODEL 1 #####
model.1.data <- data[c(
                       't2d',
                       'd__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Prevotella',
                       'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Ruminococcaceae.g__Faecalibacterium',
                       'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Coriobacteriaceae.g__Collinsella',
                       'd__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Rikenellaceae.g__Alistipes_A_871400',
                       'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Christensenellales.f__CAG.138.g__Aphodomorpha',
                       'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Lachnospirales.f__Lachnospiraceae.g__Pseudobutyrivibrio',
                       'd__Bacteria.p__Pseudomonadota.c__Alphaproteobacteria.o__RF32.f__CAG.239.g__CAG.267',
                       'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Monoglobales.f__Monoglobaceae.g__')]
model1 <- glm(t2d ~ ., family=binomial, data=model.1.data)
summary(model1)

##### MODEL 2 #####

# Alistipes_A_871400 (6)
alis.data <- data[c(
  't2d',
  'study',
  'd__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Rikenellaceae.g__Alistipes_A_871400',
  'd__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Marinifilaceae.g__Butyricimonas',
  'd__Bacteria.p__Bacillota_I.c__Bacilli_A.o__Erysipelotrichales.f__Coprobacillaceae.g__',
  'd__Bacteria.p__Bacillota_I.c__Bacilli_A.o__Erysipelotrichales.f__Erysipelotrichaceae.g__Merdibacter',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Ruminococcaceae.g__Anaerotruncus',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__UBA644.g__UBA644'
)]
alis <- glm(t2d ~ ., family=binomial, data=alis.data)
summary(alis)

# Aphodomorpha (6)
apho.data <- data[c(
  't2d',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Christensenellales.f__CAG.138.g__Aphodomorpha',
  'd__Bacteria.p__Bacillota_I.c__Bacilli_A.o__Erysipelotrichales.f__Erysipelotrichaceae.g__Faecalicoccus',
  'd__Bacteria.p__Bacillota_I.c__Bacilli_A.o__Erysipelotrichales.f__Coprobacillaceae.g__',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Christensenellales.f__Borkfalkiaceae.g__Borkfalkia',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.__.__.__',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__UBA644.g__UBA644',
  'd__Bacteria.p__Bacillota_I.c__Bacilli_A.o__Erysipelotrichales.f__Erysipelotrichaceae.g__Longicatena'
)]
apho <- glm(t2d ~ ., family=binomial, data=apho.data)
summary(apho)

# CAG267
cag267.data <- data[c(
  't2d',
  'd__Bacteria.p__Pseudomonadota.c__Alphaproteobacteria.o__RF32.f__CAG.239.g__CAG.267',
  'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Eggerthellaceae.__',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Lachnospirales.f__Lachnospiraceae.g__Pseudobutyrivibrio',
  'd__Bacteria.p__Pseudomonadota.c__Alphaproteobacteria.o__RF32.f__CAG.239.g__51.20'
)]
cag267 <- glm(t2d ~ ., family=binomial, data=cag267.data)
summary(cag267)

# Collinsella (7)
coll.data <- data[c(
  't2d',
  'region_num',
  'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Coriobacteriaceae.g__Collinsella',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Ruminococcaceae.g__Faecalibacterium',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Lachnospirales.f__Lachnospiraceae.g__Anaerobutyricum',
  'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Eggerthellaceae.__',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Lachnospirales.f__Lachnospiraceae.g__Pseudobutyrivibrio',
  'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Atopobiaceae.g__Tractidigestivibacter',
  'd__Bacteria.p__Synergistota.c__Synergistia.o__Synergistales.f__Dethiosulfovibrionaceae.g__Pyramidobacter'
)]
coll <- glm(t2d ~ ., family=binomial, data=coll.data)
summary(coll)

# Faecalibacterium (5)
faec.data <- data[c(
  't2d',
  'region_num',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Ruminococcaceae.g__Faecalibacterium',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Butyricicoccaceae.g__Agathobaculum',
  'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Coriobacteriaceae.g__Collinsella',
  'd__Bacteria.p__Fusobacteriota.c__Fusobacteriia.o__Fusobacteriales.f__Fusobacteriaceae.g__Fusobacterium_A',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__UBA1381.f__UBA1381.g__Hominilimicola'
)]
faec <- glm(t2d ~ ., family=binomial, data=faec.data)
summary(faec)

# Prevotella (5)
prev.data <- data[c(
  't2d',
  'd__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Prevotella', 
  'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Eggerthellaceae.g__Eggerthella',
  'd__Bacteria.p__Bacillota_C.c__Negativicutes.o__Selenomonadales.f__Selenomonadaceae_42771.g__Mitsuokella',
  'd__Bacteria.p__Pseudomonadota.c__Gammaproteobacteria.o__Enterobacterales_737866.f__Succinivibrionaceae.g__Succinivibrio',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Acutalibacteraceae.g__Ruminococcus_E',
  'd__Bacteria.p__Pseudomonadota.c__Alphaproteobacteria.o__RF32.f__CAG.239.g__RUG410'
)]
prev <- glm(t2d ~ ., family=binomial, data=prev.data)
summary(prev)

# Pseudobutyrivibrio (4)
pseu.data <- data[c(
  't2d',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Lachnospirales.f__Lachnospiraceae.g__Pseudobutyrivibrio',
  'd__Bacteria.p__Actinomycetota.c__Coriobacteriia.o__Coriobacteriales.f__Coriobacteriaceae.g__Collinsella',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Oscillospiraceae_88309.g__Lawsonibacter',
  'd__Bacteria.p__Bacillota_I.c__Bacilli_A.o__Erysipelotrichales.f__Erysipelotrichaceae.g__Longicatena',
  'd__Bacteria.p__Pseudomonadota.c__Alphaproteobacteria.o__RF32.f__CAG.239.g__CAG.267'
)]
pseu <- glm(t2d ~ ., family=binomial, data=pseu.data)
summary(pseu)

# unclassified Monoglobaceae
u.mono.data <- data[c(
  't2d',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Monoglobales.f__Monoglobaceae.g__',
  'd__Bacteria.p__Pseudomonadota.c__Gammaproteobacteria.o__Enterobacterales_737866.f__Enterobacteriaceae_A_725029.g__Escherichia',
  'd__Bacteria.p__Bacillota_A_368345.c__Clostridia_258483.o__Oscillospirales.f__Acutalibacteraceae.g__CAG.488',
  'd__Bacteria.p__Bacteroidota.c__Bacteroidia.o__Bacteroidales.f__Marinifilaceae.g__Butyricimonas'
)]
u.mono <- glm(t2d ~ ., family=binomial, data=u.mono.data)
summary(u.mono)
