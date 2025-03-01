# install.packages(c("glasso", "dplyr", "data.table", "rjson"), repos="http://cran.us.r-project.org")
library(glasso) 
library(dplyr)
library(data.table)
library(rjson)

set.seed(1)

param <- fromJSON(file="src/config-t2d.json")
disease <- param$name
groups <- param$cohort_names
group0 <- groups[1]
group1 <- groups[2]
transformation <- param$transformation

# healthy cohort
hc <- read.csv(sprintf("data/%s/%s_%s.csv", disease, group0, transformation))
hc.X <- select(hc, -c(SampleID))
hc.s <- var(hc.X)
hc.a <- glasso(hc.s, rho=0.01, nobs=dim(hc)[1])
hc.aa <- glasso(hc.s, rho=0.02, w.init=hc.a$w, wi.init=hc.a$wi)
fwrite(hc.aa$wi, file=sprintf("data/%s/glasso_%s_%s.csv", disease, group0, transformation))

# diseased cohort
pcos <- read.csv(sprintf("data/%s/%s_%s.csv", disease, group1, transformation))
pcos.X <- select(pcos, -c(SampleID))
pcos.s <- var(pcos.X)
pcos.a <- glasso(pcos.s, rho=0.01, nobs=dim(pcos)[1])
pcos.aa <- glasso(pcos.s, rho=0.02, w.init=pcos.a$w, wi.init=pcos.a$wi)
fwrite(pcos.aa$wi, file=sprintf("data/%s/glasso_%s_%s.csv", disease, group1, transformation))
