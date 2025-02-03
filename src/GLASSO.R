install.packages(c("glasso", "dplyr", "data.table", "rjson"), repos="http://cran.us.r-project.org")
library(glasso) 
library(dplyr)
library(data.table)
library(rjson)

set.seed(1)

param <- fromJSON(file="src/data-params.json")
disease <- param$disease
if (disease == 'pcos') {
  group0 <- 'hc'
  group1 <- 'pcos'
}

if(disease == 't2d') {
  group0 <- 'IS'
  group1 <- 'IR'
}

# healthy cohort
hc <- read.csv(sprintf("data/%s/%s.csv", disease, group0))
hc.X <- select(hc, -c(X))
hc.s <- var(hc.X)
hc.a <- glasso(hc.s, rho=0.01, nobs=dim(hc)[1])
hc.aa <- glasso(hc.s, rho=2, w.init=hc.a$w, wi.init=hc.a$wi)
fwrite(hc.aa$wi, file=sprintf("data/%s/glasso_%s.csv", disease, group0))

# diseased cohort
pcos <- read.csv(sprintf("data/%s/%s.csv", disease, group1))
pcos.X <- select(pcos, -c(X))
pcos.s <- var(pcos.X)
pcos.a <- glasso(pcos.s, rho=0.01, nobs=dim(pcos)[1])
pcos.aa <- glasso(pcos.s, rho=2, w.init=pcos.a$w, wi.init=pcos.a$wi)
fwrite(pcos.aa$wi, file=sprintf("data/%s/glasso_%s.csv", disease, group1))
