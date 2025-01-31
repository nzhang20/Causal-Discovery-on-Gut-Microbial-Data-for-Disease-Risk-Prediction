install.packages(c("glasso", "dplyr", "data.table"), repos="http://cran.us.r-project.org")
library(glasso) 
library(dplyr)
library(data.table)

set.seed(1)

# HC
hc <- read.csv("data/hc.csv")
hc.X <- select(hc, -c(X))
hc.s <- var(hc.X)
hc.a <- glasso(hc.s, rho=0.01, nobs=435)
hc.aa <- glasso(hc.s, rho=2, w.init=hc.a$w, wi.init=hc.a$wi)
fwrite(hc.aa$wi, file="data/glasso_hc.csv")

# PCOS
pcos <- read.csv("data/pcos.csv")
pcos.X <- select(pcos, -c(X))
pcos.s <- var(pcos.X)
pcos.a <- glasso(pcos.s, rho=0.01, nobs=513)
pcos.aa <- glasso(pcos.s, rho=2, w.init=pcos.a$w, wi.init=pcos.a$wi)
fwrite(pcos.aa$wi, file="data/glasso_pcos.csv")
