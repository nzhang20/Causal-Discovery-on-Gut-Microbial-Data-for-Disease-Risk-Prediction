library(CVglasso)
library(dplyr)
library(data.table)

# CVglasso example
set.seed(123)

# generate data from a sparse oracle precision matrix
# we will try to estimate this matrix from the data

# first compute the oracle covariance matrix
S = matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    S[i, j] = S[i, j]^abs(i - j)
  }
}

# print oracle precision matrix (shrinkage might be useful)
(Omega = round(qr.solve(S), 3))

# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt

# calculate sample covariance
sample = (nrow(X) - 1)/nrow(X)*cov(X)

# print sample precision matrix (perhaps a bad estimate)
round(qr.solve(cov(X)), 5)

# CVglasso (lam = 0.5)
CVglasso(S = sample, lam = 0.5)

# cross validation
(CV = CVglasso(X, trace = "none"))
