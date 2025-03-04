#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config <- args[1]

if (!require("pacman")) install.packages("pacman", repos="http://cran.us.r-project.org")
pacman::p_load(glmnet, dplyr, data.table, rjson)

# library(glmnet)
# library(dplyr)
# library(data.table)
# library(rjson)

# FEATURE SHRINKING ATTEMPT #1: LOGISTIC LASSO
param <- fromJSON(file=config)
disease <- param$name
disease_col <- param$disease_column
transformation <- param$transformation

data <- read.csv(sprintf("data/%s/filtered_otu_table_%s.csv", disease, transformation), check.names=F)
X <- select(data, -c(1))
Y <- read.csv(sprintf("data/%s/metadata.csv", disease))[disease_col]
data[disease_col] <- Y

# balanced training set
set.seed(1)
num_healthy <- sum(Y == 0)
num_diseased <- sum(Y == 1)
min_train <- min(num_healthy, num_diseased)
healthy_df <- as.matrix(data[which(Y == 0),])
healthy_df_index <- sample(1:num_healthy, size=min_train)
diseased_df <- as.matrix(data[which(Y == 1),])
diseased_df_index <- sample(1:num_diseased, size=min_train)
data_train <- rbind(healthy_df[healthy_df_index, ], diseased_df[diseased_df_index, ])
X_train <- data.table(data_train[, -c(1, dim(data_train)[2])])
X_train <- apply(X_train, 2, function(x) as.numeric(x))
Y_train <- as.numeric(data_train[, disease_col])

cv.fit_train <- cv.glmnet(as.matrix(X_train), Y_train, family="binomial", type.measure="auc")
png(sprintf("plots/%s/LASSO_AUC_%s.png", disease, transformation))
plot(cv.fit_train)
dev.off()
best_lambda_train <- cv.fit_train$lambda.min

coef_train <- coef(cv.fit_train, s = "lambda.min")
coef_train <- coef_train[coef_train[,1] != 0,] 
fwrite(as.list(coef_train), file=sprintf("data/%s/lasso_covariates_%s.txt", disease, transformation))
