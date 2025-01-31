install.packages(c("glmnet", "dplyr", "data.table"), repos="http://cran.us.r-project.org")
library(glmnet)
library(dplyr)
library(data.table)

# FEATURE SHRINKING ATTEMPT #1: LOGISTIC LASSO
data <- read.csv("data/filter_rare.csv")
X <- select(data, -c(X, group))
Y <- data$group

# balanced training set (435 HC, 435 random PCOS)
set.seed(1)
data_train <- rbind(sample_n(data[which(data$group == 0),], 435), sample_n(data[which(data$group == 1),], 435))
X_train <- select(data_train, -c(X, group))
Y_train <- data_train$group

cv.fit_train <- cv.glmnet(as.matrix(X_train), Y_train, family="binomial", type.measure="auc")
png("plots/LASSO_AUC.png")
plot(cv.fit_train)
dev.off()
best_lambda_train <- cv.fit_train$lambda.min

coef_train <- coef(cv.fit_train, s = "lambda.min")
coef_train <- coef_train[coef_train[,1] != 0,] 
fwrite(as.list(coef_train), file="data/lasso_covariates.txt")
