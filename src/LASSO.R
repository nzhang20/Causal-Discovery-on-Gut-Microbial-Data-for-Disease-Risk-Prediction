install.packages(c("glmnet", "dplyr", "data.table"), repos="http://cran.us.r-project.org")
library(glmnet)
library(dplyr)
library(data.table)

# FEATURE SHRINKING ATTEMPT #1: LOGISTIC LASSO
data <- read.csv("data/clean.csv")
X <- select(data, -c(X, group))
Y <- data$group

fit <- glmnet(X, Y)

# cross validation
# set.seed(1)
# cv.fit <- cv.glmnet(as.matrix(X), Y, family="binomial")
# plot(cv.fit)
# best_lambda <- cv.fit$lambda.min

# try with the recommended training set (400 random PCOS, 400 random PCOS)
set.seed(1)
data_train <- rbind(sample_n(data[which(data$group == 0),], 100), sample_n(data[which(data$group == 1),], 100))
X_train <- select(data_train, -c(X, group))
Y_train <- data_train$group

cv.fit_train <- cv.glmnet(as.matrix(X_train), Y_train, family="binomial", type.measure="auc")
png("plots/LASSO_AUC.png")
plot(cv.fit_train)
dev.off()
best_lambda_train <- cv.fit_train$lambda.min

coef_train <- coef(cv.fit_train, s = "lambda.min")
coef_train <- coef_train[coef_train[,1] != 0,] # shrank to 48 features
coef_train

# FEATURE SHRINKING ATTEMPT #2: SURE SCREENING
# sis <- sapply(X, function(x) cor(x, Y, method="spearman"))
# sis <- sort(abs(sis), decreasing=T)
# hist(sis, breaks=20)
# sis_theoretical_threshold <- sqrt(log(1128)/948)
# sis <- sis[sis > sis_theoretical_threshold]

##### RESULTS & COMMENTARY
# Since Logistic Regression LASSO yielded 48 features 
# and sure screening via Spearman correlations yielded (theoretically) 57 features, 
# we will go ahead and run our algorithms via the logistic regression LASSO pruned features.
fwrite(as.list(coef_train), file="src/LASSO_covariates.txt")
