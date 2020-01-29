install.packages("glmnet_2.0-16")
BiocManager::install("glmnet")
require(devtools)
install_version("glmnet", version = "2.0-16", repos = "http://cran.us.r-project.org")

# Loaging the library
library(glmnet)

# Loading the data
data(swiss)

x_vars <- model.matrix(Fertility~. , swiss)[,-1]
y_var <- swiss$Fertility
data <- cbind(x_vars,y_var)


# Splitting the data into test and train
set.seed(100) 

index = sample(1:nrow(data), 0.7*nrow(data)) 

train = data[index,] # Create the training data 
test = data[-index,] # Create the test data

dim(train)
dim(test)
#### Scaling the Numeric Features
library(caret)
cols = c("Agriculture", "Examination", "Education", "Catholic", "Infant.Mortality")

pre_proc_val <- preProcess(train[,cols], method = c("center", "scale"))

train[,cols] = predict(pre_proc_val, train[,cols])
test[,cols] = predict(pre_proc_val, test[,cols])
train <- as.data.frame(train)
summary(train)
## Linear Regression
lr = lm(y_var ~ Agriculture + Examination + Education + Catholic + Infant.Mortality, data = train)
summary(lr)

####### Regularization ############
cols_reg = c("Agriculture", "Examination", "Education", "Catholic", "Infant.Mortality", "y_var")

dummies <- dummyVars(y_var ~ ., data = data[,cols_reg])

train_dummies = predict(dummies, newdata = train[,cols_reg])

test_dummies = predict(dummies, newdata = test[,cols_reg])

print(dim(train_dummies)); print(dim(test_dummies))

### Ridge Regression
library("glmnet")

x = as.matrix(train_dummies)
y_train = train$y_var

x_test = as.matrix(test_dummies)
y_test = test$y_var

lambdas <- 10^seq(2, -3, by = -.1)
ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)

summary(ridge_reg)


cv_ridge <- cv.glmnet(x, y_train, alpha = 0, lambda = lambdas)
optimal_lambda <- cv_ridge$lambda.min
optimal_lambda



### Lasso Regression
lambdas <- 10^seq(2, -3, by = -.1)

# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(x, y_train, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)

# Best 
lambda_best <- lasso_reg$lambda.min 
lambda_best


lasso_model <- glmnet(x, y_train, alpha = 1, lambda = lambda_best, standardize = TRUE)
coef(lasso_model)
predictions_train <- predict(lasso_model, s = lambda_best, newx = x)
eval_results(y_train, predictions_train, train)

predictions_test <- predict(lasso_model, s = lambda_best, newx = x_test)
eval_results(y_test, predictions_test, test)