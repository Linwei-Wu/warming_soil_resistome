
library(randomForest)
library(dplyr)

#### perform random forest based on envs 
# dat <- read.csv("dataset for random forest.csv", row.names = 1)
## 1. construct the model using training dataset 
# year 2009 to 2018
dat_train <- dat[dat$year < 2019, ]
dat_test <- dat[dat$year >= 2019, ]

models <- dat_train %>% 
  group_by(year) %>% 
  do(model = randomForest(ARG_cov_gb ~ ., data = ., importance=T))


saveRDS(models, file="randomForest_models_09_18_envs.RData")

## 2. model prediction on test dataset 

# model prediction on each year-based model
predictions <- lapply(models$model, predict, newdata = dat_test)

# get the average of models
pre_df <- do.call(cbind, predictions)
avg_predictions <- rowMeans(pre_df)

mnames <- paste0("model_", 2009:2018)
colnames(pre_df) <- mnames

sum(row.names(dat_test) == row.names(pre_df))
sum(row.names(dat_test) == names(avg_predictions))

res <- data.frame(ARG_covGB_observed = dat_test$ARG_cov_gb,
                  avg_predictions, pre_df)


plot(res$avg_predictions, res$ARG_covGB_observed)
summary(lm(res$avg_predictions ~ res$ARG_covGB_observed))
write.csv(res, "model predictions_envs.csv")

## 3. model prediction on training dataset 
models <- readRDS("randomForest_models_09_18_envs.RData")

# model prediction on each year-based model
predictions <- lapply(models$model, predict, newdata = dat_train)

# get the average of models
pre_df <- do.call(cbind, predictions)
avg_predictions <- rowMeans(pre_df)

mnames <- paste0("model_", 2009:2018)
colnames(pre_df) <- mnames

sum(row.names(dat_train) == row.names(pre_df))
sum(row.names(dat_train) == names(avg_predictions))

res <- data.frame(ARG_covGB_observed = dat_train$ARG_cov_gb,
                  avg_predictions, pre_df)


plot(res$avg_predictions, res$ARG_covGB_observed)
summary(lm(res$avg_predictions ~ res$ARG_covGB_observed))

write.csv(res, "model predictions_envs_training dataset.csv")

##4. get the importance of predictors 

importances <- sapply(models$model, function(x) {
  imp <- x$importance 
  incmse <- imp[, 1]
  incmse
})
colnames(importances) <- mnames
write.csv(importances, "predictor_importance_envs.csv")

#sort(rowMeans(importances, na.rm = T), decreasing = T)




