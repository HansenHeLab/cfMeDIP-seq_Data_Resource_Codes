library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(NMF)
library(Rtsne)
library(caret)
library(randomForest)
library(glmnet)
library(limma)

############################# Load the input data ###############################
# Data: a dataframe or matrix with dimensions n_sample x n_feature
# y: a vector of length n_sample with class labels (e.g., "positive", "negative")

############################# Classifier Optimization ###############################

# Check algorithm defined
if (!exists("algorithm") || is.null(algorithm)) {
  stop("ERROR: 'algorithm' is not defined or is NULL. Please specify an algorithm.")
}

if (!is.character(algorithm)) {
  algorithm <- as.character(algorithm)
}

cat("DEBUG: Selected algorithm is:", algorithm, "\n")


# Define a caret trainControl for repeated 10-fold cross-validation
Features.CVparam <- trainControl(method = "repeatedcv", 
                                 number = 10, 
                                 repeats = 1, 
                                 verboseIter = TRUE, 
                                 returnData = FALSE, 
                                 classProbs = TRUE, 
                                 savePredictions = "final")

# Define tuning grids for models that require hyperparameter tuning
gbmGrid <- expand.grid(interaction.depth = c(2, 3, 4),
                       n.trees = c(100, 150, 200),
                       shrinkage = c(0.01, 0.1),
                       n.minobsinnode = c(3, 5))

xgbGrid <- expand.grid(nrounds = c(100, 200),
                       max_depth = c(3, 6),
                       eta = c(0.01, 0.1),
                       gamma = c(0, 1),
                       colsample_bytree = c(0.5, 0.8),
                       min_child_weight = c(1, 3),
                       subsample = c(0.7, 1))

# In caret, Lasso is implemented via glmnet with alpha=1
lassoGrid <- expand.grid(alpha = 1,
                         lambda = seq(0.0001, 0.1, length = 10))

plsGrid <- expand.grid(ncomp = 1:5)  # Number of components for PLS

# Split the data into 10 folds.
# (Note: downSample() returns a list containing balanced training data in $x and $y.)
Splits <- list()
kappa <- c()

for (i in 1:10) {
  # Here, we use downSample() to generate balanced folds.
  Splits[[i]] <- downSample(Data, y, list = TRUE)
}

# Loop through each fold and train the model with the chosen algorithm.
for (j in 1:length(Splits)) {
  
  TrainData <- Splits[[j]][["x"]]
  TrainPheno <- Splits[[j]][["y"]]
  
  if (algorithm == "gbm") {
    Model <- train(x = TrainData, y = TrainPheno, 
                   trControl = Features.CVparam, 
                   method = "gbm", 
                   metric = "Kappa", 
                   tuneGrid = gbmGrid)
  
  } else if (algorithm == "xgbTree") {
    Model <- train(x = TrainData, y = TrainPheno, 
                   trControl = Features.CVparam, 
                   method = "xgbTree", 
                   metric = "Kappa", 
                   tuneGrid = xgbGrid)

  } else if (algorithm == "glmnet") {  # for Lasso use glmnet with alpha = 1
    Model <- train(x = TrainData, y = TrainPheno, 
                   trControl = Features.CVparam, 
                   method = "glmnet", 
                   metric = "Kappa", 
                   tuneGrid = lassoGrid)

  } else if (algorithm == "pls") {
    Model <- train(x = TrainData, y = TrainPheno, 
                   trControl = Features.CVparam, 
                   method = "pls", 
                   metric = "Kappa", 
                   tuneGrid = plsGrid)

  } else if (algorithm %in% c("rf", "glm", "svmRadial", "knn", "bayesglm", "lda", "qda")) {
    # For these models, we use caret’s default tuning grid.
    Model <- train(x = TrainData, y = TrainPheno, 
                   trControl = Features.CVparam, 
                   method = algorithm, 
                   metric = "Kappa")
  } else {
    stop(paste("ERROR: Model", algorithm, "is not supported."))
  }
  
  # Append the Kappa performance metric from the current fold
  kappa <- c(kappa, Model[["results"]][["Kappa"]])
}

