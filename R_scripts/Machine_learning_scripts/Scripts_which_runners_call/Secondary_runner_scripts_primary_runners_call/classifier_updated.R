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

############################# load the input data here ###############################
# Data = n_sample * n_feature dataframe # the features are the bins
# y = vector of length n_sample, each element represents the label of the sample (positive or negative)

############################# classifer starts here###############################

Features.CVparam <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = TRUE, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE)

gbmGrid <- expand.grid(interaction.depth = 3,
                       n.trees = 150,
                       shrinkage = 0.1,
                       n.minobsinnode = 3)

# Split the data into 10 folds
Splits10 <- list()
All.kFold <- list()
#CHARM.kfold <- list()

if(algorithm == "gbm") {
  for (j in 1:20) {
    # split the data for 20 times, result of each split is saved in Splits10
    Splits10[[j]] <- SplitkFold(Data, y, 10)
    
    #10-fold cross validation to calculate the probability (score) of each sample being cancer
    kFold.list <- list()
    kFold.charm <- list()
    
    for(i in 1:10) {
      Indices = Splits10[[j]]$samples[[i]]
      classes.df = Splits10[[j]]$df
      
      TrainData <- Splits10[[j]][["data"]][Indices, ]
      TrainPheno <- classes.df[Indices,]
      
      TestData <- Splits10[[j]][["data"]][!(row.names(Splits10[[j]][["data"]]) %in% row.names(TrainData)), ]
      TestPheno <- classes.df[classes.df$ID %in% row.names(TestData), ]
      
      #Train Model and predict test set
      Model <- train(x = TrainData, y = TrainPheno$Classes, trControl = Features.CVparam, method = algorithm, metric = "Kappa", tuneGrid = gbmGrid)
      Prediction.classProbs <- predict(Model, newdata = TestData, type = "prob") %>% data.frame
      
      Prediction.classProbs$ActualClass <- TestPheno$Classes
      Prediction.classProbs$PredictedClass <- predict(Model, newdata = TestData, type = "raw")
      Prediction.classProbs$sample <- row.names(TestData)
      
      kFold.list[[i]] <- Prediction.classProbs
      
      #Predict CHARM data
   #   Prediction.charm <- predict(Model, newdata = Data_charm, type = "prob") %>% data.frame
      
    #  Prediction.charm$PredictedClass <- predict(Model, newdata = Data_charm, type = "raw")
    #  Prediction.charm$sample <- row.names(Data_charm)
      
    #  kFold.charm[[i]] <- Prediction.charm
    }
    
    All.kFold[[j]] <- kFold.list
 #   CHARM.kfold[[j]] <- kFold.charm
  }
} else {
  for (j in 1:20) {
    # split the data for 20 times, result of each split is saved in Splits10
    Splits10[[j]] <- SplitkFold(Data, y, 10)
    
    #10-fold cross validation to calculate the probability (score) of each sample being cancer
    kFold.list <- list()
#    kFold.charm <- list()
    
    for(i in 1:10) {
      Indices = Splits10[[j]]$samples[[i]]
      classes.df = Splits10[[j]]$df
      
      TrainData <- Splits10[[j]][["data"]][Indices, ]
      TrainPheno <- classes.df[Indices,]
      
      TestData <- Splits10[[j]][["data"]][!(row.names(Splits10[[j]][["data"]]) %in% row.names(TrainData)), ]
      TestPheno <- classes.df[classes.df$ID %in% row.names(TestData), ]
      
      #Train Model and predict test set
      Model <- train(x = TrainData, y = TrainPheno$Classes, trControl = Features.CVparam, method = algorithm, metric = "Kappa")
      Prediction.classProbs <- predict(Model, newdata = TestData, type = "prob") %>% data.frame
      
      Prediction.classProbs$ActualClass <- TestPheno$Classes
      Prediction.classProbs$PredictedClass <- predict(Model, newdata = TestData, type = "raw")
      Prediction.classProbs$sample <- row.names(TestData)
      
      kFold.list[[i]] <- Prediction.classProbs
      
      #Predict CHARM data - this is if want to verify with new dataset. Can be a mix of M4 files and new healthy controls in the same
      # Format as previous dataframe with the columns and the fragment ratio
   #   Prediction.charm <- predict(Model, newdata = Data_charm, type = "prob") %>% data.frame
      
    #  Prediction.charm$PredictedClass <- predict(Model, newdata = Data_charm, type = "raw")
     # Prediction.charm$sample <- row.names(Data_charm)
      
      #kFold.charm[[i]] <- Prediction.charm
    }
    
    All.kFold[[j]] <- kFold.list
    print(paste0("Iteration number", j, "has completed"))
 #   CHARM.kfold[[j]] <- kFold.charm
  }
}
