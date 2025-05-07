SplitkFold <- function(Mat, Classes, K){
  
  require(dplyr)
  require(caret)
  require(glmnet)
  
  ## Downsample
  downsample <- downSample(Mat, Classes, list = TRUE)
  row.names(downsample[["x"]]) <- downsample[["x"]]$sample
  downsample_data <- downsample[["x"]]$sample <- NULL
  
  df <- data.frame(ID = row.names(downsample[["x"]]), Classes = downsample[["y"]])
  samples <- createFolds(df$Classes, k = K, list = TRUE, returnTrain=TRUE)
  
  return(list(df = df, samples = samples, data = downsample[["x"]]))
}

SplitFunction <- function(Mat, Classes) {
  
  require(dplyr)
  require(caret)
  require(glmnet)
  
  ##Split into training and test sets
  
  df <- data.frame(ID = rownames(Mat), Classes = Classes)
  samples <- createDataPartition(df$Classes, p = 0.8, times = 100)
  
  return(list(df = df, samples = samples))
}
