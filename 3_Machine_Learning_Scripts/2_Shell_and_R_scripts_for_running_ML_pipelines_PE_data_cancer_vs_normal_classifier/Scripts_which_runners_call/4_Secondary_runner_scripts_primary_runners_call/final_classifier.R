################################################################################
# final_classifier.R
#
# This script implements a robust classification approach that:
#   1) Automatically decides between LOOCV (if N < 10) and Nested CV (otherwise).
#   2) Dynamically picks the outer and inner folds based on sample size:
#       - 10 <= N <= 100: Nested CV  (5×5)
#       - 100 <  N <= 250: Nested CV (10×5)
#       - N > 250: Nested CV (10×10)
#   3) Includes hyperparameter tuning in the inner loop for algorithms like GBM.
#   4) Performs external validation predictions if a `data.external.validation`
#      object is provided (and not NULL).
#   5) Outputs key objects in the global environment for easy saving:
#       - All.kFold        (internal CV predictions)
#       - External.kfold   (external predictions per fold, if any)
#       - final_model      (model trained on the entire dataset)
#       - external_results (final model predictions on the external data)
#       - nested_cv_results (detailed summary of folds, confusion matrices, etc.)
#
# Requirements in your environment:
#   - Data: a matrix or data frame (row = samples, col = features)
#   - y: a factor vector of length nrow(Data) with class labels
#   - algorithm: a string, e.g. "gbm", "rf", "glm", "svmRadial", "knn", etc.
#   - data.external.validation (optional): external dataset with same columns as Data
#
# Usage:
#   source("final_classifier.R")
#   # Then, All.kFold, External.kfold, final_model, external_results, nested_cv_results
#   # will be created in the global environment.
#
# Debugging:
#   - Look for messages of the form: "DEBUG:" in the script output to see progress.
################################################################################

# (A) Load libraries
cat("DEBUG: Loading required packages...\n")
library(caret)
library(dplyr)
library(doParallel)
library(purrr)

# (B) Global debugging message
cat("DEBUG: Checking if 'Data', 'y', 'algorithm' are in the environment...\n")
if (!exists("Data")) stop("ERROR: 'Data' is not defined in the environment.")
if (!exists("y")) stop("ERROR: 'y' is not defined in the environment.")
if (!exists("algorithm")) stop("ERROR: 'algorithm' is not defined in the environment.")

# (C) Determine sample size
cat("DEBUG: Determining sample size of 'Data'...\n")
n_samples <- nrow(Data)
cat(paste("DEBUG: n_samples =", n_samples, "\n"))

# (D) Decide LOOCV or Nested CV, and fold sizes
cat("DEBUG: Setting cross-validation strategy based on sample size...\n")
use_loocv <- FALSE
outer_folds <- NA
inner_folds <- NA

if (n_samples < 10) {
  use_loocv <- TRUE
  cat("DEBUG: N < 10; Will use LOOCV.\n")
} else {
  # Nested CV
  if (n_samples <= 100) {
    outer_folds <- 5
    inner_folds <- 5
    cat("DEBUG: 10 <= N <= 100; Will use nested CV (5×5).\n")
  } else if (n_samples <= 250) {
    outer_folds <- 10
    inner_folds <- 5
    cat("DEBUG: 100 < N <= 250; Will use nested CV (10×5).\n")
  } else {
    outer_folds <- 10
    inner_folds <- 10
    cat("DEBUG: N > 250; Will use nested CV (10×10).\n")
  }
}

# (E) Register parallel backend
cat("DEBUG: Setting up parallel backend...\n")
num_cores <- max(parallel::detectCores() - 1, 1)
cat(paste("DEBUG: Number of cores detected =", parallel::detectCores(), "=> using", num_cores, "cores\n"))
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# (F) Define an inner trainControl for hyperparameter tuning
cat("DEBUG: Defining trainControl for hyperparameter tuning (inner loop)...\n")
innerControl <- trainControl(
  method         = "repeatedcv",  # repeated cross-validation for tuning
  number         = ifelse(use_loocv, 2, inner_folds), 
  # (If LOOCV, inner_folds is not used, so we default to 2 for caretaker reason.)
  repeats        = 3,
  classProbs     = TRUE,
  search         = "grid",
  summaryFunction= multiClassSummary, # or twoClassSummary for binary
  verboseIter    = FALSE
)

# (G) Example grid for GBM tuning; adapt for other algos if desired
cat("DEBUG: Defining a sample hyperparameter grid for different algorithms...\n")

# GBM Grid
cat("DEBUG: Defining GBM grid...\n")
gbmGrid <- expand.grid(
  interaction.depth = c(2, 3, 4),
  n.trees           = c(100, 150, 200),
  shrinkage         = c(0.01, 0.1),
  n.minobsinnode    = c(3, 5)
)
cat(paste("DEBUG: GBM grid defined with", nrow(gbmGrid), "combinations.\n"))

# XGBoost Grid
cat("DEBUG: Defining XGBoost grid...\n")
xgbGrid <- expand.grid(
  nrounds           = c(100, 200),
  max_depth         = c(3, 6),
  eta               = c(0.01, 0.1),
  gamma             = c(0, 1),
  colsample_bytree  = c(0.5, 0.8),
  min_child_weight  = c(1, 3),
  subsample         = c(0.7, 1)
)
cat(paste("DEBUG: XGBoost grid defined with", nrow(xgbGrid), "combinations.\n"))

# Lasso (GLMNet) Grid
cat("DEBUG: Defining Lasso (glmnet) grid...\n")
lassoGrid <- expand.grid(
  alpha  = 1,  # Lasso (Elastic Net would use alpha between 0 and 1)
  lambda = seq(0.0001, 0.1, length = 10)
)
cat(paste("DEBUG: Lasso grid defined with", nrow(lassoGrid), "combinations.\n"))

# PLS Grid
cat("DEBUG: Defining PLS grid...\n")
plsGrid <- expand.grid(
  ncomp = 1:5  # Number of components for PLS
)
cat(paste("DEBUG: PLS grid defined with", nrow(plsGrid), "combinations.\n"))

cat("DEBUG: Hyperparameter grid definitions completed.\n")


# (H) Helper function for training with the inner CV
cat("DEBUG: Creating helper function 'train_with_inner_cv'...\n")
train_with_inner_cv <- function(trainX, trainY, method, inner_ctrl) {
  cat("DEBUG: Entering 'train_with_inner_cv' function...\n")
  cat(paste("DEBUG: Training method:", method, "\n"))
  cat(paste("DEBUG: Training set dimensions: (", nrow(trainX), " samples,", ncol(trainX), " features)\n"))
  
  modelFit <- NULL  # Initialize modelFit
  
  if (method == "gbm") {
    cat("DEBUG: Training GBM with custom grid...\n")
    cat(paste("DEBUG: GBM Grid size:", nrow(gbmGrid), "parameter combinations\n"))
    modelFit <- train(
      x = trainX,
      y = trainY,
      method = "gbm",
      metric = "Accuracy",   # or "Kappa"
      trControl = inner_ctrl,
      tuneGrid = gbmGrid
    )
    
  } else if (method == "xgbTree") {
    cat("DEBUG: Training XGBoost with custom grid...\n")
    cat(paste("DEBUG: XGBoost Grid size:", nrow(xgbGrid), "parameter combinations\n"))
    modelFit <- train(
      x = trainX,
      y = trainY,
      method = "xgbTree",
      metric = "Accuracy",   # or "Kappa"
      trControl = inner_ctrl,
      tuneGrid = xgbGrid
    )
    
  } else if (method == "glmnet") {  # for Lasso, use glmnet with alpha=1
    cat("DEBUG: Training glmnet (Lasso) with custom grid...\n")
    cat(paste("DEBUG: Lasso Grid size:", nrow(lassoGrid), "parameter combinations\n"))
    modelFit <- train(
      x = trainX,
      y = trainY,
      method = "glmnet",
      metric = "Accuracy",   # or "Kappa"
      trControl = inner_ctrl,
      tuneGrid = lassoGrid
    )
    
  } else if (method == "pls") {
    cat("DEBUG: Training PLS with custom grid...\n")
    cat(paste("DEBUG: PLS Grid size:", nrow(plsGrid), "parameter combinations\n"))
    modelFit <- train(
      x = trainX,
      y = trainY,
      method = "pls",
      metric = "Accuracy",   # or "Kappa"
      trControl = inner_ctrl,
      tuneGrid = plsGrid
    )
    
  } else if (method %in% c("rf", "glm", "svmRadial", "knn", "bayesglm", "lda", "qda")) {
    cat(paste("DEBUG: Training", method, "using caret defaults (no custom grid)...\n"))
    modelFit <- train(
      x = trainX,
      y = trainY,
      method = method,
      metric = "Accuracy",   # or "Kappa"
      trControl = inner_ctrl
    )
    
  } else {
    stop(paste("ERROR: Model", method, "is not supported.\n"))
  }
  
  cat("DEBUG: Model training completed for", method, "\n")
  
  # Print best parameters if applicable
  if (!is.null(modelFit) && !is.null(modelFit$bestTune)) {
    cat("DEBUG: Best hyperparameters selected:\n")
    print(modelFit$bestTune)
  }
  
  return(modelFit)
}


# (I) Prepare placeholders for final outputs
cat("DEBUG: Initializing output data structures...\n")
All.kFold        <- list()  # per-fold internal CV predictions
External.kfold   <- list()  # per-fold external predictions
nested_cv_results <- list() # consolidated info

# (J) If we are using LOOCV
if (use_loocv) {
  cat("DEBUG: Starting LOOCV procedure...\n")
  all_loocv_predictions <- list()
  
  for (i in seq_len(n_samples)) {
    cat(paste("DEBUG: LOOCV fold", i, "of", n_samples, "...\n"))
    train_idx <- setdiff(seq_len(n_samples), i)
    test_idx  <- i
    
    # Subset data
    trainX <- Data[train_idx, , drop = FALSE]
    trainY <- y[train_idx]
    testX  <- Data[test_idx, , drop = FALSE]
    testY  <- y[test_idx]
    
    # Train model with inner CV (though for LOOCV, it's somewhat trivial)
    modelFit <- train_with_inner_cv(trainX, trainY, algorithm, innerControl)
    
    # Predict on hold-out
    pred_raw  <- predict(modelFit, newdata = testX)
    pred_prob <- predict(modelFit, newdata = testX, type = "prob")
    
    # Store
    fold_res <- data.frame(
      sample       = rownames(testX),
      ActualClass  = testY,
      Predicted    = pred_raw,
      stringsAsFactors = FALSE
    )
    fold_res <- cbind(fold_res, pred_prob)
    all_loocv_predictions[[i]] <- fold_res
    
    # External predictions (optional)
    if (exists("data.external.validation") && !is.null(data.external.validation)) {
      cat(paste("DEBUG: Generating external predictions for fold", i, "...\n"))
      ext_prob <- predict(modelFit, newdata = data.external.validation, type = "prob") %>%
        as.data.frame()
      ext_prob$PredictedClass <- predict(modelFit, newdata = data.external.validation, type = "raw")
      ext_prob$sample <- rownames(data.external.validation)
      External.kfold[[i]] <- ext_prob
    }
  }
  
  # Combine LOOCV predictions
  all_loocv_df <- do.call(rbind, all_loocv_predictions)
  nested_cv_results$predictions <- all_loocv_df
  
  # Place them in All.kFold similarly to your original code
  # We'll just store them as a single "iteration" with N folds
  All.kFold[[1]] <- all_loocv_predictions
  
  # Confusion matrix for overall LOOCV
  cat("DEBUG: Computing confusion matrix for LOOCV...\n")
  cm <- confusionMatrix(
    factor(all_loocv_df$Predicted, levels = levels(y)),
    factor(all_loocv_df$ActualClass, levels = levels(y))
  )
  nested_cv_results$performance <- cm
  cat("DEBUG: LOOCV complete.\n")
  
} else {
  # (K) Otherwise, Nested CV
  cat("DEBUG: Starting nested cross-validation...\n")
  cat(paste("DEBUG: Outer folds =", outer_folds, "| Inner folds =", inner_folds, "\n"))
  # Create outer folds
  set.seed(123)
  cat("DEBUG: Creating outer fold indices with caret::createFolds...\n")
  outerFoldIndices <- createFolds(y, k = outer_folds, list = TRUE, returnTrain = TRUE)
  
  # We'll store predictions from each outer fold
  outer_predictions <- list()
  
  for (fold_idx in seq_along(outerFoldIndices)) {
    cat(paste("DEBUG: Outer fold", fold_idx, "of", outer_folds, "...\n"))
    train_indices <- outerFoldIndices[[fold_idx]]
    test_indices  <- setdiff(seq_len(n_samples), train_indices)
    
    # Subset data
    trainX <- Data[train_indices, , drop = FALSE]
    trainY <- y[train_indices]
    testX  <- Data[test_indices, , drop = FALSE]
    testY  <- y[test_indices]
    
    # Overwrite the innerControl number for the inner folds
    # (We've already done innerControl but let's be explicit.)
    cat("DEBUG: Re-defining innerControl for the inner loop...\n")
    innerControl_fold <- trainControl(
      method         = "repeatedcv",
      number         = inner_folds,
      repeats        = 3,
      classProbs     = TRUE,
      search         = "grid",
      summaryFunction= multiClassSummary,
      verboseIter    = FALSE
    )
    
    # Train with inner CV
    cat("DEBUG: Commencing inner CV for hyperparameter tuning...\n")
    modelFit <- train_with_inner_cv(trainX, trainY, algorithm, innerControl_fold)
    
    # Predict on outer test
    cat("DEBUG: Predicting on outer test fold...\n")
    pred_raw  <- predict(modelFit, newdata = testX)
    pred_prob <- predict(modelFit, newdata = testX, type = "prob")
    
    fold_res <- data.frame(
      sample       = rownames(testX),
      ActualClass  = testY,
      Predicted    = pred_raw,
      stringsAsFactors = FALSE
    )
    fold_res <- cbind(fold_res, pred_prob)
    
    # Store
    outer_predictions[[fold_idx]] <- fold_res
    All.kFold[[fold_idx]] <- list(fold_res)
    
    # External predictions if applicable
    if (exists("data.external.validation") && !is.null(data.external.validation)) {
      cat(paste("DEBUG: Generating external predictions for outer fold", fold_idx, "...\n"))
      ext_prob <- predict(modelFit, newdata = data.external.validation, type = "prob") %>%
        as.data.frame()
      ext_prob$PredictedClass <- predict(modelFit, newdata = data.external.validation, type = "raw")
      ext_prob$sample <- rownames(data.external.validation)
      
      External.kfold[[fold_idx]] <- ext_prob
    }
  }
  
  cat("DEBUG: Combining outer fold predictions...\n")
  all_outer_preds <- do.call(rbind, outer_predictions)
  nested_cv_results$outer_fold_predictions <- all_outer_preds
  
  cat("DEBUG: Computing confusion matrix for entire outer folds...\n")
  cm <- confusionMatrix(
    factor(all_outer_preds$Predicted, levels = levels(y)),
    factor(all_outer_preds$ActualClass, levels = levels(y))
  )
  nested_cv_results$outer_performance <- cm
  cat("DEBUG: Nested CV complete.\n")
}

# (L) Train final model on all data with best approach from the inner loop
cat("DEBUG: Training final model on ALL data...\n")
final_fit <- train_with_inner_cv(Data, y, algorithm, innerControl)
final_model <- final_fit
assign("final_model", final_model, envir = .GlobalEnv)
cat("DEBUG: Final model training completed.\n")

# (M) Final external predictions, if applicable
external_results <- NULL
if (exists("data.external.validation") && !is.null(data.external.validation)) {
  cat("DEBUG: Now predicting final external results with 'final_model'...\n")
  ext_raw  <- predict(final_model, newdata = data.external.validation)
  ext_prob <- predict(final_model, newdata = data.external.validation, type = "prob")
  
  ext_df <- data.frame(
    sample    = rownames(data.external.validation),
    Predicted = ext_raw,
    stringsAsFactors = FALSE
  )
  ext_df <- cbind(ext_df, ext_prob)
  
  external_results <- ext_df
  assign("external_results", external_results, envir = .GlobalEnv)
  cat("DEBUG: Finished external predictions.\n")
} else {
  cat("DEBUG: 'data.external.validation' is NULL or does not exist; skipping external predictions.\n")
}

# (N) Place final objects in the global environment
cat("DEBUG: Storing All.kFold, External.kfold, nested_cv_results in .GlobalEnv...\n")
if (length(External.kfold) == 0) {
  External.kfold <- NULL
}
assign("All.kFold", All.kFold, envir = .GlobalEnv)
assign("External.kfold", External.kfold, envir = .GlobalEnv)
assign("nested_cv_results", nested_cv_results, envir = .GlobalEnv)

# (O) Stop parallel cluster
cat("DEBUG: Stopping parallel cluster...\n")
stopCluster(cl)
cat("DEBUG: All done with final_classifier.R!\n")
