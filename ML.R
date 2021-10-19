library(caret)
library(e1071)


ml_all <- function(df_list, indices=NULL, filename, df_test_list = NULL){
  
  # The function is generalized to either accept a list of dataframes and the test indices,
  # or a list of train dataframes (in place of df_list) and list of test dataframes without the indices
  
  outputs <- list()
  
  for (i in 1:length(df_list)){
    df_name <- str_extract(names(df_list)[[i]], "[a-z]+")
    print(df_name)
    full_name <- names(df_list)[[i]]
    
    if(is.null(df_test_list)){
      df <- df_list[[i]]
      index <- indices[[df_name]]
      
      df_test <- datasets[[df_name]][index, ]
      df_train <- df[-index, ]
    } else {
      df_train <- df_list[[i]]
      df_test <- df_test_list[[df_name]]
    }
    
    print(full_name)
    print(df_name)
    print(dim(df_train))
    
    print(sum(is.na(df_train)))
    print(sum(is.na(df_test)))
    
    df_train <- rename_col(df_train, "df.Output", "Output")
    
    if(dim(df_train)[1] == 0){
      output <- paste(full_name, "-", "-", sep = "_")
    } else if(is.factor(df_train$Output)){
      tic <- Sys.time()
      result <- paste(unlist(ml_classification(df_train, df_test)), collapse = "_")
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      output <- paste(full_name, result, "classification", round(time, 4), sep = "_")
    } else {
      tic <- Sys.time()
      result <- paste(unlist(ml_regression(df_train, df_test)), collapse = "_")
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      output <- paste(full_name, result, "regression", round(time, 4), sep = "_")
    }
    print(output)
    output <- str_split(output, "_")
    outputs <- c(outputs, output)
    
    write.table(t(unlist(output)), paste(filename, "csv", sep = "."), sep = ", ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
  return(outputs)
}


ml_classification <- function(train.data, test.data) {
  library(randomForest)
  library(caret)
  library(psych)
  library(doParallel)
  
  train_method = "rf"
  
  n_cores <- detectCores()
  cl <- makeCluster(n_cores%/%2)
  registerDoParallel(cl)
  print(paste("Number of cores in use:", getDoParWorkers(), sep = " "))
  
  fit <- train(
    y = train.data$Output,
    x = select(train.data,-c(Output)),
    data = train.data,
    method = train_method,
    trControl = trainControl("cv", number = 10, allowParallel = T)
  )
  
  stopCluster(cl)
  
  predTest <- predict(fit, test.data)
  
  # Checking classification accuracy
  acc <- mean(predTest == test.data$Output)
  
  # Compute kapa
  res.k <- cohen.kappa(data.frame(predTest, test.data$Output))$kappa
  return(list(res1 = acc, res2 = res.k))
}

ml_regression <- function(train.data, test.data){
  library(rpart)
  library(doParallel)
  
  # print(tail(test.data), 5)
  train_method <- "lm"
  n <- names(train.data)
  f <- as.formula(paste("Output ~", paste(n[!n %in% "Output"], collapse = " + ")))
  #f <- reformulate(setdiff(colnames(train.data), "Output"), response="Output")
  
  n_cores <- detectCores()
  cl <- makeCluster(n_cores%/%2)
  registerDoParallel(cl)
  print(paste("Number of cores in use:", getDoParWorkers(), sep = " "))
  fit <- train(
    f, data = train.data, method = train_method,
    trControl = trainControl("cv", number = 10, allowParallel = T)
  )
  
  stopCluster(cl)
  
  predictions <- fit %>% predict(test.data)
  
  rmse <- RMSE(predictions, test.data$Output)
  r2 <- R2(predictions, test.data$Output)
  
  return(list(res1 = rmse, res2 = r2))
}


ensemble_ml <- function(df_list, test_dfs, filename){
  library(mice)
  library(randomForest)
  library(caret)
  library(psych)
  library(doParallel)
  
  for (i in 1:length(df_list)){
    df_name <- str_extract(names(df_list)[[i]], "[a-z]+")
    full_name <- names(df_list)[[i]]
    print(paste(full_name, " processing...", sep = ""))
    imp <- df_list[[i]]
    
    df <- complete(imp, 1)
    df <- rename_col(df, "df.Output", "Output")
    classification <- is.factor(df$Output)
    
    df_test <- test_dfs[[df_name]]
    
    predictions <- list()
    predictions_reg <- {}
    r2_list <- list()
    rmse_list <- list()
    tic <- Sys.time()
    
    for (j in 1:imp$m){
      df_train <- complete(imp, j)
      df_train <- rename_col(df_train, "df.Output", "Output")
      if(classification){
        n_cores <- detectCores()
        cl <- makeCluster(n_cores%/%2)
        registerDoParallel(cl)
        print(paste("Number of cores in use:", getDoParWorkers(), sep = " "))
        
        fit <- train (
          y = df_train$Output,
          x = select(df_train,-c(Output)),
          data = df_train,
          method = 'rf',
          trControl = trainControl("cv", number = 10, allowParallel = TRUE)
        )
        
        stopCluster(cl)
        
        predictions[[j]] <- predict(fit, df_test)
        print(paste("Prediction ", j, ": ", sep = ""))
        #print(predictions[[j]][1:400])
      } else {
        n <- names(df_train)
        f <- as.formula(paste("Output ~", paste(n[!n %in% "Output"], collapse = " + ")))
        
        n_cores <- detectCores()
        cl <- makeCluster(n_cores%/%2)
        registerDoParallel(cl)
        print(paste("Number of cores in use:", getDoParWorkers(), sep = " "))
        
        fit <- train(
          f, data = df_train, method = "lm",
          trControl = trainControl("cv", number = 10, allowParallel = TRUE),
          tuneLength = 10
        )
        
        stopCluster(cl)
        
        predictions_reg[[j]] <- fit %>% predict(df_test)
      }
    }
    
    if(classification){
      predictions <- as.data.frame(matrix(unlist(predictions), nrow = length(unlist(predictions[1]))))
      
      predictions <- as.data.frame(predictions, stringsAsFactors = TRUE)
      ensemble_pred <- {}
      for (i in 1:nrow(predictions)) {
        tt <- table(as.matrix(predictions[i, ]))
        ensemble_pred[[i]] <- names(tt[which.max(tt)])
      }
      ensemble_pred <- as.factor(unlist(ensemble_pred))
      accuracy <- mean(ensemble_pred == df_test$Output)
      kappa <- cohen.kappa(data.frame(ensemble_pred, df_test$Output))$kappa
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      output <- paste(full_name, accuracy, kappa, round(time, 4), "classification", sep = "_")
      print(paste("Classification:", accuracy, kappa, round(time, 4), sep = " "))
    } else {
      predictions_reg <- as.data.frame(matrix(unlist(predictions_reg), nrow = length(unlist(predictions_reg[1]))))
      
      ensemble_reg <- {}
      for (i in 1:nrow(predictions_reg)) {
        ensemble_reg[[i]] <- mean(as.matrix(predictions_reg[i, ]))
      }
      ensemble_reg <- as.numeric(unlist(ensemble_reg))
      rmse <- RMSE(ensemble_reg, df_test$Output)
      r2 <- R2(ensemble_reg, df_test$Output)
      
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      print(paste("Regression:", rmse, r2, round(time, 4), sep = " "))
      output <- paste(full_name, rmse, r2, round(time, 4), "regression", sep = "_")
    }
    write.table(t(unlist(output)), paste(filename, "csv", sep = "."), sep = ", ", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}