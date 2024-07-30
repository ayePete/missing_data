ampute_dfs <- function(dfs){
  log_output <- paste(Sys.time(), "------------------------------- AMPUTATION: BEGIN -------------------------------")
  write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  #write.csv(c(Sys.time(), "AMPUTATION: BEGIN"), "log.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  print("AMPUTATION: BEGIN")
  library(mice)
  amputed_dfs <- list()
  names <- vector()
  for (i in 1:length(dfs)) {
    for (prop in props) {
      for (mech in mechs){
        df <- data.frame(dfs[[i]])
        df_name <- paste(names(dfs)[i], mech, prop, sep = "_")
        names <- c(names, df_name)
        to_ampute <- df[ , setdiff(names(df), c("Output"))]
        #to_ampute <- df
        amped_df <- ampute(data = to_ampute, prop = prop, bycases = FALSE, patterns = diag(ncol(to_ampute)), mech = mech)$amp
        amped_df <- data.frame(amped_df, df$Output)
        names(amped_df)[names(amped_df) == 'df.Output'] <- 'Output'
        #str(amped_df)
        if(is.factor(df$Output))
          amped_df$Output <- as.factor(amped_df$Output)
        amputed_dfs <- c(amputed_dfs, list(amped_df))
        log_output <- paste(Sys.time(), df_name)
        write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
        print(df_name)
      }
    }
  }
  names(amputed_dfs) <- names
  log_output <- paste(Sys.time(), "---------------------------- AMPUTATION: END ------------------------------")
  write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  return(amputed_dfs)
}

impute_dfs <- function(dfs, option){
  library(mice)
  library(VIM)
  library(tictoc)
  
  imputed_dfs <- list()
  
  
  if(option == "cca" || option == "all"){
    log_output <- paste(Sys.time(), "--------------- CCA: BEGIN ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    print("_________________________")
    print("Complete Case Analysis")
    print("_________________________")
    print("  ")
    
    imputed_cca <- list()
    
    for (i in 1:length(dfs)) {
      df <- dfs[[i]]
      
      name_list <- str_split(names(dfs)[i], "_")
      print(name_list)
      tic <- Sys.time()
      imp <- df[complete.cases(df), ]
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      if(is.factor(df$Output))
        imp$Output <- as.factor(imp$Output)
      imputed_cca <- c(imputed_cca, list(imp))
      time.dat <- paste(names(dfs)[i], round(time, 4), option, sep = "_")
      time.dat <- str_split(time.dat, "_")
      
      write.table(t(unlist(time.dat)), "timings.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      log_output <- c(paste(Sys.time()), time.dat)
      write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    if (option != "all"){
      return(imputed_cca)
    } else {
      imputed_dfs["cca"] <- imputed_cca
    }
    
    log_output <- paste(Sys.time(), "--------------- CCA: END ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
  if(option == "mean" || option == "all"){
    log_output <- paste(Sys.time(), "--------------- MEAN IMPUTATION: BEGIN ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    print("  ")
    print("_________________________")
    print("Mean Imputation")
    print("_________________________")
    print("  ")
    
    imputed_mean <- list()
    for (i in 1:length(dfs)) {
      
      
      df <- dfs[[i]]
      
      name_list <- str_split(names(dfs)[i], "_")
      print(name_list)
      tic <- Sys.time()
      imp <- mean_impute(df)
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      if(is.factor(df$Output))
        imp$Output <- as.factor(imp$Output)
      imputed_mean <- c(imputed_mean, list(imp))
      time.dat <- paste(names(dfs)[i], round(time, 4), option, sep = "_")
      time.dat <- str_split(time.dat, "_")
      
      
      
      write.table(t(unlist(time.dat)), "timings.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      log_output <- c(paste(Sys.time()), time.dat)
      write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    if (option != "all"){
      return(imputed_mean)
    } else {
      imputed_dfs[["mean"]] <- imputed_mean
    }
    
    log_output <- paste(Sys.time(), "--------------- MEAN IMPUTATION: END ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
  if(option == "hotdeck" || option == "all"){
    log_output <- paste(Sys.time(), "--------------- HOTDECK IMPUTATION: BEGIN ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    print("  ")
    print("_________________________")
    print("Hotdeck Imputation")
    print("_________________________")
    print("  ")
    
    imputed_hotdeck <- list()
    for (i in 1:length(dfs)) {
      
      
      df <- data.frame(dfs[[i]])
      
      
      name_list <- str_split(names(dfs)[i], "_")
      print(name_list)
      tic <- Sys.time()
      #imp <- hot.deck(df, m = 1)
      imp <- hotdeck(df, imp_var = FALSE)
      
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      if(is.factor(df$Output))
        imp$Output <- as.factor(imp$Output)
      imputed_hotdeck <- c(imputed_hotdeck, list(imp))
      time.dat <- paste(names(dfs)[i], round(time, 4), option, sep = "_")
      time.dat <- str_split(time.dat, "_")
      
      write.table(t(unlist(time.dat)), "timings.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      log_output <- c(paste(Sys.time()), time.dat)
      write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      
    }
    if (option != "all"){
      return(imputed_hotdeck)
    } else {
      imputed_dfs[["hotdeck"]] <- imputed_hotdeck
    }
    
    log_output <- paste(Sys.time(), "--------------- HOTDECK IMPUTATION: END ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
  if(option == "regression" || option == "all"){
    log_output <- paste(Sys.time(), "--------------- REGRESSION IMPUTATION: BEGIN ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    tic()
    print("  ")
    print("_________________________")
    print("Regression Imputation")
    print("_________________________")
    print("  ")
    
    imputed_regression <- list()
    for (i in 1:length(dfs)) {
      
      df <- dfs[[i]]
      
      name_list <- str_split(names(dfs)[i], "_")
      print(name_list)
      tic <- Sys.time()
      imp <- complete(mice(df, m = 1, defaultMethod = c("norm","logreg","polyreg","polr")), 1)
      #imp <-  complete(mice(df, m = 1, defaultMethod = c("cart","logreg","polyreg","polr"), predictorMatrix = quickpred(df, mincor = 0.25, minpuc = 0.25)), 1)
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      if(is.factor(df$Output))
        imp$Output <- as.factor(imp$Output)
      imputed_regression <- c(imputed_regression, list(imp))
      time.dat <- paste(names(dfs)[i], round(time, 4), option, sep = "_")
      time.dat <- str_split(time.dat, "_")
      
      write.table(t(unlist(time.dat)), "timings.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      log_output <- c(paste(Sys.time()), time.dat)
      write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    if (option != "all"){
      return(imputed_regression)
    } else {
      imputed_dfs[["regression"]] <- imputed_regression
    }
    
    log_output <- paste(Sys.time(), "--------------- REGRESSION IMPUTATION: END ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
  if(option == "mi" || option == "all"){
    
    log_output <- paste(Sys.time(), "--------------- MULTIPLE IMPUTATION: BEGIN ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
    print("  ")
    print("_________________________")
    print("Multiple Imputation")
    print("_________________________")
    print("  ")
    
    imputed_mi <- list()
    for (i in 1:length(dfs)) {
      
      
      df <- dfs[[i]]
      name_list <- str_split(names(dfs)[i], "_")
      print(name_list)
      tic <- Sys.time()
      imp <- mice(df, m = 5)
      #imp <- mice(df, m = 5, defaultMethod = c("cart","logreg","polyreg","polr"), predictorMatrix = quickpred(df, mincor = 0.25, minpuc = 0.25) )
      time <- difftime(Sys.time(), tic, units = "secs")[[1]]
      imputed_mi <- c(imputed_mi, list(imp))
      time.dat <- paste(names(dfs)[i], round(time, 4), option, sep = "_")
      time.dat <- str_split(time.dat, "_")
      
      write.table(t(unlist(time.dat)), "imputation_timings.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      log_output <- c(paste(Sys.time()), time.dat)
      write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    if (option != "all"){
      return(imputed_mi)
    } else {
      imputed_dfs[["_mi"]] <- imputed_mi
    }
    
    log_output <- paste(Sys.time(), "--------------- MULTIPLE IMPUTATION: END ----------------")
    write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  return(imputed_dfs)
}

generate_test_indices <- function(df_list){
  index_list <- list()
  for (i in 1:length(df_list)){
    x <- sample(nrow(df_list[[i]]), round(0.3 * nrow(df_list[[i]]), 0))
    index_list[[i]] <- x
  }
  names(index_list) <- names(df_list)
  return(index_list)
}

train_test_split <- function(df_list, index_list=NULL){
  train <- list()
  test <- list()
  
  # pendigits <- cbind(pendigits, split1)
  # data_train <- pendigits[split1 == 0,]
  # data_test <- pendigits[split1 == 1,]
  for (i in 1:length(df_list)){
    split1<- sample(c(rep(0, 0.7 * nrow(df_list[[i]])), rep(1, 0.3 * nrow(df_list[[i]]))))
    if(length(split1) != nrow(df_list[[i]])){
      split1 <- c(split1,0)
    }
    df_list[[i]] <- cbind(df_list[[i]], split1)
    train[[i]] <- df_list[[i]][split1 == 0,]
    test[[i]] <- df_list[[i]][split1 == 1,]
    df_list[[i]] <- select(df_list[[i]], -split1)
    #df_list[[i]]$split1 <- NULL
    #df_list[[i]] <- data.frame(df_list[[i]])
    #test[[i]] <- data.frame(df_list[[i]][index_list[[i]], ])
    #train[[i]] <- data.frame(df_list[[i]][-index_list[[i]], ])
  }
  
  names(test) <- names(df_list)
  names(train) <- names(df_list)
  return(list(train = train, test = test))
}

print_true_missing_ratio <- function(dfs){
  names <- list()
  for (i in 1:length(dfs)){
    ratio <- mean(is.na(dfs[[i]]))
    name <- paste(names(dfs)[i], round(ratio, 4), sep = "_")
    name <- t(unlist(str_split(name, "_")))
    names[[i]] <- paste(name, sep = ",")
    print(names[[i]])
    write.table(t(names[[i]]), "missingness_ratios_3.csv", sep = ",", 
                quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  return(names)
}


rename_col <- function(df, oldname, newname){
  names(df)[names(df) == oldname] <- newname
  return(df)
}


mean_impute <- function(df){
  df <- data.frame(df)
  
  for(i in 1:ncol(df)){
    
    if(is.numeric(df[[i]])){
      #print("Numeric")
      #print(df[[i]])
      df[is.na(df[,i]), i] <- mean(df[,i], na.rm = TRUE)
    } else {
      #print("Factor")
      df[is.na(df[,i]), i] <- Mode(df[,i], na.rm = TRUE)
    }
  }
  
  return(df)
}


Mode <- function(x, na.rm=FALSE) {
  
  # // Source
  # // https://stackoverflow.com/questions/55212746/rcpp-fast-statistical-mode-function-with-vector-input-of-any-type
  # // Author: Ralf Stubner, Joseph Wood
  
  #if(!is.atomic(x) | is.matrix(x)) stop("Mode supports only atomic vectors. Use sapply(*, Mode) instead.")
  
  if (na.rm) 
    x <- x[!is.na(x)]
  
  if (anyNA(x)) 
    # there are NAs, so no mode exist nor frequency
    return(structure(NA_real_, freq = NA_integer_))
  
  if(length(x) == 1L)
    # only one value in x, x is the mode
    return(structure(x, freq = 1L)) 
  
  # we don't have NAs so far, either there were then we've already stopped
  # or they've been stripped above
  res <- calculate_mode(x)
  
  #print(res)
  
  return(res)
}


calculate_mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

save_plots_metric <- function(mechanisms, names){
  for (mech in mechanisms){
    for (name in names){
      h <- gold_standard[gold_standard$dataset == name, 3]
      print(gold_standard[gold_standard$dataset == name, 4])
      if(gold_standard[gold_standard$dataset == name, 4] == "regression"){
        metric <- "R-Squared"
      } else {
        metric <- "Kappa"
      }
      combined_results %>% filter(dataset == name & Mechanism == mech) %>% ggplot(., aes(x = real_miss_rate, y = Metric_2)) + 
        geom_line(aes(color = imputation), size = 1) + 
        geom_hline(yintercept = h, size = 1) + 
        #geom_text(aes(0,h,label = "G S", vjust = -1, hjust = 1, size = 18)) +
        labs(title  = toupper(name), x = "Missingness Rate", y = metric) + 
        theme_bw() + 
        theme(plot.title = element_text(size = 24, hjust = 0.5, face = 'bold'), legend.position = "none", 
              axis.text = element_text(size = 18, face = 'bold'), axis.title = element_text(size = 18, face = 'bold')) #, legend.position = "none")
      filename <- paste(name, mech, sep = "-")
      ggsave(paste(filename, ".png", sep = ""), width = 9, height = 7)
    }
  }
}


save_plots_time <- function(mechanisms, names){
  for (mech in mechanisms){
    for (name in names){
      df_real_timings %>% filter(dataset == name & mechanism == mech) %>% ggplot(., aes(x = real, y = time)) + 
        geom_line(aes(color = imputation), size = 1) + 
        labs(title  = toupper(name), x = "Missingness Rate", y = "Time (s)") + 
        theme_bw() + 
        theme(plot.title = element_text(size = 24, hjust = 0.5, face = 'bold'), legend.position = "none", 
              axis.text = element_text(size = 18, face = 'bold'), axis.title = element_text(size = 18, face = 'bold'))
      filename <- paste(name, mech, "time", sep = "-")
      ggsave(paste(filename, ".png", sep = ""), width = 9, height = 7)
    }
  }
}


save_plots <- function(mechanisms, names){
  for (mech in mechanisms){
    for (name in names){
      if(combined_results$dataset )
      h <- gold_standard[gold_standard$dataset == name, 3]
      combined_results %>% filter(dataset == name & Mechanism == mech) %>% ggplot(., aes(x = real_miss_rate, y = Metric_2)) + 
        geom_line(aes(color = imputation, linetype = imputation)) + 
        geom_hline(yintercept = h) +
        geom_text(aes(0,h,label = "g.s.", vjust = 0, hjust = -11.5)) +
        theme_bw() + ggtitle(toupper(name)) + 
        theme(plot.title = element_text(size = 24, hjust = 0.5, face = 'bold'), legend.position = "none", 
              axis.text = element_text(size = 18, face = 'bold'), axis.title = element_text(size = 18, face = 'bold'))
      filename <- paste(name, mech, sep = "-")
      ggsave(paste(filename, ".emf", sep = ""), width = 9, height = 7)
    }
  }
}



compute_errors <- function(imputations_ls, gs_df){
  imputations_err = {}
  print(names(imputations_ls))
  #gs_df <- as.data.frame(gs_ls)
  #names(gs_df) <- c('data','rmse','r2','ml_type','time(s)')
  for (i in 1:length(imputations_ls)) {
    ls_name <- names(imputations_ls)[i]
    print(paste("List:", ls_name))
    computed = as.data.frame(imputations_ls[[i]])
    print(dim(computed))
    names(computed) <- c('dataset', 'Mechanism', 'Perc_Missing', 'Metric_1', 'Metric_2','ML_Type','Time')
    join_df <- merge(y = gs_df, x = computed, by = "dataset", all.x = TRUE)
    print(str(join_df))
    join_df$Err_1 <- join_df$M_1 - join_df$Metric_1
    join_df$Err_2 <- join_df$M_2 - join_df$Metric_2
    #view(join_df)
    imputations_err[[i]] <- join_df
  }
  names(imputations_err) <- c("mean","cca","hotdeck","regression","mi")
  return(imputations_err)
}