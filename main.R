# Dependencies
library(stringr)
library(tidyverse)
library(mice)

## Load datasets

log_output <- paste(Sys.time(), "INIT...")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

dataset_names <- c("asn","ccp","concrete","egs","magic","avila","pageblocks","pendigits","spambase","superconductor")
datasets <- list()
for (i in 1:10){
  df_name <- paste("Datasets/", dataset_names[i], ".csv", sep = "")
  
  datasets[[dataset_names[i]]] <- read.csv(df_name, header = TRUE)
}

#Convert target variables to factor for classification datasets
datasets[["pageblocks"]]$Output <- as.factor(datasets[["pageblocks"]]$Output)
datasets[["avila"]]$Output <- as.factor(datasets[["avila"]]$Output)
datasets[["pendigits"]]$Output <- as.factor(datasets[["pendigits"]]$Output)
datasets[["spambase"]]$Output <- as.factor(datasets[["spambase"]]$Output)
datasets[["magic"]]$Output <- as.factor(datasets[["magic"]]$Output)


## Perform Amputation and Imputation

props <- c(0.05, 0.25, 0.5, 0.75)
mechs <- c("MCAR","MAR","MNAR")

source("Imputation_ML.R")

#test_indices <- generate_test_indices(datasets)
train_test <- train_test_split(datasets)

# Remove spli1 column from all datasets

train_test$train <- lapply(train_test$train, function(x){
  x$split1 <- NULL
  return(x)
})

train_test$test <- lapply(train_test$test, function(x){
  x$split1 <- NULL
  return(x)
})

#Ampute Datasets
amputed_datasets <- ampute_dfs(train_test$train)
#amputed_datasets <- ampute_dfs(datasets)

#imputed_datasets <- impute_dfs(amputed_datasets, "all")
#names(imputed_datasets) <- names(amputed_datasets)

# Remove spli1 column from all datasets
amputed_datasets <- lapply(amputed_datasets, function(x){
  x$split1 <- NULL
  return(x)
})

train_test$train <- lapply(train_test$train, function(x){
  x$split1 <- NULL
  return(x)
})

train_test$test <- lapply(train_test$test, function(x){
  x$split1 <- NULL
  return(x)
})


log_output <- paste(Sys.time(), "------------------------------- IMPUTATION: BEGIN -------------------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
imputed_datasets_cca <- impute_dfs(amputed_datasets, "cca")
names(imputed_datasets_cca) <- names(amputed_datasets)
save.image("C:\\Users\\User\\OneDrive - University of KwaZulu-Natal\\PhD\\Edited\\workspace_simple_time.RData")

imputed_datasets_mean <- impute_dfs(amputed_datasets, "mean")
names(imputed_datasets_mean) <- names(amputed_datasets)
save.image("C:\\Users\\User\\OneDrive - University of KwaZulu-Natal\\PhD\\Edited\\workspace_simple_time.RData")

imputed_datasets_hotdeck <- impute_dfs(amputed_datasets, "hotdeck")
names(imputed_datasets_hotdeck) <- names(amputed_datasets)
save.image("C:\\Users\\User\\OneDrive - University of KwaZulu-Natal\\PhD\\Edited\\workspace_simple_time.RData")

imputed_datasets_regression <- impute_dfs(amputed_datasets, "regression")
names(imputed_datasets_regression) <- names(amputed_datasets)
save.image("C:\\Users\\User\\OneDrive - University of KwaZulu-Natal\\PhD\\Edited\\workspace_simple_time.RData")

imputed_datasets_mi <- impute_dfs(amputed_datasets, "mi")
names(imputed_datasets_mi) <- names(amputed_datasets)
log_output <- paste(Sys.time(), "------------------------------- IMPUTATION: END -------------------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
save.image("C:\\Users\\User\\OneDrive - University of KwaZulu-Natal\\PhD\\Edited\\workspace_simple_time.RData")

imputed_datasets <- c(imputed_datasets_cca, imputed_datasets_mean, imputed_datasets_hotdeck, imputed_datasets_regression, imputed_datasets_mi)


## Machine Learning
source("ML.R")
log_output <- paste(Sys.time(), "------------------------------- MACHINE LEARNING: BEGIN -------------------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

log_output <- paste(Sys.time(), "------------ MACHINE LEARNING: CCA  ----------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
ml_all(df_list = imputed_datasets_cca, df_test_list = train_test$test, filename = "ml_cca")
log_output <- paste(Sys.time(), "------------ MACHINE LEARNING: MEAN  ----------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
ml_all(df_list = imputed_datasets_mean, df_test_list = train_test$test, filename = "ml_mean")
log_output <- paste(Sys.time(), "------------ MACHINE LEARNING: HOTDECK  ----------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
ml_all(df_list = imputed_datasets_hotdeck, df_test_list = train_test$test, filename = "ml_hotdeck")
log_output <- paste(Sys.time(), "------------ MACHINE LEARNING: REGRESSION  ----------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
ml_all(df_list = imputed_datasets_regression, df_test_list = train_test$test, filename = "ml_regression")
log_output <- paste(Sys.time(), "------------ MACHINE LEARNING: MULTIPLE IMPUTATION  ----------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
ensemble_ml(df_list = imputed_datasets_mi, test_dfs = train_test$test, filename = "ml_mi")

log_output <- paste(Sys.time(), "------------------------------- MACHINE LEARNING: END -------------------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

log_output <- paste(Sys.time(), "------------------------------- MACHINE LEARNING (GS): BEGIN -------------------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
ml_all(df_list = train_test$train, df_test_list = train_test$test, filename = "gold_standard_ml")

log_output <- paste(Sys.time(), "------------------------------- MACHINE LEARNING (GS): END -------------------------")
write.table(t(unlist(log_output)), "log.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)


## Compute Deviations

print_true_missing_ratio(amputed_datasets)

results_mean <- read.csv("ml_mean.csv", header = FALSE, na.strings = '-')
names(results_mean) <- c("dataset","Mechanism","Perc_Missing","Metric_1","Metric_2","ML_Type","Time")
results_cca <- read.csv("ml_cca.csv", header = FALSE, na.strings = '-')
names(results_cca) <- c("dataset","Mechanism","Perc_Missing","Metric_1","Metric_2","ML_Type","Time")
results_hotdeck <- read.csv("ml_hotdeck.csv", header = FALSE, na.strings = '-')
names(results_hotdeck) <- c("dataset","Mechanism","Perc_Missing","Metric_1","Metric_2","ML_Type","Time")
results_regression <- read.csv("ml_regression.csv", header = FALSE, na.strings = '-')
names(results_regression) <- c("dataset","Mechanism","Perc_Missing","Metric_1","Metric_2","ML_Type","Time")
results_mi <- read.csv("ml_mi.csv", header = FALSE, na.strings = '-', sep = "_")
names(results_mi) <- c("dataset","Mechanism","Perc_Missing","Metric_1","Metric_2","Time","ML_Type")

results <- list(mean = results_mean, cca = results_cca, hotdeck = results_hotdeck, regression = results_regression, mi = results_mi)

gold_standard <- read.csv("gold_standard_ml.csv", header = FALSE)
names(gold_standard) <- c("dataset","M_1","M_2","ML_Type","Time")

results <- compute_errors(results, gold_standard)

results_r <- data.frame(results$mean, results)

# Get real missingness ratios
df_miss_rates <- read.csv("missingness_ratios_3.csv", header = FALSE)
names(df_miss_rates) <- c("dataset","mechanism","ideal","real")
results <- lapply(results, function(x) cbind(x, df_miss_rates$real))
results <- lapply(results, rename_col, "df_miss_rates$real","real_miss_rate")

imputation_timings <- read.csv("imputation_timings.csv", header = FALSE)
names(imputation_timings) <- c("dataset","mechanism","ideal","time","imputation")
df_real_timings <- unique(left_join(imputation_timings, df_miss_rates))

results$cca$imputation <- "cca"
results$mean$imputation <- "mean"
results$hotdeck$imputation <- "hotdeck"
results$regression$imputation <- "regression"
results$mi$imputation <- "mi"

#results$mi <- subset(results$mi, select=c(dataset:Metric_2, ML_Type, Time, Metric_1_err:imputation)) # rearrange MI columns to match others

library(data.table)
combined_results <- as.data.frame(rbindlist(results))
combined_results$imputation <- as.factor(combined_results$imputation)



# Generate Plots
results$cca %>% filter(ML_Type == "regression") %>% ggplot(., aes(x = real_miss_rate, y = Metric_2_err)) + 
  geom_line(aes(color = Mechanism, linetype = Mechanism)) + theme_bw() + ggtitle("Time")

df_real_timings %>% filter(dataset == "egs") %>% ggplot(., aes(x = real, y = time)) + 
  geom_line(aes(color = imputation)) + theme_bw() + ggtitle("Time")

h <- gold_standard[gold_standard$dataset == "concrete", 3]
combined_results %>% filter(dataset == "concrete" & Mechanism == "MCAR") %>% ggplot(., aes(x = real_miss_rate, y = Metric_2)) + 
  geom_line(aes(color = imputation, linetype = imputation)) + 
  geom_hline(yintercept = h) +
  geom_text(aes(0,h,label = "g.s.", vjust = 0, hjust = -11.5)) +
  theme_bw() + ggtitle("CONCRETE")
ggsave("ccp-mcar.png")


mechanisms <- c("MCAR","MAR","MNAR")
names <- names(datasets)

save_plots_time(mechanisms, names)

save_plots_metric(mechanisms, names)
