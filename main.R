## Load datasets

dataset_names <- c("asn", "ccp", "concrete", "egs", "magic", "onp", "pageblocks", "pendigits", "spambase", "superconductor")
datasets <- list()
for (i in 1:10){
  df_name <- paste("Datasets/", dataset_names[i], ".csv", sep = "")
  
  datasets[[dataset_names[i]]] <- read.csv(df_name, header = TRUE)
}

#Convert target variables to factor
datasets[["pageblocks"]]$Output <- as.factor(datasets[["pageblocks"]]$Output)
datasets[["onp"]]$Output <- as.factor(datasets[["onp"]]$Output)
datasets[["pendigits"]]$Output <- as.factor(datasets[["pendigits"]]$Output)
datasets[["spambase"]]$Output <- as.factor(datasets[["spambase"]]$Output)
datasets[["magic"]]$Output <- as.factor(datasets[["magic"]]$Output)


## Perform Amputation and Imputation

source("Imputation_ML.R")
source("logreg.R")
test_indices <- generate_test_indices(datasets)
train_test <- train_test_split(datasets, test_indices)
ampute_init()
amputed_datasets <- ampute_dfs(train_test$train)
source("ML.R")
imputed_datasets <- impute_dfs(amputed_datasets, "all")
names(imputed_datasets) <- names(amputed_datasets)


## Machine Learning

ml_all(df_list = imputed_datasets, df_test_list = train_test$test, filename = "ml_cca")
ml_all(df_list = imputed_datasets, df_test_list = train_test$test, filename = "ml_mean")
ml_all(df_list = imputed_datasets, df_test_list = train_test$test, filename = "ml_hotdeck")
ml_all(df_list = imputed_datasets, df_test_list = train_test$test, filename = "ml_regression")
ensemble_ml(df_list = imputed_datasets, df_test_list = train_test$test, filename = "ml_mean")


ml_all(df_list = train_test$train, df_test_list = train_test$test, filename = "gold_standard_ml")



## Compute Deviations

print_true_missing_ratio(amputed_datasets)

results_mean <- read.csv("ml_mean.csv", header = FALSE)
names(results_mean) <- c("dataset", "Mechanism", "Perc_Missing", "Metric_1", "Metric_2", "ML_Type", "Time")
results_cca <- read.csv("ml_cca.csv", header = FALSE)
names(results_cca) <- c("dataset", "Mechanism", "Perc_Missing", "Metric_1", "Metric_2", "ML_Type", "Time")
results_hotdeck <- read.csv("ml_hotdeck.csv", header = FALSE)
names(results_hotdeck) <- c("dataset", "Mechanism", "Perc_Missing", "Metric_1", "Metric_2", "ML_Type", "Time")
results_regression <- read.csv("ml_regression.csv", header = FALSE)
names(results_regression) <- c("dataset", "Mechanism", "Perc_Missing", "Metric_1", "Metric_2", "ML_Type", "Time")
results_mi <- read.csv("ml_mi.csv", header = FALSE)
names(results_mi) <- c("dataset", "Mechanism", "Perc_Missing", "Metric_1", "Metric_2", "Time", "ML_Type")

results <- list(mean = results_mean, cca = results_cca, hotdeck = results_hotdeck, regression = results_regression, mi = results_mi)

gold_standard <- read.csv("gold_standard_ml.csv", header = FALSE)
names(gold_standard) <- c("dataset", "Metric_1", "Metric_2", "ML_Type", "Time")

results <- lapply(results, compute_errors, gold_standard)

results_r <- data.frame(results$mean, results)

# Get real missingness ratios
df_miss_rates <- read.csv("missingness_ratios.csv", header = FALSE)
names(df_miss_rates) <- c("dataset", "mechanism", "ideal", "real")
results <- lapply(results, function(x) cbind(x, df_miss_rates$real))
results <- lapply(results, rename_col, "df_miss_rates$real", "real_miss_rate")

imputation_timings <- read.csv("imputation_timings.csv", header = FALSE)
names(imputation_timings) <- c("dataset", "mechanism", "ideal", "time", "imputation")
df_real_timings <- unique(left_join(imputation_timings, df_miss_rates))

results$cca$imputation <- "cca"
results$mean$imputation <- "mean"
results$hotdeck$imputation <- "hotdeck"
results$regression$imputation <- "regression"
results$mi$imputation <- "mi"

results$mi <- subset(results$mi, select=c(dataset:Metric_2, ML_Type, Time, Metric_1_err:imputation)) # rearrange MI columns to match others

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
  geom_text(aes(0,h,label = "g s", vjust = -1, hjust = 1)) +
  theme_bw() + ggtitle("CONCRETE")
ggsave("ccp-mcar.png")


mechanisms <- c("MCAR", "MAR", "MNAR")
names <- names(datasets)

save_plots_time(mechanisms, names)
