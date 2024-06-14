# logistic regression to predict scz status

# load packages ----

library(tidyverse) # data manipulation
library(matrixStats) # for rowSds
library(ggplot2) # for plotting
library(RColorBrewer) # colors to make heatmaps
library(pheatmap) # for heatmaps
library(dplyr) # for data manipulation
library(heatmaply) # interactive heatmaps
library(plotly) # for interactive plots
library(edgeR) # for differential expression
library(limma) # for differential expression
library(genefilter) # for filtering genes
library(gt) # for interactive tables
library(DT) # for interactive tables
library(cowplot) # multiple plot in the same figure
library(readxl) # read excel files
library(glmnet) # for LASSO
library(pROC) # for ROC curve

# Functions ----

predict_status <- function(expr_data) {
    # Convert the status to a binary factor
    expr_data$status <- ifelse(expr_data$status == "scz", 1, 0)

    print(table(expr_data$status))

    # Split the data into predictors and response
    predictors <- expr_data[, -ncol(expr_data)] # Exclude the last column (status)
    response <- expr_data[, ncol(expr_data)] # The status column

    # Prepare the matrix of predictors and response vector
    X <- as.matrix(predictors)
    y <- as.numeric(response) - 1 # Ensure binary status is 0 and 1

    # Fit the LASSO model
    set.seed(123) # For reproducibility
    cv_lasso <- cv.glmnet(X, y, family = "binomial", alpha = 1)

    # Optimal lambda that gives minimum mean cross-validated error
    opt_lambda <- cv_lasso$lambda.min

    # Fit the final model using the optimal lambda
    lasso_model <- glmnet(X, y, family = "binomial", alpha = 1, lambda = opt_lambda)

    return(list(lasso_model = lasso_model, opt_lambda = opt_lambda))
}

evaluate_the_model <- function(geo_database, expr_data_train, sample_info) {
    path <- file.path(
        "/Users", "herutdor", "Documents", "Research", "Schizoprenia_biomarkers",
        paste0(geo_database, "_data"), geo_database
    )

    expr_data_test <- readRDS(file = paste0(path, "_expr_data.rds"))
    sample_info_test <- readRDS(file = paste0(path, "_sample_info.rds"))

    # print the nummber of scz and control in the test data set
    print(table(sample_info_test$status))

    ### keep only genes that exsist in the test data set ----
    expr_data_train <- expr_data_train[rownames(expr_data_train) %in% rownames(expr_data_test), ]
    expr_data_test <- expr_data_test[rownames(expr_data_test) %in% rownames(expr_data_train), ]

    # check that the 2 data set have the same genes
    all(rownames(expr_data_train) == rownames(expr_data_test))

    ### scale the data ----
    expr_data_train <- as.data.frame((apply(expr_data_train, 1, function(x) (x - mean(x)) / sd(x))))
    expr_data_test <- as.data.frame((apply(expr_data_test, 1, function(x) (x - mean(x)) / sd(x))))

    ### Add the status column to the data ----
    sample_info_train <- sample_info[match(rownames(expr_data_train), rownames(sample_info)), ]

    sample_info_test <- sample_info_test[match(rownames(expr_data_test), rownames(sample_info_test)), ]

    expr_data_train <- cbind(expr_data_train, status = sample_info_train$status)

    model_res <- predict_status(expr_data_train)
    lasso_model <- model_res$lasso_model
    opt_lambda <- model_res$opt_lambda

    # check the model performance on a second data set ----

    # Predict the status using the final model
    predictions_prob <- predict(lasso_model,
        newx = data.matrix(expr_data_test),
        type = "response",
        s = opt_lambda
    )

    # Convert probabilities to binary outcomes based on a threshold
    predictions <- ifelse(predictions_prob > 0.5, 1, 0)
    # predictions <- ifelse(predictions_prob > n_control / (n_scz + n_control), 1, 0)
    response <- ifelse(sample_info_test$status == "scz", 1, 0)

    # Confusion matrix
    conf_matrix <- table(response, predictions)
    print(conf_matrix)

    return(conf_matrix)
}

# Function to perform the binomial test
calculate_binom_p_value <- function(tp, tp_fp, chance_level) {
    # Perform the binomial test
    test_result <- binom.test(tp, tp_fp, p = chance_level, alternative = "greater")

    # Return the p-value
    return(test_result$p.value)
}

set.seed(123) # for reproducibility

# load data ----
geo_database <- "GSE38484"
path <- file.path(
    "/Users", "herutdor", "Documents", "Research", "Schizoprenia_biomarkers",
    paste0(geo_database, "_data"), geo_database
)

expr_data <- readRDS(file = paste0(path, "_expr_data.rds"))
sample_info <- readRDS(file = paste0(path, "_sample_info.rds"))
genes_info <- readRDS(file = paste0(path, "_genes_info.rds"))
ribosomal_genes_id <- readRDS(file = paste0(path, "_ribosomal_genes_id.rds"))
ubl_genes_id <- readRDS(file = paste0(path, "_ubl_genes_id.rds"))

# Train the model based on the ubl and ribosome genes found in the first data set ----

expr_data <- expr_data[rownames(expr_data) %in% c(ribosomal_genes_id, ubl_genes_id), ]

### import the test data set ----

conf_matrix_temp <- matrix(0, nrow = 2, ncol = 2)
conf_matrix_all <- matrix(0, nrow = 2, ncol = 2)
geo_databases <- c("GSE27383", "GSE38481", "GSE18312", "GSE48072")

conf_matrix_ppv <- matrix(0, nrow = length(geo_databases), ncol = 2)
rownames(conf_matrix_ppv) <- geo_databases
colnames(conf_matrix_ppv) <- c("TP", "FP")

for (geo_database in geo_databases) {
    conf_matrix_temp <- evaluate_the_model(geo_database, expr_data, sample_info)
    print(conf_matrix_temp)

    conf_matrix_ppv[geo_database, "TP"] <- conf_matrix_temp[2, 2]
    conf_matrix_ppv[geo_database, "FP"] <- conf_matrix_temp[1, 2]
    conf_matrix_all <- conf_matrix_all + conf_matrix_temp
}

print(conf_matrix_all)

# accuracy
accuracy <- sum(diag(conf_matrix_all)) / sum(conf_matrix_all)
print(accuracy)

# calculate PPV
ppv <- conf_matrix_all[2, 2] / sum(conf_matrix_all[, 2])
print(ppv)

# Convert conf_matrix_ppv to a data frame
conf_matrix_ppv <- as.data.frame(conf_matrix_ppv)

# Add a new row with the sum of TP and FP
conf_matrix_ppv["Combined", ] <- c(sum(conf_matrix_ppv$TP), sum(conf_matrix_ppv$FP))

conf_matrix_ppv$Total <- conf_matrix_ppv$TP + conf_matrix_ppv$FP

# Calculate the PPV (Positive Predictive Value) for each dataset
conf_matrix_ppv$PPV <- conf_matrix_ppv$TP / (conf_matrix_ppv$TP + conf_matrix_ppv$FP)
conf_matrix_ppv$PPV_random <- c(43 / (43 + 29), 15 / (15 + 22), 13 / (13 + 8), 35 / (35 + 31), 106 / (106 + 90))

# Apply the function to each dataset
conf_matrix_ppv$p_value <- mapply(calculate_binom_p_value, conf_matrix_ppv$TP, conf_matrix_ppv$Total, chance_level = conf_matrix_ppv$PPV_random)

saveRDS(conf_matrix_ppv, file = paste0(path, "_conf_matrix.rds"))
