# load packages and functions ----
source("source_code.R")

# Color palette: https://coolors.co/palette/ef476f-ffd166-06d6a0-118ab2-073b4c

# import data ----
geo_database <- "GSE38484"
path <- file.path(
    "/Users", "herutdor", "Documents", "Research", "Schizoprenia_biomarkers",
    paste0(geo_database, "_data"), geo_database
)

expr_data <- readRDS(file = paste0(path, "_expr_data.rds"))
sample_info <- readRDS(file = paste0(path, "_sample_info.rds"))
genes_info <- readRDS(file = paste0(path, "_genes_info.rds"))
diff_scz_data <- readRDS(file = paste0(path, "_diff_scz_data.rds"))
cluster_per_sample <- readRDS(file = paste0(path, "_cluster_per_sample.rds"))
si_score_by_gene <- readRDS(file = paste0(path, "_si_score_by_gene.rds"))
ribosomal_genes_id <- readRDS(file = paste0(path, "_ribosomal_genes_id.rds"))
ubl_genes_id <- readRDS(file = paste0(path, "_ubl_genes_id.rds"))
conf_matrix <- readRDS(file = paste0(path, "_conf_matrix.rds"))
si_threshold <- 0.33

control_samples <- rownames(sample_info[sample_info$status == "control", ])
scz_samples <- rownames(sample_info[sample_info$status == "scz", ])

score_threshold <- quantile(si_score_by_gene$score, si_threshold)
clustering_genes <- rownames(si_score_by_gene[si_score_by_gene$score <= score_threshold, , drop = FALSE])

expr_data_clustering <-
    expr_data[rownames(expr_data) %in% rownames(si_score_by_gene), ]


# create label data frame containing the label for each sample:
# 0 for healthy control, 1 for cluster 1, 2 for cluster 2
label_per_sample <- data.frame(replicate(1, rep(0, nrow(sample_info))))
rownames(label_per_sample) <- rownames(sample_info)
colnames(label_per_sample) <- "label"
label_per_sample$label[rownames(label_per_sample) %in% control_samples] <- 0
label_per_sample$label[rownames(label_per_sample) %in%
    rownames(cluster_per_sample[cluster_per_sample$num_cluster == 1, , drop = FALSE])] <- 1
label_per_sample$label[rownames(label_per_sample) %in%
    rownames(cluster_per_sample[cluster_per_sample$num_cluster == 2, , drop = FALSE])] <- 2

expr_data_ribosome <- expr_data[rownames(expr_data) %in% ribosomal_genes_id, ]
expr_data_ubl <- expr_data[rownames(expr_data) %in% ubl_genes_id, ]

# get the demographic information in each data set ----
geo_databases <- c("GSE38484", "GSE27383", "GSE38481", "GSE18312", "GSE48072")
for (geo_database in geo_databases) {
    temp_path <- file.path(
        "/Users", "herutdor", "Documents", "Research", "Schizoprenia_biomarkers",
        paste0(geo_database, "_data"), geo_database
    )

    # number of scz and control in the data set
    sample_info <- readRDS(file = paste0(temp_path, "_sample_info.rds"))

    # print the demographic information
    cat("\nGEO Database:", geo_database, "\n")

    cat("\nSample Status Table:")
    print(table(sample_info$status))

    # Check if "age" column exists and output mean and standard deviation
    if ("age" %in% colnames(sample_info)) {
        mean_age <- round(mean(sample_info$age), 2)
        sd_age <- round(sd(sample_info$age), 2)

        cat("\nAge Information:\n")
        cat("Mean Age:", mean_age, "\n")
        cat("Standard Deviation of Age:", sd_age, "\n")
    }
}

# Graphs ----

### Fig. 1 ----
# A heatmap of all the clustering genes annotated by label (healthy control, cluster 1 and cluster 2).
# This will demonstrate that one cluster has up-regulation of genes and the second cluster has gene expression
# patterns similar to that of the healthy control.
expr_data_high_si <-
    expr_data[rownames(expr_data) %in% clustering_genes, ]

gene_groups <- data.frame(replicate(1, rep(0, nrow(si_score_by_gene))))
rownames(gene_groups) <- rownames(si_score_by_gene)
colnames(gene_groups) <- "group"
gene_groups$group[rownames(gene_groups) %in% ribosomal_genes_id] <- 1
gene_groups$group[rownames(gene_groups) %in% ubl_genes_id] <- 2

heatmap_clustering_genes(expr_data_high_si, label_per_sample, gene_groups)

### Fig. 2a ----
# PCA classification of the samples with labels indicating the cluster of this sample - ribosome genes
pca_plot(expr_data_ribosome, label_per_sample, "Ribosomal", type = "2D")

### Fig. 2b ----
# PCA classification of the samples with labels indicating the cluster of this sample - ubl genes
pca_plot(expr_data_ubl, label_per_sample, "UPS", type = "2D")

### Fig. 3a ----
# A heatmap of the enrichment genes annotated by clusters - ribosome genes
heatmap_gene_analysis(expr_data_ribosome, cluster_per_sample, "ribosome")
heatmap_gene_analysis_control(expr_data_ribosome, label_per_sample, "ribosome")

### Fig. 3b ----
# A heatmap of the enrichment genes annotated by clusters - ubl genes
heatmap_gene_analysis(expr_data_ubl, cluster_per_sample, "UPS")

### Fig 4. ----
# Barplot of the truly identify scz patients out of the total scz patients in each dataset

conf_matrix$TP_percent <- conf_matrix$TP * 100

# Create a new column for the dataset names from row names
conf_matrix$Dataset <- rownames(conf_matrix)

# Create the bar plot with textures for TP and FP
bar_plot(conf_matrix)

# Supplementary Graphs ----

### Supplementary Fig. 1 ----
info_by_cluster <- merge(sample_info, cluster_per_sample, by = "row.names")
rownames(info_by_cluster) <- info_by_cluster$Row.names
info_by_cluster <- subset(info_by_cluster, select = -c(Row.names, status))

expr_data_ribosome_scz <- expr_data_ribosome[, scz_samples]

lr_age_ribosome <- age_gene_lr(info_by_cluster, expr_data_ribosome_scz)
lr_age_ubl <- age_gene_lr(info_by_cluster, expr_data_ubl)

# scatter plot of the regression coefficients
lr_coefficients_plot(lr_age_ribosome, "ribosome")
lr_coefficients_plot(lr_age_ubl, "UPS")
