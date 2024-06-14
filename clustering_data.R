# import libraries and function from the source_code.R file
source("source_code.R")

# load processed data ----
geo_database <- "GSE38484"
path <- file.path(
    "/Users", "herutdor", "Documents", "Research", "Schizoprenia_biomarkers",
    paste0(geo_database, "_data"), geo_database
)

expr_data <- readRDS(file = paste0(path, "_expr_data.rds"))
sample_info <- readRDS(file = paste0(path, "_sample_info.rds"))
genes_info <- readRDS(file = paste0(path, "_genes_info.rds"))
si_threshold <- 0.33

# Clustering ----

### Find the differentially expressed genes ----
log_fc_threshold <- 0.1
p_val_threshold <- 0.05

deg_res <- find_differntially_expression_genes(
    expr_data, sample_info, p_val_threshold, log_fc_threshold
)
deg_stat <- deg_res$deg_stat
significant_deg <- deg_res$significant_deg
volcano_plot(deg_stat, log_fc_threshold, p_val_threshold)
deg_genes <- filter_the_deg_genes(significant_deg, deg_stat)

### Create the data for clustering ----
# The SCZ data
scz_data <- t(expr_data[, which(sample_info$status == "scz")])

# scz data that incluse only the DEG genes
diff_scz_data <- scz_data[, as.vector(significant_deg[, 1]) != 0]

### Cluster the data ----
set.seed(123)

### assessing clustering tendency by compute Hopkins statistic.
# result higher than 0.75 indicates a clustering tendency
# at the 90% confidence level.
hopkins_test <- get_clust_tendency(diff_scz_data,
    n = nrow(diff_scz_data) - 1, graph = FALSE
)
hopkins_test$hopkins_stat

# kmenas
cluster_data_res <- cluster_data(diff_scz_data, cluster_method = "kmeans", number_of_clusters = 2, show_plot = TRUE)
cluster_per_sample <- data.frame(cluster_data_res$cluster_res$cluster)
colnames(cluster_per_sample) <- "num_cluster"

control_samples <- rownames(sample_info[sample_info$status == "control", ])
scz_samples <- rownames(sample_info[sample_info$status == "scz", ])

label_per_sample <- data.frame(replicate(1, rep(0, nrow(sample_info))))
rownames(label_per_sample) <- rownames(sample_info)
colnames(label_per_sample) <- "label"
label_per_sample$label[rownames(label_per_sample) %in% control_samples] <- 0
label_per_sample$label[rownames(label_per_sample) %in%
    rownames(cluster_per_sample[cluster_per_sample$num_cluster == 1, , drop = FALSE])] <- 1
label_per_sample$label[rownames(label_per_sample) %in%
    rownames(cluster_per_sample[cluster_per_sample$num_cluster == 2, , drop = FALSE])] <- 2

### The most significant genes for the clustering ----
si_score_by_gene <- find_significant_features(diff_scz_data)
hist_info <- hist(si_score_by_gene$score, breaks = 20, col = "lightblue", border = "black")

expr_data_clustering <-
    expr_data[rownames(expr_data) %in% rownames(si_score_by_gene), ]

### Heatmaps ----
heatmap_cluster(diff_scz_data, cluster_per_sample, si_score_by_gene)

# Interpretation of the most important genes ----

# Find the score that 90% of the genes are above it
score_threshold <- quantile(si_score_by_gene$score, si_threshold)
clustering_genes <- rownames(si_score_by_gene[si_score_by_gene$score <= score_threshold, , drop = FALSE])

### save the genes lists for DAVID pathway analysis ----
write.table(clustering_genes, "clustering_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Because DAVID not identify the gene symbol, we'll convert the genes symbol to entrez_gene_ID
clustering_genes_entrez <- genes_info[genes_info$ILMN_Gene %in% clustering_genes, "Entrez_Gene_ID"]
write.table(clustering_genes_entrez, "clustering_genes_entrez.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Pathway analysis of the important genes to clustering ----
genes_id <- read_david_output(genes_info)
ribosomal_genes_id <- genes_id$ribosomal_genes_id
ubl_genes_id <- genes_id$ubl_genes_id

expr_data_ribosome <- expr_data[rownames(expr_data) %in% ribosomal_genes_id, ]
expr_data_ubl <- expr_data[rownames(expr_data) %in% ubl_genes_id, ]
# heatmap of the group genes expression level in the scz samples divided by the level of expression in the control samples
heatmap_gene_analysis(expr_data_ribosome, cluster_per_sample, "ribosome")
heatmap_gene_analysis(expr_data_ubl, cluster_per_sample, "ubiquitine")

# Statistical test for significantly different between the 2 clusters
clustering_wilcoxon_results <-
    wilcox_test_gene_cluster(expr_data_clustering, label_per_sample,
        cluster1 = "0", cluster2 = "2"
    )
ribosome_wilcoxon_results <-
    wilcox_test_gene_cluster(expr_data_ribosome, label_per_sample,
        cluster1 = "0", cluster2 = "1"
    )
ubl_wilcoxon_results <-
    wilcox_test_gene_cluster(expr_data_ubl, label_per_sample,
        cluster1 = "0", cluster2 = "1"
    )
ribosome_wilcoxon_results <-
    wilcox_test_gene_cluster(expr_data_ribosome, label_per_sample,
        cluster1 = "0", cluster2 = "2"
    )
ubl_wilcoxon_results <-
    wilcox_test_gene_cluster(expr_data_ubl, label_per_sample,
        cluster1 = "0", cluster2 = "2"
    )

deg_stat_ribosome <- deg_stat[deg_stat$gene_name %in% ribosomal_genes_id, ]
deg_stat_ubl <- deg_stat[deg_stat$gene_name %in% ubl_genes_id, ]

# Count the number of up-regulated genes & down-regulated genes in the ribosomal and ubl genes
up_regulated_ribosome <- deg_stat_ribosome[deg_stat_ribosome$logFC > 0 & deg_stat_ribosome$adj.P.Val < 0.05, ]
print(paste("The number of up-regulated genes in the ribosomal genes is", nrow(up_regulated_ribosome)))
down_regulated_ribosome <- deg_stat_ribosome[deg_stat_ribosome$logFC < 0 & deg_stat_ribosome$adj.P.Val < 0.05, ]
print(paste("The number of down-regulated genes in the ribosomal genes is", nrow(down_regulated_ribosome)))

up_regulated_ubl <- deg_stat_ubl[deg_stat_ubl$logFC > 0 & deg_stat_ubl$adj.P.Val < 0.05, ]
print(paste("The number of up-regulated genes in the ubl genes is", nrow(up_regulated_ubl)))
down_regulated_ubl <- deg_stat_ubl[deg_stat_ubl$logFC < 0 & deg_stat_ubl$adj.P.Val < 0.05, ]
print(paste("The number of down-regulated genes in the ubl genes is", nrow(down_regulated_ubl)))

# Count the number of up-regulated genes & down-regulated genes in the ribosomal and ubl genes for each cluster


# Save processed data ----
saveRDS(diff_scz_data, file = paste0(path, "_diff_scz_data.rds"))
saveRDS(cluster_per_sample, file = paste0(path, "_cluster_per_sample.rds"))
saveRDS(ribosomal_genes_id, file = paste0(path, "_ribosomal_genes_id.rds"))
saveRDS(ubl_genes_id, file = paste0(path, "_ubl_genes_id.rds"))
saveRDS(si_score_by_gene, file = paste0(path, "_si_score_by_gene.rds"))

# save the ribosome and UPS genes to txt file
write.table(ribosomal_genes_id, "ribosomal_genes_id.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(ubl_genes_id, "ubl_genes_id.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Statistical difference between the 2 clusters ----
# Merge the data frame
info_by_cluster <- merge(sample_info, cluster_per_sample, by = "row.names")
rownames(info_by_cluster) <- info_by_cluster$Row.names
info_by_cluster <- subset(info_by_cluster, select = -c(Row.names, status))

# Split the data into two groups based on the cluster number
cluster1_info <- info_by_cluster[info_by_cluster$num_cluster == 1, ]
cluster2_info <- info_by_cluster[info_by_cluster$num_cluster == 2, ]

### Age comparison ----

expr_data_ribosome_scz <- expr_data_ribosome[, scz_samples]

# print the avergae and the standard deviation of the age in each cluster
print(sprintf(
    "The average age in cluster 1 is %.2f with a standard deviation of %.2f",
    mean(cluster1_info$age), sd(cluster1_info$age)
))
print(sprintf(
    "The average age in cluster 2 is %.2f with a standard deviation of %.2f",
    mean(cluster2_info$age), sd(cluster2_info$age)
))

# use a t-test because the distribution are normal (checked using Shapiro-Wilk test)
t_test_age <- t.test(cluster1_info$age, cluster2_info$age)
print(t_test_age) # The result's ahow statistically significant difference between the 2 clusters so we cheack if this diffrence is also clinically significant

# correlation between age and the each gene expression
cor_age_ribosome <- age_gene_correlation(info_by_cluster, expr_data_ribosome_scz)
# scatter plot of the correlation
plot(cor_age_ribosome, main = "Correlation between age and ribosome genes expression", ylab = "Correlation")

# linear regession with the level of expression as the depandent variable
# and the age and cluster as the indepandant variables
threshold <- 0.1

lr_age_ribosome <- age_gene_lr(info_by_cluster, expr_data_ribosome_scz)
# scatter plot of the regression coefficents
plot(lr_age_ribosome[, "age_coeffient"],
    ylim = range(-1, 1),
    main = "Linear regression between age and ribosome genes expression", ylab = "Coefficients", col = "red"
) # age
points(lr_age_ribosome[, "cluster_coeffient"], col = "blue") # cluster

# Evaluate the size of the coefficients
small_age_coefficients <- abs(lr_age_ribosome[, "age_coeffient"]) < threshold
print(paste(
    "The percantage of genes with small age coefficients is",
    sum(small_age_coefficients) / length(small_age_coefficients) * 100, "%"
))
small_cluster_coefficients <- abs(lr_age_ribosome[, "cluster_coeffient"]) < threshold
print(paste(
    "The percantage of genes with small cluster coefficients is",
    sum(small_cluster_coefficients) / length(small_cluster_coefficients) * 100, "%"
))

lr_age_ubl <- age_gene_lr(info_by_cluster, expr_data_ubl)
# scatter plot of the regression coefficents
plot(lr_age_ubl[, "age_coeffient"],
    ylim = range(-1, 1),
    main = "Linear regression between age and ribosome genes expression", ylab = "Coefficients", col = "red"
) # age
points(lr_age_ubl[, "cluster_coeffient"], col = "blue") # cluster

# Evaluate the size of the coefficients
small_age_coefficients <- abs(lr_age_ubl[, "age_coeffient"]) < threshold
print(paste(
    "The percantage of genes with small age coefficients is",
    sum(small_age_coefficients) / length(small_age_coefficients) * 100, "%"
))
small_cluster_coefficients <- abs(lr_age_ubl[, "cluster_coeffient"]) < threshold
print(paste(
    "The percantage of genes with small cluster coefficients is",
    sum(small_cluster_coefficients) / length(small_cluster_coefficients) * 100, "%"
))

### Gender comparison ----
# Use a chi-squared test or Fisher's exact test depending on the size of the data

# Table of counts for chi-squared test
gender_table <- table(info_by_cluster$num_cluster, info_by_cluster$gender)

# If the expected counts are too low for chi-squared test, use Fisher's exact test
fisher_test_gender <- fisher.test(gender_table)

# Print the results
print(fisher_test_gender)

### Anti-psychotic drugs ----

# AP_change <- read.table(file=paste0(path, "gene_expression_change.txt"), header = TRUE)
AP_change <- read.table(file = paste0(path, "gene_symbols_worst_responser.txt"), header = TRUE)
AP_change <- read.table(file = paste0(path, "gene_symbols_best_responser.txt"), header = TRUE)
common_genes_si_AP <- intersect(clustering_genes, AP_change$Gene_symbol)

score_threshold <- quantile(si_score_by_gene$score, si_threshold)
clustering_genes <- rownames(si_score_by_gene[si_score_by_gene$score <= score_threshold, , drop = FALSE])
