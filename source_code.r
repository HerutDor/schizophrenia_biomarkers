# source code that contains the libarary and function to be used in the main
# code

# load packages ----
library(GEOquery) # download data from GEO
library(tidyverse) # data manipulation
library(matrixStats) # for rowSds
library(ggplot2) # for plotting
library(plotly) # for interactive plots
library(RColorBrewer) # colors to make heatmaps
library(pheatmap) # for heatmaps
library(dplyr) # for data manipulation
library(heatmaply) # interactive heatmaps
library(edgeR) # for differential expression
library(limma) # for differential expression
library(genefilter) # for filtering genes
library(gt) # for interactive tables
library(DT) # for interactive tables
library(factoextra) # for visualise clusters
library(fpc) # for visualise clusters
library(cluster) # for visualise clusters
library(NbClust) # for visualise clusters
library(cowplot) # multiple plot in the same figure
library(progress) # progress bar
library(readxl) # read excel files
library(grid) # for grid.arrange

# Functions ----

extract_geo_data <- function(geo_database) {
    # Download data from GEO
    gse <- getGEO(geo_database, GSEMatrix = TRUE)

    # Extract expression matrix
    expr_data <- exprs(gse[[1]]) %>%
        as.data.frame()

    # Extract annotations
    genes_info <- fData(gse[[1]])

    # create the sample info data frame
    sample_info <- pData(gse[[1]])

    # remove columns that are identical for all samples
    common_cols <- lengths(lapply(sample_info, unique))
    common_data <- sample_info[1, common_cols == 1]

    # Printing the common data about this dataset
    for (i in seq_along(common_data)) {
        cat(names(common_data)[i], ":", common_data[1, i], "\n")
    }

    # Return the data frames in a list
    return(list(expr_data = expr_data, genes_info = genes_info, sample_info = sample_info))
}

# Function to extract the gene symbol - this function is needed only for GSE18312
extract_gene_symbol <- function(text) {
    # Use regular expression to extract the string between the first two "//"
    gene_symbol <- sub(".*?// (.*?) //.*", "\\1", text)
    return(gene_symbol)
}

probe_to_gene <- function(expr_data, probes_id, genes_symbol) {
    # This function return the expr_data with additional column of genes symbol

    mm <- match(rownames(expr_data), probes_id)
    if (any(is.na(mm))) {
        # Keep only expression rows that has a matching gene
        rows_to_keep <- rownames(expr_data) %in% probes_id
        expr_data <- expr_data[rows_to_keep, ]
        mm <- match(rownames(expr_data), probes_id)
    }
    expr_data_gene <- cbind(expr_data, as.character(genes_symbol[mm]))
    colnames(expr_data_gene)[ncol(expr_data_gene)] <- "gene_name"

    return(expr_data_gene)
}

remove_duplicate_probes <- function(expr_data) {
    # expr_data: rownames are probe_id. the last column is gene_name

    # Handling duplicate probe using average expression since the variance dose not change much with the level of expression

    avg_expr <- rowMeans(expr_data[, -ncol(expr_data)], na.rm = TRUE)
    std_expr <- rowSds(as.matrix(expr_data[, -ncol(expr_data)], na.rm = TRUE))

    avg_expr_std <- as_tibble(cbind(avg_expr, std_expr))

    p <- ggplot(avg_expr_std, aes(avg_expr, std_expr)) +
        geom_point(size = 0.5) +
        geom_smooth(se = TRUE) +
        xlab(~mu) +
        ylab(~sigma) +
        ggtitle("The std as function of the avg expression of genes")
    print(p)

    # Print the number of unique and total probes
    cat("Total number of probes:", length(expr_data[, "gene_name"]), "\n")
    cat("Number of unique probes:", length(unique(expr_data[, "gene_name"])), "\n")

    # Ensure gene_name is the last column
    if (!"gene_name" %in% colnames(expr_data) || tail(colnames(expr_data), 1) != "gene_name") {
        stop("'gene_name' should be the last column of expr_data")
    }

    # Calculate the average for rows with the same gene_name
    averaged_data <- expr_data %>%
        group_by(gene_name) %>%
        summarise_all(mean) %>%
        as.data.frame()

    # Assign gene names as row names
    rownames(averaged_data) <- averaged_data$gene_name
    averaged_data$gene_name <- NULL # remove the gene_name column

    return(averaged_data)
}

filter_by_variability <- function(expr_data, sds_threshold) {
    # Filter out genes with low variability

    total_genes <- dim(expr_data)[1]

    # The overall variability of each gene
    sds <- apply(expr_data, 1, sd)
    sds0 <- sort(sds)
    # plot the distribution of the standard deviation
    plot(1:length(sds0), sds0, main = "Distribution of variability for all genes", sub = "Horizontal line the threshold standard deviation for filtering", xlab = "geneID index (from least to most variable)", ylab = "Standard deviation")
    abline(h = sds_threshold, col = "red", lty = "dashed")

    # remove gene with standard deviation less than the threshold
    genes_to_keep <- which(sds > sds_threshold)
    expr_data <- expr_data[genes_to_keep, ]

    # Print the number of genes removed
    genes_removed <- total_genes - length(genes_to_keep)
    cat("Total number of genes removed:", genes_removed, "\n")
    cat("Percantage of genes removed:", round(genes_removed / total_genes, 2), "\n")

    return(expr_data)
}

violin_plot <- function(expr_data, number_of_samples) {
    # expr_data: data frame with genes in rows and samples in columns.

    random_cols <- sample(ncol(expr_data), number_of_samples)
    random_data <- expr_data[, random_cols]

    random_data$genes_name <- rownames(expr_data)

    # Convert expression data to long format for visualization
    long_data <- random_data %>%
        pivot_longer(-genes_name, names_to = "Sample", values_to = "Expression")

    ggplot(long_data, aes(x = Sample, y = Expression, fill = Sample)) +
        geom_violin(trim = FALSE, show.legend = FALSE) +
        stat_summary(
            fun = "median",
            geom = "point",
            shape = 95,
            size = 10,
            color = "black",
            show.legend = FALSE
        ) +
        ggtitle("Violin plot of 10 Random Samples") +
        xlab("Sample Name") +
        ylab("Expression Level") +
        theme(plot.title = element_text(hjust = 0.5, size = 16))
}

histogram_expression <- function(expr_data) {
    # Histogram for the expression value of random sample
    random_col <- sample(1:ncol(expr_data), 1)
    hist(expr_data[, random_col], breaks = 100, main = "Histogram of gene Expression for random sample", xlab = "Expression Value")
}

heatmap_annotation <- function(expr_data, sample_info) {
    # display.brewer.all()
    expr_data_scale <- t(scale(t(expr_data), center = TRUE, scale = TRUE))

    # extract the "status" column
    sample_order <- order(sample_info$status)

    # order the columns of 'data' based on 'status' value
    ordered_data <- expr_data_scale[, sample_order]

    # Plot the heatmap with annotations
    pheatmap(ordered_data,
        # scale = "row", # a character indicating if the values should be centered and scaled.
        annotation_col = sample_info,
        main = paste("Heatmap dor data", geo_database),
        show_rownames = FALSE, show_colnames = FALSE,
        cluster_row = FALSE, cluster_cols = FALSE
    ) # This line ensures that columns are not re-clustered
}

find_differntially_expression_genes <- function(expr_data, sample_info, pVal_threshold, logFC_threshold) {
    # Set up your design matrix
    group <- sample_info$status %>%
        factor()
    design <- model.matrix(~ 0 + sample_info$status)
    colnames(design) <- levels(group)

    # fit a linear model to your data
    fit <- lmFit(expr_data, design)
    head(fit$coefficients)

    # extract the linear model fit
    contrasts <- makeContrasts(scz - control, levels = design)
    ebFit <- contrasts.fit(fit, contrasts)

    # apply the empirical Bayes’ step to get our differential expression statistics and p-values.
    ebFit <- eBayes(ebFit)

    # Extract a table of the top-ranked genes from a linear model fit. (set 'number=40000' to capture all genes)
    deg_stat <- topTable(ebFit, adjust = "BH", coef = 1, number = 40000, sort.by = "logFC") %>%
        as_tibble(rownames = "gene_name")

    # TopTable (from Limma) outputs a few different stats:
    # logFC, AveExpr, and P.Value should be self-explanatory
    # adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
    # B statistic is the log-odds that that gene is differentially expressed.
    # t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...because of Bayesian approach)

    # decideTests to pull out the DEGs
    significant_deg <- decideTests(ebFit, method = "global", adjust.method = "BH", p.value = pVal_threshold, lfc = logFC_threshold) # lcf - minimum absolute log2-fold-change required.

    return(list(deg_stat = deg_stat, significant_deg = significant_deg))
}

volcano_plot <- function(deg_stat, logFC_threshold, pVal_threshold) {
    # Initialize the color column
    deg_stat$color <- "Unchanged"

    # Define the color for each point based on logFC and p-value
    deg_stat$color[deg_stat$adj.P.Val < -log10(pVal_threshold) & deg_stat$logFC < -logFC_threshold] <- "Down-regulated"
    deg_stat$color[deg_stat$adj.P.Val < -log10(pVal_threshold) & deg_stat$logFC > logFC_threshold] <- "Up-regulated"

    vplot <- ggplot(deg_stat) +
        aes(x = logFC, y = -log10(adj.P.Val), text = paste("Symbol:", gene_name)) +
        geom_point(aes(color = color), size = 2 / 5) +
        xlab("log2(FC)") +
        ylab("-log10(FDR)") +
        scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
        guides(colour = guide_legend(override.aes = list(size = 1.5))) +
        ggtitle("Volcano plot: ration vs p-value")

    # Now make the volcano plot above interactive with plotly
    ggplotly(vplot)
}

filter_the_deg_genes <- function(significant_deg, deg_stat) {
    # take a look at what the results of decideTests looks like
    # head(significant_deg)
    print(summary(significant_deg))

    # retrieve the differentially expressed genes
    deg_genes <- deg_stat[as.vector(significant_deg[, 1]) != 0, ]

    # Create interactive tables to display your DEGs
    # gt(deg_stat[1:10,])
    t <- datatable(deg_genes,
        extensions = c("KeyTable", "FixedHeader"),
        caption = "Table 1: DEGs in schizoprenia",
        options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))
    ) %>%
        formatRound(columns = c(2:4, 7))
    print(t)

    return(deg_genes)
}

cluster_data <- function(data_to_cluster, cluster_method, number_of_clusters, show_plot) {
    # Allow input for cluster_method:
    # 1. kmeans
    # 2. pam. PAM is similar to k-means but more robust to outliers.
    # 3. hclust. Agglomerative hierarchical clustering

    if (cluster_method == "kmeans" || "pam") {
        # find the best k value by using the Elbow method. possible methods: "silhouette", "wss", "gap_stat"
        p <- fviz_nbclust(data_to_cluster, FUNcluster = kmeans, method = "wss") +
            geom_vline(xintercept = 2, linetype = 2)
        print(p)

        # Compute the clustering with k = number_of_clusters
        cluster_res <- eclust(data_to_cluster, cluster_method, k = number_of_clusters, nstart = 25, graph = FALSE)
    } else if (cluster_method == "hclust") {
        dist_method <- "euclidean"
        clust_method <- "ward.D" # Allowed values are those accepted by the function dist() [including “euclidean”, “manhattan”, “maximum”, “canberra”, “binary”, “minkowski”] and correlation based distance measures [“pearson”, “spearman” or “kendall”].

        # create the hierarchical tree
        cluster_res <- eclust(data_to_cluster, "hclust",
            k = 2,
            hc_metric = dist_method,
            clust_method = clust_method, graph = FALSE
        )

        # visualize hierarchical clustering as dendrogram
        cluster_dend <- fviz_dend(cluster_res,
            cex = 0.45, k = 2,
            main = paste("Dendrogram - ", clust_method),
            xlab = "Schizoprenia samples", ylab = "Distance", sub = "",
            rect = TRUE
        )
    } else {
        return("Not a known cluster method")
    }

    # visualize the clustering
    cluster_plot <- fviz_cluster(cluster_res, data_to_cluster,
        ellipse.type = "norm",
        geom = "point",
        main = paste(cluster_method, "clustering"),
        ggtheme = theme_minimal()
    )
    cluster_sil <- fviz_silhouette(cluster_res, ggtheme = theme_minimal())

    if (show_plot) {
        print(plot_grid(plotlist = list(cluster_plot, cluster_sil), ncol = 2, labels = c("A", "B")))
    }

    return(list(cluster_res = cluster_res, cluster_plot = cluster_plot, cluster_sil = cluster_sil))
}

find_significant_features <- function(mat) {
    # Check if input is matrix
    if (!is.matrix(mat)) {
        stop("Input should be a matrix!")
    }

    # Function to compute average silhouette width.
    # silhouette coefficient (Si) measures how similar an object i is to the the other objects in its own cluster versus those in the neighbor cluster. A value of Si close to 1 indicates that the object is well clustered.
    avg_silhouette <- function(data_to_cluster) {
        km.res <- eclust(data_to_cluster, "kmeans", k = 2, nstart = 25, graph = FALSE)
        silhouette_res <- silhouette(km.res$cluster, dist(data_to_cluster))
        return(mean(silhouette_res[, 3])) # return the average silhouette width
    }

    # Initialize data frame to store silhouette scores
    silhouette_scores <- data.frame(score = numeric(length(colnames(mat))))
    rownames(silhouette_scores) <- colnames(mat)

    num_iterations <- dim(mat)[2]
    pb <- progress_bar$new(total = num_iterations, format = "[:bar] :percent")

    # For each gene, remove it, cluster, and store silhouette score
    for (gene in colnames(mat)) {
        data_to_cluster <- mat[, !(colnames(mat) == gene)]
        silhouette_scores[gene, "score"] <- avg_silhouette(data_to_cluster)

        pb$tick()
    }

    return(silhouette_scores)
}

heatmap_cluster <- function(expr_data, cluster_per_sample, si_score_by_gene) {
    if (nrow(expr_data) == nrow(si_score_by_gene)) {
        expr_data <- t(expr_data)
    }
    # expr_data <- scale(expr_data, center = TRUE, scale = TRUE)
    col_name <- colnames(cluster_per_sample)[1]

    # extract the order of the samples by the cluster number
    sample_order <- order(cluster_per_sample[, col_name])
    genes_order <- order(si_score_by_gene$score)

    # convert the cluster number to string, so the annotation will presented better
    cluster_per_sample[, col_name] <- as.character(cluster_per_sample[, col_name])

    # order the columns of 'data' based on 'status' value
    ordered_data <- expr_data[sample_order, genes_order]

    # Plot the heatmap with annotations
    p <- pheatmap(ordered_data,
        # scale = "row", # a character indicating if the values should be centered and scaled.
        annotation_row = cluster_per_sample,
        annotation_col = si_score_by_gene,
        main = paste("Heatmap for data", geo_database),
        show_rownames = FALSE, show_colnames = FALSE,
        cluster_row = FALSE, cluster_cols = FALSE
    ) # This line ensures that columns are not re-clustered
    print(p)
}

# not used
heatmap_correlation_matrix <- function(data_to_cluster, cluster_per_sample, si_score_by_gene) {
    col_name <- colnames(cluster_per_sample)[1]
    # extract the order of the samples by the cluster number
    sample_order <- order(cluster_per_sample[, col_name])

    data_to_cluster <- cor(t(data_to_cluster))
    # order the columns of 'data' based on 'status' value
    ordered_data <- data_to_cluster[, sample_order]
    ordered_data <- ordered_data[sample_order, ]

    # Plot the heatmap with annotations
    p <- pheatmap(ordered_data,
        # scale = "row", # a character indicating if the values should be centered and scaled.
        annotation_col = cluster_per_sample,
        annotation_row = cluster_per_sample,
        main = paste("Heatmap for correlation matrix between samples in data", geo_database),
        show_rownames = FALSE, show_colnames = FALSE,
        cluster_row = FALSE, cluster_cols = FALSE
    ) # This line ensures that columns are not re-clustered
    print(p)
}

heatmap_clustering_genes <- function(expr_data, label_per_sample, gene_groups) {
    # Order the columns of 'expr_data' based on 'sample_order' order
    sample_order <- order(label_per_sample$label)
    ordered_data <- expr_data[, sample_order]

    # Convert the cluster number to string, so the annotation will be presented better
    label_per_sample$label <- as.character(label_per_sample$label)
    label_per_sample$label <- factor(label_per_sample$label)

    # Define the colors for the annotation
    annotation_colors_list <- c("#06D6A0", "#EF476F", "#FFD166")
    names(annotation_colors_list) <- c(0, 1, 2)

    # Create a custom annotation for the rows based on cluster assignment
    annotation_col <- data.frame(label = label_per_sample$label[sample_order])
    rownames(annotation_col) <- rownames(label_per_sample)[sample_order]
    annotation_col$label <- factor(annotation_col$label)

    # Order the rows of 'expr_data' based on 'gene_groups' order
    gene_order <- order(gene_groups$group)
    ordered_data <- ordered_data[gene_order, ]

    # Convert the cluster number to string, so the annotation will be presented better
    gene_groups$group <- factor(gene_groups$group)

    # Define the colors for the annotation
    annotation_colors_list_row <- c("#fbfbfb", "#118AB2", "#073B4C")
    names(annotation_colors_list_row) <- c("n", "r", "u")

    annotation_colors <- list(group = annotation_colors_list_row, label = annotation_colors_list)

    # Create a custom annotation based on cluster assignment and gene groups
    annotation_col <- data.frame(label = label_per_sample$label[sample_order])
    rownames(annotation_col) <- rownames(label_per_sample)[sample_order]
    annotation_col$label <- factor(annotation_col$label)

    annotation_row <- data.frame(group = gene_groups$group[gene_order])
    rownames(annotation_row) <- rownames(gene_groups)[gene_order]
    annotation_row$group <- factor(annotation_row$group)

    # Plot the heatmap with annotations
    p <- pheatmap(ordered_data,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        annotation_colors = annotation_colors,
        show_rownames = FALSE, show_colnames = FALSE,
        cluster_rows = FALSE, cluster_cols = FALSE,
        main = "Heatmap for the clustering genes",
        gaps_col = cumsum(c(96, 33, 73))
    )
    print(p)
}

pca_plot <- function(expr_data, label_per_sample, name) {
    # Perform PCA and keep the first three principal components
    pca_results <- prcomp(scale(t(expr_data)))
    pca_data <- as.data.frame(pca_results$x[, 1:3])

    # Add column for the cluster information (1 or 2), 0 indicates healthy control
    pca_data$label <- factor(label_per_sample$label,
        levels = c(0, 1, 2),
        labels = c("healthy controls", "cluster I", "cluster II")
    )

    # Define color mapping consistent with the previous plots
    color_mapping <- c("healthy controls" = "#06D6A0", "cluster I" = "#EF476F", "cluster II" = "#FFD166")

    # Plot the PCA with annotation of the clusters
    plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point(aes(color = label, fill = label), size = 6, shape = 21, alpha = 0.8) +
        scale_color_manual(values = color_mapping) +
        scale_fill_manual(values = color_mapping) +
        theme_minimal() +
        labs(
            x = "PC1",
            y = "PC2"
        ) +
        theme(
            panel.background = element_rect(fill = "white", color = "white"), # Set background to white
            plot.background = element_rect(fill = "white", color = "white"), # Set entire plot background to white
            legend.text = element_text(size = 12), # Size the legend text
            legend.title = element_blank(), # Remove the legend title
            legend.position = c(.87, .87), # Position the legend in the top right
            axis.title = element_text(size = 14), # Increase size of axis titles
            axis.text.x = element_text(size = 12), # Increase size of x-axis text
            axis.text.y = element_text(size = 12), # Increase size of y-axis text
            legend.background = element_rect(color = "black", size = 0.5) # Black border around the legend
        )

    # Save plot as high-resolution TIFF file
    ggsave(
        filename = "/Users/herutdor/Library/Mobile Documents/com~apple~CloudDocs/Herut/R/GSE38484/Paper/figure for submission/Fig_1.tiff",
        plot = plot,
        device = "tiff",
        width = 10, height = 10, units = "in",
        dpi = 300 # High resolution for publication
    )

    return(plot)
}

heatmap_gene_analysis <- function(expr_data, cluster_per_sample, gene_group_name) {
    # Prepare expression data
    samples_id <- colnames(expr_data)
    scz_samples <- intersect(samples_id, rownames(cluster_per_sample))
    control_samples <- setdiff(samples_id, scz_samples)

    expr_data_control <- expr_data[, colnames(expr_data) %in% control_samples]
    expr_data_scz <- expr_data[, colnames(expr_data) %in% scz_samples]

    mean_ribosome_control <- rowMeans(expr_data_control)
    expr_data_relative <- sweep(expr_data_scz, 1, mean_ribosome_control, `-`)

    # Extract the order of the samples by the cluster number
    sample_order <- order(cluster_per_sample$num_cluster)

    # Adjust the cluster numbers to labels "I" and "II"
    cluster_per_sample$num_cluster <- factor(cluster_per_sample$num_cluster, levels = c(1, 2), labels = c("I", "II"))

    # Order the columns of 'expr_data' based on 'sample_order'
    ordered_data <- expr_data_relative[, sample_order]

    # Define colors for the annotations
    annotation_colors <- list(cluster = c("I" = "#EF476F", "II" = "#FFD166"))

    # Create custom annotation data frame
    annotation_col <- data.frame(cluster = cluster_per_sample$num_cluster[sample_order])
    rownames(annotation_col) <- rownames(cluster_per_sample)[sample_order]

    genes_font_size <- 10
    if (nrow(expr_data) > 100) {
        genes_font_size <- 6
    }

    # Plot the heatmap
    p <- pheatmap(ordered_data,
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        fontsize_row = genes_font_size, # Set the font size for row names
        fontsize = 14,
        fontsize_number = 14,
        show_rownames = TRUE,
        show_colnames = FALSE,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = paste("Heatmap for the", gene_group_name, "genes"),
        annotation_legend_param = list(
            grid_height = unit(1, "cm"), # Increase the height of the legend
            grid_width = unit(1, "cm"), # Increase the width of the legend
            labels_gp = gpar(fontsize = 15) # Increase the font size of the legend labels
        ),
        legend = TRUE,
        gaps_col = cumsum(c(33))
    )

    print(p)
}

heatmap_gene_analysis_control <- function(expr_data, label_per_sample, gene_group_name) {
    # Prepare expression data
    control_samples <- rownames(label_per_sample)[label_per_sample$label == 0]
    expr_data_control <- expr_data[, colnames(expr_data) %in% control_samples]
    mean_control <- rowMeans(expr_data_control)
    expr_data_relative <- sweep(expr_data, 1, mean_control, `-`)

    # Adjust the cluster numbers to labels "I" and "II"
    label_per_sample$label <- factor(label_per_sample$label, levels = c(1, 2, 0), labels = c("I", "II", "Control  "))

    # Extract the order of the samples by the cluster number
    sample_order <- order(label_per_sample$label)

    # Order the columns of 'expr_data' based on 'sample_order'
    ordered_data <- expr_data_relative[, sample_order]

    # Define colors for the annotations
    annotation_colors <- list(group = c("Control  " = "#06D6A0", "I" = "#EF476F", "II" = "#FFD166"))

    # Create custom annotation data frame
    annotation_col <- data.frame(group = label_per_sample$label[sample_order])
    rownames(annotation_col) <- rownames(label_per_sample)[sample_order]

    genes_font_size <- 10
    if (nrow(expr_data) > 100) {
        genes_font_size <- 6
    }


    tiff(
        file = "/Users/herutdor/Library/Mobile Documents/com~apple~CloudDocs/Herut/R/GSE38484/Paper/figure for submission/Fig_2.tiff", # The directory you want to save the file in
        width = 10, # The width of the plot in inches
        height = 10, # The height of the plot in inches
        units = "in", # Units for width and height
        res = 300 # Resolution in DPI for high quality
    )

    # Plot the heatmap
    p <- pheatmap(ordered_data,
        show_rownames = TRUE,
        show_colnames = FALSE,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        fontsize_row = genes_font_size, # Set the font size for row names
        fontsize = 14,
        fontsize_number = 14,
        annotation_legend_param = list(
            grid_height = unit(1, "cm"), # Increase the height of the legend
            grid_width = unit(1, "cm"), # Increase the width of the legend
            labels_gp = gpar(fontsize = 15) # Increase the font size of the legend labels
        ),
        legend = TRUE,
        gaps_col = cumsum(c(33, 73, 96))
    )

    dev.off() # Close the graphics device
}

# Function to plot the confusion matrix data
bar_plot <- function(conf_matrix) {
    # Adding a column for fill to help with legend
    conf_matrix$Fill <- rep("All schizophrenia samples", nrow(conf_matrix))
    conf_matrix$Fill_TP <- "Schizophrenia samples the model predicted"

    # Create the bar plot
    plot <- ggplot() +
        geom_bar(
            data = conf_matrix, aes(x = Dataset, y = Total, fill = "All schizophrenia samples"),
            stat = "identity"
        ) +
        geom_bar(
            data = conf_matrix, aes(x = Dataset, y = TP, fill = "Schizophrenia samples the model predicted"),
            stat = "identity"
        ) +
        geom_text(
            data = conf_matrix,
            aes(x = Dataset, y = Total + 2, label = sprintf("%.0f%% (p-value = %.3f)", round(PPV * 100), round(p_value, 3))),
            size = 12, color = "#073B4C" # Adjust size for readability
        ) +
        scale_fill_manual(
            name = "Legend",
            values = c("All schizophrenia samples" = "#118AB2", "Schizophrenia samples the model predicted" = "#073B4C"),
            labels = c("All schizophrenia samples", "Schizophrenia samples the model predicted")
        ) +
        labs(
            x = "Dataset Name",
            y = "Number of Samples"
        ) +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "white", color = "white"), # Set panel background to white
            plot.background = element_rect(fill = "white", color = "white"), # Set entire plot background to white
            axis.title = element_text(size = 34), # Increase size of axis titles
            axis.text.x = element_text(size = 32), # Increase x-axis labels size (Dataset names)
            axis.text.y = element_text(size = 32), # Increase size of y-axis text
            legend.position = c(.73, .87), # Move legend to the top inside plot area
            legend.text = element_text(size = 34), # Size the legend text
            legend.title = element_blank(), # Remove the legend title
            legend.background = element_rect(color = "black", size = 2) # Black border around the legend
        )

    # Save the plot as a high-resolution TIFF file
    ggsave(
        filename = "/Users/herutdor/Library/Mobile Documents/com~apple~CloudDocs/Herut/R/GSE38484/Paper/figure for submission/Fig_3.tiff",
        plot = plot,
        device = "tiff",
        width = 30, height = 30, units = "in",
        dpi = 300
    ) # High resolution for publication

    return(plot)
}

# Statistical test for significantly different between the 2 clusters
wilcox_test_gene_cluster <- function(expr_data, label_per_sample, cluster1, cluster2) {
    # Initialize an empty data frame to store the results
    wilcoxon_results <- data.frame(W = numeric(), P_Value = numeric(), Direction = numeric(), stringsAsFactors = FALSE)

    # Loop over the genes
    for (gene_name in rownames(expr_data)) {
        # Extract the expression data for the current gene
        gene_expression <- expr_data[gene_name, ]

        # Split the gene expression data into two groups based on cluster assignments
        group1_expression <- as.numeric(unlist(gene_expression[label_per_sample$label == cluster1]))
        group2_expression <- as.numeric(unlist(gene_expression[label_per_sample$label == cluster2]))

        # Perform the Wilcoxon rank-sum test
        test_result <- wilcox.test(group1_expression, group2_expression)

        # Determine direction of the difference. Here, we use median, but mean could be used if appropriate for your data
        direction <- ifelse(median(group1_expression) > median(group2_expression), 1, 2)

        # Add the results to the data frame
        wilcoxon_results <- rbind(wilcoxon_results, data.frame(W = test_result$statistic, P_Value = test_result$p.value, Direction = direction))
    }

    # Set the row names of the results data frame to the gene names
    rownames(wilcoxon_results) <- rownames(expr_data)

    # correcte the p-value for multiple testing
    wilcoxon_results$P_Value <- p.adjust(wilcoxon_results$P_Value, method = "fdr")

    # print the results
    # Calculate the number of genes with p-value <= 0.05
    significant_genes <- wilcoxon_results[wilcoxon_results$P_Value <= 0.05, ]
    num_significant_genes <- nrow(significant_genes)

    # Calculate the percentage of significant genes
    percentage_significant <- (num_significant_genes / nrow(wilcoxon_results)) * 100

    # Determine the number of significant genes with higher expression in each cluster
    num_higher_in_cluster1 <- nrow(significant_genes[significant_genes$Direction == 1, ])
    num_higher_in_cluster2 <- nrow(significant_genes[significant_genes$Direction == 2, ])

    # Print the summary
    cat("Summary of Wilcoxon Test Results between", cluster1, "and", cluster2, ":\n")
    cat("Total number of significant genes (p-value <= 0.05):", num_significant_genes, "\n")
    cat("Percentage of significant genes:", sprintf("%.2f%%", percentage_significant), "\n")
    cat("Number of significant genes with higher expression in cluster", cluster1, ":", num_higher_in_cluster1, "\n")
    cat("Number of significant genes with higher expression in cluster", cluster2, ":", num_higher_in_cluster2, "\n")

    return(wilcoxon_results)
}

read_david_output <- function(genes_info) {
    # Read the specified range from the Excel file
    ribosomal_genes <- read_excel("cluster_analysis.xlsx", range = "F2:F12")
    ubl_genes <- read_excel("cluster_analysis.xlsx", range = "F15:F18")

    # Split each cell by commas and unlist into a single vector
    ribosomal_genes_numbers <-
        unlist(strsplit(ribosomal_genes$Genes, split = ",", fixed = TRUE))
    ubl_genes_numbers <-
        unlist(strsplit(ubl_genes$Genes, split = ",", fixed = TRUE))

    # Trim whitespace and convert to numeric
    ribosomal_genes_numbers <-
        as.numeric(trimws(ribosomal_genes_numbers)) %>% unique()
    ubl_genes_numbers <-
        as.numeric(trimws(ubl_genes_numbers)) %>% unique()

    # convert the genes numbers to genes ID
    ribosomal_genes_id <-
        genes_info[genes_info$Entrez_Gene_ID %in% ribosomal_genes_numbers, "ILMN_Gene"]
    ubl_genes_id <-
        genes_info[genes_info$Entrez_Gene_ID %in% ubl_genes_numbers, "ILMN_Gene"]

    return(list(ribosomal_genes_id = ribosomal_genes_id, ubl_genes_id = ubl_genes_id))
}

age_gene_correlation <- function(info_by_cluster, expr_data) {
    # Align the info data with the expression data
    aligned_ages <- as.numeric(t(info_by_cluster[colnames(expr_data), "age"]))

    # Initialize a vector to store the correlation results
    correlation_results <- matrix(nrow = nrow(expr_data), ncol = 1)
    colnames(correlation_results) <- "correlation"
    rownames(correlation_results) <- rownames(expr_data)

    # Loop through each gene and calculate the correlation with age
    for (gene in rownames(expr_data)) {
        gene_expression <- as.numeric(expr_data[gene, ])
        correlation_results[gene, "correlation"] <- cor(gene_expression, aligned_ages)
    }

    return(correlation_results)
}

age_gene_lr <- function(info_by_cluster, expr_data) {
    # Align the age data with the expression data
    aligned_info <- t(info_by_cluster[colnames(expr_data), ])

    # Initialize a vector to store the correlation results
    lr_results <- matrix(nrow = nrow(expr_data), ncol = 2)
    colnames(lr_results) <- c("age_coeffient", "cluster_coeffient")
    rownames(lr_results) <- rownames(expr_data)

    for (gene in rownames(expr_data)) {
        gene_expression <- as.numeric(expr_data[gene, ])
        age_numeric <- as.numeric(aligned_info["age", ]) / 10
        cluster_numeric <- as.numeric(aligned_info["num_cluster", ])

        lm_result <- lm(gene_expression ~ age_numeric + cluster_numeric)
        lr_results[gene, ] <- lm_result$coefficients[2:3]
    }

    return(lr_results)
}

# Function to plot the coefficients of the linear regression model
lr_coefficients_plot <- function(lr_age_ribosome, gene_group_name, text_subplot) {
    pdf(
        file = paste("/Users/herutdor/Library/Mobile Documents/com~apple~CloudDocs/Herut/R/GSE38484/Paper/figure for submission/Fig_", text_subplot, ".pdf", sep = ""),
        width = 10, # Width of the plot in inches
        height = 10, # Height of the plot in inches
        pointsize = 12 # Base font size
    )

    # scatter plot of the regression coefficients
    plot(lr_age_ribosome[, "age_coeffient"],
        ylim = range(-1, 1),
        xlab = "Genes", # X-axis label
        ylab = "Coefficients", # Y-axis label
        main = paste("Magnitude of the coefficients in linear regression model predicting level of expression - ", gene_group_name), # Main title
        col = "red", # Initial color for age coefficient
        pch = 15, # Square shape for points
        cex.main = 2.0, # Font size for the title, where 1.0 is the default size and 2.0 is twice as large
        cex.lab = 1.4, # Font size for x and y labels
        cex.axis = 1.4 # Font size for axis numbering
    )

    # Add text (S1a) in the top left corner
    text(x = 1, y = 0.95, labels = paste("(", text_subplot, ")", sep = ""), pos = 4, cex = 1.6, font = 2)

    # Add points for cluster coefficients
    points(lr_age_ribosome[, "cluster_coeffient"], col = "blue", pch = 17) # Triangle shape for points

    # Adding legend
    legend("topright", # Position of the legend
        legend = c("Age", "Cluster"), # Text in the legend
        col = c("red", "blue"), # Colors
        pch = c(15, 17), # Shapes
        title = "Coefficient Type", # Title of the legend
        cex = 1.4 # Font size for text within the legend
    )

    dev.off()
}
