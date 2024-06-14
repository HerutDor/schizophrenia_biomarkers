# import libraries and function from the source_code.R file
source("source_code.R")

# Extract data from GEO ----

geo_databases <- c("GSE38484", "GSE27383", "GSE38481", "GSE18312", "GSE48072")

for (geo_database in geo_databases) {
    data <- extract_geo_data(geo_database)
    expr_data <- data$expr_data
    genes_info <- data$genes_info
    sample_info <- data$sample_info

    path <- file.path(
        "/Users", "herutdor", "Documents", "Research", "Schizoprenia_biomarkers",
        paste0(geo_database, "_data"), geo_database
    )

    # keep the sample IDs, age and status of each sample

    # some of the dataset lack age information
    if (geo_database == "GSE27383") {
        sample_info <- sample_info %>%
            select("disease state:ch1")
        colnames(sample_info) <- c("status")
        sample_info["status"] <- ifelse(sample_info$status == "acutely admitted, severely psychotic schizophrenia patient", "scz", "control")
    } else if (geo_database == "GSE38481" || geo_database == "GSE38484") {
        sample_info <- sample_info %>%
            select("age:ch1", "gender:ch1", "status:ch1")
        colnames(sample_info) <- c("age", "gender", "status")
        sample_info["status"] <- ifelse(sample_info$status == "schizophrenia (SCZ)", "scz", "control")
    } else if (geo_database == "GSE18312") {
        sample_info <- sample_info %>%
            select("age:ch1", "gender:ch1", "diagnosis:ch1")
        colnames(sample_info) <- c("age", "gender", "status")
        # remove diagnosis with "bipolar disorder"
        sample_info <- sample_info[sample_info$status != "Bipolar Disorder", ]
        sample_info["status"] <- ifelse(sample_info$status == "Schizophrenia", "scz", "control")
    } else if (geo_database == "GSE48072") {
        sample_info <- sample_info %>%
            select("gender:ch1", "affection status:ch1")
        colnames(sample_info) <- c("gender", "status")
        sample_info["status"] <- ifelse(sample_info$status == "CASE", "scz", "control")
        # convert the geneder to lower case
        sample_info$gender <- tolower(sample_info$gender)
    } else {
        print("The dataset is not recognized.")
    }

    # convert the age to numeric
    if ("age" %in% colnames(sample_info)) {
        sample_info$age <- as.numeric(sample_info$age)
    }

    # print the demographic information
    print(table(sample_info$status))
    if ("age" %in% colnames(sample_info)) {
        print(round(mean(sample_info$age), 2))
        print(round(sd(sample_info$age), 2))
    }

    probes_id <- rownames(expr_data)
    if (geo_database == "GSE27383") {
        genes_symbol <- genes_info$"Gene Symbol"
    } else if (geo_database == "GSE18312") {
        # get the gene symbol from the gene_assignment acording to the structure ... // gene_symbol // ...
        genes_symbol <- sapply(genes_info$gene_assignment, extract_gene_symbol, USE.NAMES = FALSE)
    } else {
        genes_symbol <- genes_info$ILMN_Gene
    }

    expr_data <- probe_to_gene(expr_data, probes_id, genes_symbol)
    expr_data <- remove_duplicate_probes(expr_data)

    sds_threshold <- 0.05
    expr_data <- filter_by_variability(expr_data, sds_threshold)

    # visualize the data
    number_of_samples <- 10
    violin_plot(expr_data, number_of_samples)
    histogram_expression(expr_data)
    heatmap_annotation(expr_data, sample_info)

    # save the data
    saveRDS(expr_data, file = paste0(path, "_expr_data.rds"))
    saveRDS(genes_info, file = paste0(path, "_genes_info.rds"))
    saveRDS(sample_info, file = paste0(path, "_sample_info.rds"))
}
