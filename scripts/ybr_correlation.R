suppressMessages(library(tidyverse))
coexp_mat_file <- "/home/oma21/coexpression/data/interim/coexpression_matrix_rho_spqn.csv" #
obs_file <- "/home/oma21/coexpression/data/interim/pairwise_obs.csv" #



coexpression_matrix <- data.table::fread(coexp_mat_file)
obs <- data.table::fread(obs_file)

mat_matrix <- as.matrix(coexpression_matrix[, 2:ncol(coexpression_matrix)])
obs_matrix <- as.matrix(obs[, 2:ncol(coexpression_matrix)])

mat_matrix[obs_matrix < 400] <- NA

rownames(mat_matrix) <- colnames(mat_matrix)

con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
# translation <- DBI::dbGetQuery(con, "select * from aaron.scer_orfs_translation_info")
coexpression_orf_list <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")


annotated_orfs <- coexpression_orf_list %>%
    filter(is_canonical == "canonical") %>%
    dplyr::select(transcript, gene)



annotated_matrix <- mat_matrix[annotated_orfs$transcript[annotated_orfs$transcript %in% rownames(mat_matrix)], annotated_orfs$transcript[annotated_orfs$transcript %in% rownames(mat_matrix)]]


ybr <- "chr2_614024"
other_orfs <- rownames(annotated_matrix)[rownames(annotated_matrix) != ybr]
cor_data <- data.frame(matrix(ncol = 2, nrow = length(other_orfs)))
colnames(cor_data) <- c("orf_name", "correlation")
cor_data$orf_name <- other_orfs
for (orf_name in other_orfs) {
    corr_with_ybr <- cor(annotated_matrix[ybr, ], annotated_matrix[orf_name, ], use = "pairwise.complete.obs")

    cor_data[cor_data$orf_name == orf_name, "correlation"] <- corr_with_ybr
}
cor_data$rho <- annotated_matrix[ybr, colnames(annotated_matrix) %in% cor_data$orf_name]

ggplot(cor_data, aes(x = rho, y = correlation, color = orf_name == "chr5_138918")) +
    geom_point(alpha = .6) +
    theme_bw() +
    xlab("Coexpression") +
    ylab("Correlation of coexpression") +
    scale_color_discrete(name = "GCN4")
ggsave("foo.png", width = 7, height = 7)
