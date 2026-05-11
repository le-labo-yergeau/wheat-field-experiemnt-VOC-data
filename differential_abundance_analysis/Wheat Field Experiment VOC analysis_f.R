# =========================
# LOAD LIBRARIES
# =========================
library(dplyr)
library(readxl)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(forcats)
library(stringr)
library(patchwork)
library(officer)
library(rvg)
library(pheatmap)
# =========================
# LOAD DATA
# =========================
counts_data <- read_excel(
  "Field_exp_wheat_alignment0.5_sorted_nist_metabo.xlsx", 
  sheet = "Field_exp_wheat_alignment0.5_so"
)

metadata <- read_excel("meta.xlsx", sheet = "meta")

# =========================
# PREPARE COUNT MATRIX
# =========================
# Move feature column to rownames
data_org_n <- counts_data %>% column_to_rownames(var = colnames(data_org_n)[1])

# Convert to numeric
data_org_n <- data_org_n %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

count_table <- data_org_n
# =========================
# FILTER FEATURES
# present in ≥3 replicates in ≥1 group
# =========================
group_names <- sub("_[0-9]+$", "", colnames(count_table))
unique_groups <- unique(group_names)

keep_features <- apply(count_table, 1, function(x) {
  any(sapply(unique_groups, function(g) {
    idx <- which(group_names == g)
    sum(x[idx] > 0, na.rm = TRUE) >= 3
  }))
})

filtered_table <- count_table[keep_features, ]

# =========================
# HALF-MIN IMPUTATION
# =========================
data_imputed <- t(apply(filtered_table, 1, function(x) {
  min_val <- min(x[x > 0], na.rm = TRUE)
  
  if (is.infinite(min_val)) return(x)
  
  x[is.na(x)] <- min_val / 2
  return(x)
}))

data_imputed <- as.data.frame(data_imputed)
colnames(data_imputed) <- colnames(filtered_table)
rownames(data_imputed) <- rownames(filtered_table)

# =========================
# PREPARE METADATA
# =========================
# Keep only matching samples
meta <- meta %>% filter(ID %in% colnames(data_imputed))

# Reorder metadata to match count matrix
meta <- meta[match(colnames(data_imputed), meta$ID), ]

# Check alignment
stopifnot(all(meta$ID == colnames(data_imputed)))

# Convert to factors
meta$Treatment   <- factor(meta$Treatment, levels = c("100","75","50","25"))
meta$Compartment <- factor(meta$Compartment, levels = c("DS","DT"))  # DS = reference
meta$Block       <- factor(meta$Block)

# =========================
# CREATE DESeq2 OBJECT
# =========================
dds <- DESeqDataSetFromMatrix(
  countData = round(data_imputed),
  colData   = meta,
  design    = ~ Block + Treatment + Compartment + Treatment:Compartment
)

# =========================
# RUN DESeq2
# =========================
dds <- DESeq(dds)

# =========================
# SAVE IMPUTED TABLE
# =========================
write.csv(
  data_imputed,
  "01_imputed_table.csv",
  row.names = TRUE
)

# =========================
# DESeq2 NORMALIZATION
# =========================

# VST (PCA / heatmap)
vsd <- vst(dds, blind = FALSE)
write.csv(
  assay(vsd),
  "02_deseq2_vst_normalized.csv",
  row.names = TRUE
)

# Normalized counts (raw scaling)
write.csv(
  counts(dds, normalized = TRUE),
  "03_deseq2_normalized_counts.csv",
  row.names = TRUE
)


# =========================
# PCA (VST data)
# =========================
pca_data <- plotPCA(vsd, intgroup = c("Treatment", "Compartment"), returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Treatment, shape = Compartment)) +
  geom_point(size = 4, alpha = 0.9) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

print(pca_plot)


# =========================
# HEATMAP (Top 20 features)
# =========================
select <- order(rowVars(assay(vsd)), decreasing = TRUE)[1:50]

top_features <- rownames(vsd)[select]

heatmap_matrix <- assay(vsd)[select, ]

pheatmap(
  heatmap_matrix,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = FALSE
)

# Optional: save feature list
write.csv(data.frame(Feature = top_features),
          "top20_heatmap_features.csv",
          row.names = FALSE)
# =========================
# MA PLOT
# =========================
res <- results(dds)

# Base MA plot
plotMA(res, ylim = c(-5, 5))

# Optional: ggplot MA plot (clean publication version)
ma_data <- as.data.frame(res) %>%
  mutate(significance = ifelse(pvalue < 0.05, "Significant", "NS"))

ma_plot <- ggplot(ma_data, aes(x = baseMean, y = log2FoldChange, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_x_log10() +
  scale_color_manual(values = c("Significant" = "red", "NS" = "grey")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(x = "Mean Expression (log scale)", y = "Log2 Fold Change")

print(ma_plot)

# =========================
# SUMMARY FUNCTION (pvalue )
# =========================
get_sig_summary <- function(res, name) {
  data.frame(
    Comparison = name,
    Increased = sum(res$log2FoldChange > 0 & res$pvalue <= 0.05, na.rm = TRUE),
    Decreased = sum(res$log2FoldChange < 0 & res$pvalue <= 0.05, na.rm = TRUE),
    Non_significant = sum(res$pvalue > 0.05 | is.na(res$pvalue), na.rm = TRUE)
  )
}


# =========================
# CONTRASTS
# =========================

# DS
res_25_DS <- results(dds, contrast = c("Treatment", "25", "100"))
res_50_DS <- results(dds, contrast = c("Treatment", "50", "100"))
res_75_DS <- results(dds, contrast = c("Treatment", "75", "100"))

# DT
res_25_DT <- results(dds, contrast = list(c("Treatment_25_vs_100", "Treatment25.CompartmentDT")))
res_50_DT <- results(dds, contrast = list(c("Treatment_50_vs_100", "Treatment50.CompartmentDT")))
res_75_DT <- results(dds, contrast = list(c("Treatment_75_vs_100", "Treatment75.CompartmentDT")))

keep_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue")

# =========================
# ORDER RESULTS
# =========================
cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue")
res_25_DS <- as.data.frame(res_25_DS)[order(res_25_DS$pvalue), cols]
res_50_DS <- as.data.frame(res_50_DS)[order(res_50_DS$pvalue), cols]
res_75_DS <- as.data.frame(res_75_DS)[order(res_75_DS$pvalue), cols]

res_25_DT <- as.data.frame(res_25_DT)[order(res_25_DT$pvalue), cols]
res_50_DT <- as.data.frame(res_50_DT)[order(res_50_DT$pvalue), cols]
res_75_DT <- as.data.frame(res_75_DT)[order(res_75_DT$pvalue), cols]

# =========================
# SAVE ALL RESULTS
# =========================
write.csv(as.data.frame(res_25_DS), "DS_25_vs_100.csv")
write.csv(as.data.frame(res_50_DS), "DS_50_vs_100.csv")
write.csv(as.data.frame(res_75_DS), "DS_75_vs_100.csv")

write.csv(as.data.frame(res_25_DT), "DT_25_vs_100.csv")
write.csv(as.data.frame(res_50_DT), "DT_50_vs_100.csv")
write.csv(as.data.frame(res_75_DT), "DT_75_vs_100.csv")

# =========================
# SUMMARY TABLE
# =========================
summary_table <- do.call(rbind, list(
  get_sig_summary(res_25_DS, "DS: 25 vs 100"),
  get_sig_summary(res_50_DS, "DS: 50 vs 100"),
  get_sig_summary(res_75_DS, "DS: 75 vs 100"),
  get_sig_summary(res_25_DT, "DT: 25 vs 100"),
  get_sig_summary(res_50_DT, "DT: 50 vs 100"),
  get_sig_summary(res_75_DT, "DT: 75 vs 100")
))

# Save summary
write.csv(summary_table, "summary_table.csv", row.names = FALSE)

# Print
print(summary_table)

# =========================
# CREATE MASTER LONG FILE (pvalue <= 0.05)
# =========================

res_list <- list(
  "DS_25_vs_100" = res_25_DS,
  "DS_50_vs_100" = res_50_DS,
  "DS_75_vs_100" = res_75_DS,
  "DT_25_vs_100" = res_25_DT,
  "DT_50_vs_100" = res_50_DT,
  "DT_75_vs_100" = res_75_DT
)

master_list <- lapply(names(res_list), function(name) {
  
  df <- as.data.frame(res_list[[name]])
  
  # keep feature names
  df$Feature <- rownames(df)
  
  # filter pvalue
  df <- df[df$pvalue <= 0.05, ]
  
  # create Treatment column (DT25, DS50, etc.)
  treatment_label <- sub("_vs_.*", "", name)
  treatment_label <- gsub("_", "", treatment_label)
  df$Treatment <- treatment_label
  
  # add comparison column
  df$Comparison <- name
  
  # split Feature into mz_rt and compound_info
  df <- df %>%
    mutate(
      mz_rt = str_extract(Feature, "^[^\\(]+") %>% str_trim(),
      compound_info = str_extract(Feature, "\\(([^\\)]+)\\)") %>%
        str_replace_all("[\\(\\)]", "") %>%
        str_trim()
    )
  
  # reorder columns
  df <- df[, c("Treatment", "Comparison", "mz_rt", "compound_info",
               setdiff(colnames(df), c("Treatment", "Comparison", "mz_rt", "compound_info")))]
  
  return(df)
})

# combine all
master_df <- do.call(rbind, master_list)

# save
write.csv(master_df, "Master_long_significant_pvalue_split.csv", row.names = FALSE)

# check
head(master_df)
dim(master_df)


# -------------------------
#  FIGURE 7, Import Dot Plot Data (C)
# -------------------------
Figure7_ABC_Data_C <- read_excel("Figure7_ABC_Data.xlsx", 
                                    sheet = "C")

plot_data <- Figure7_ABC_Data %>%
  filter(Status == "Production") %>%
  mutate(Facet = case_when(
    str_starts(Treatment, "DS") ~ "DS",
    str_starts(Treatment, "DT") ~ "DT",
    TRUE ~ "Other"
  )) %>%
  group_by(compound_info) %>%
  mutate(total_count = sum(Count)) %>%
  ungroup() %>%
  mutate(compound_info = fct_reorder(compound_info, total_count))

# Fix dot plot color
status_colors <- c("Production" = "#DC143C")  # keep single color

plot_dot <- ggplot(plot_data, aes(x = Treatment, y = compound_info, size = Count, color = Status)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = status_colors, guide = "none") +  # keep color
  scale_size_continuous(range = c(3, 12), name = "Count") +
  facet_wrap(~Facet, nrow = 1, scales = "free_x") +
  labs(x = "Treatment", y = "") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1, "lines")
  )

# -------------------------
# Import Bar Plot (A)
# -------------------------
Figure7_ABC_Data_A <- read_excel("Figure7_ABC_Data.xlsx", 
                               sheet = "A") %>%
  mutate(
    Status = str_trim(Status),
    Status = case_when(
      Status %in% c("Emission increased", "Production") ~ "Production",
      Status %in% c("Emission decreased", "Consumption") ~ "Consumption",
      TRUE ~ Status
    ),
    Status = factor(Status, levels = c("Production", "Consumption"))
  )

plot_bar <- ggplot(Figure7_ABC_Data_A, aes(x = Treatment, y = VOCs, fill = Status)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  geom_text(aes(label = VOCs),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("Production" = "#DC143C", "Consumption" = "#1E90FF")) +
  labs(x = "Treatment", y = "Number of VOCs", fill = "Status") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

# -------------------------
# Import Line Graph (B)
# -------------------------
Figure7_ABC_Data_B <- read_excel("Figure7_ABC_Data.xlsx", 
                                 sheet = "B") %>%
  mutate(Treatment_num = str_extract(Treatment, "\\d+"),
         Treatment_num = factor(Treatment_num, levels = c("25","50","75")))

plot_line <- ggplot(Figure7_ABC_Data_B, aes(x = Treatment_num, y = Percent_VOCs, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = c("DS" = "#006400", "DT" = "#FF8C00")) +
  labs(x = "Treatment", y = "% DA VOCs Produced", color = "Genotype") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "top"
  )

# -------------------------
#Save SVG files
# -------------------------
ggsave("Figure7A_Bar.svg", plot_bar, width = 7, height = 5)
ggsave("Figure7B_Line.svg", plot_line, width = 7, height = 5)
ggsave("Figure7C_Dot.svg", plot_dot, width = 10, height = 7)

# Create heatmap for features with log2 fold change ≥ 3
# =========================
# FILTER FEATURES
# =========================

master_df_filtered <- master_df %>%
  filter(log2FoldChange >= 3)


# =========================
# FILTER HEATMAP MATRIX
# =========================

common_features <- intersect(
  rownames(vsd),
  master_df_filtered$Feature
)

heatmap_matrix <- assay(vsd)[common_features, ]


# =========================
# CREATE GROUP MEANS
#    (average of 6 replicates)
# =========================

sample_groups <- gsub(
  "(\\d+)_(DS|DT)_.*",
  "\\1_\\2",
  colnames(heatmap_matrix)
)

unique_groups <- unique(sample_groups)

heatmap_mean <- sapply(unique_groups, function(g) {
  rowMeans(
    heatmap_matrix[, sample_groups == g, drop = FALSE],
    na.rm = TRUE
  )
})

heatmap_mean <- as.matrix(heatmap_mean)

# desired order
desired_order <- c(
  "25_DS",
  "50_DS",
  "75_DS",
  "100_DS",
  "25_DT",
  "50_DT",
  "75_DT",
  "100_DT"
)

heatmap_mean <- heatmap_mean[, desired_order]


# =========================
# COLUMN ANNOTATION
# =========================

annotation_col <- data.frame(
  Treatment = desired_order
)

rownames(annotation_col) <- desired_order


# =========================
# TREATMENT COLORS
# =========================

treatment_colors <- c(
  "25_DS"  = "#006d2c",
  "50_DS"  = "#31a354",
  "75_DS"  = "#74c476",
  "100_DS" = "#c7e9c0",
  "25_DT"  = "#8c2d04",
  "50_DT"  = "#d95f0e",
  "75_DT"  = "#fe9929",
  "100_DT" = "#feedde"
)

ann_colors <- list(
  Treatment = treatment_colors
)


# =========================
# ROW ANNOTATION (LFC)
# =========================

row_annotation <- master_df_filtered %>%
  select(Feature, log2FoldChange)

rownames(row_annotation) <- make.unique(
  as.character(row_annotation$Feature)
)

row_annotation$Feature <- NULL

rownames(heatmap_mean) <- make.unique(
  rownames(heatmap_mean)
)

row_annotation <- row_annotation[
  rownames(heatmap_mean),
  ,
  drop = FALSE
]


# =========================
# LFC COLORS
# =========================

lfc_colors <- colorRampPalette(
  c("white", "red")
)(100)

ann_colors$log2FoldChange <- lfc_colors


# =========================
# HEATMAP COLORS
# =========================

scaled_matrix <- t(scale(t(heatmap_mean)))

break_limit <- 2.5

breaks <- seq(
  -break_limit,
  break_limit,
  length.out = 101
)

colors <- colorRampPalette(
  c("blue", "yellow", "red")
)(100)


# =========================
# ORDER FEATURES BY LFC
# =========================

feature_order <- master_df_filtered %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "DS25",
        "DS50",
        "DS75",
        "DS100",
        "DT25",
        "DT50",
        "DT75",
        "DT100"
      )
    )
  ) %>%
  arrange(
    Treatment,
    desc(log2FoldChange)
  ) %>%
  pull(Feature)

# keep only features present in matrix
feature_order <- intersect(
  feature_order,
  rownames(heatmap_mean)
)

# reorder heatmap rows
heatmap_mean <- heatmap_mean[
  feature_order,
]

# =========================
# 9. HEATMAP
# =========================

Fig_7n <-pheatmap(
  heatmap_mean,
  
  scale = "row",
  
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  
  annotation_colors = ann_colors,
  
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  
  show_colnames = FALSE,
  show_rownames = FALSE,
  
  # split between DS and DT
  gaps_col = 4,
  
  border_color = NA,
  
  color = colors,
  breaks = breaks
)


# =========================
# 10. SAVE HEATMAP MATRIX
# =========================

write.csv(
  as.data.frame(heatmap_mean),
  file = "heatmap_mean_matrix.csv",
  row.names = TRUE
)


# =========================
# MERGE MATRIX + METADATA
# =========================

heatmap_df <- as.data.frame(heatmap_mean)

heatmap_df$Feature <- rownames(heatmap_df)

merged_heatmap <- left_join(
  heatmap_df,
  master_df_filtered %>%
    select(
      Feature,
      Treatment,
      pvalue,
      log2FoldChange
    ),
  by = "Feature"
)

merged_heatmap <- merged_heatmap %>%
  relocate(
    Feature,
    Treatment,
    log2FoldChange,
    pvalue
  )


# =========================
# SAVE MERGED TABLE
# =========================

write.csv(
  merged_heatmap,
  file = "heatmap_mean_matrix_with_metadata.csv",
  row.names = FALSE
)


# =========================
# SAVE HEATMAP AS SVG
# =========================

svg(
  "heatmap_Fig_7n.svg",
  width = 8,
  height = 8
)

grid::grid.newpage()
grid::grid.draw(Fig_7n$gtable)

dev.off()
