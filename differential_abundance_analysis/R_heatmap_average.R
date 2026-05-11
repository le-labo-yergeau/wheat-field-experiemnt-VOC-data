library(dplyr)
library(pheatmap)


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
