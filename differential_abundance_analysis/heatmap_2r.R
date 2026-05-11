library(dplyr)
library(pheatmap)

# =========================
# 1. FILTER FEATURES
# =========================

master_df_filtered <- master_df %>%
  filter(log2FoldChange >= 3)


# =========================
# 2. FILTER HEATMAP MATRIX
# =========================

common_features <- intersect(
  rownames(vsd),
  master_df_filtered$Feature
)

heatmap_matrix <- assay(vsd)[common_features, ]


# =========================
# 3. ORDER SAMPLES
#    NO COLUMN CLUSTERING
# =========================

desired_order <- c(
  "DS25",
  "DS50",
  "DS75",
  "DS100",
  "DT25",
  "DT50",
  "DT75",
  "DT100"
)

sample_treatment <- gsub(
  "(\\d+)_(DS|DT)_.*",
  "\\2\\1",
  colnames(heatmap_matrix)
)

sample_order_df <- data.frame(
  Sample = colnames(heatmap_matrix),
  Treatment = sample_treatment
)

sample_order_df$Treatment <- factor(
  sample_order_df$Treatment,
  levels = desired_order
)

sample_order_df <- sample_order_df[
  order(sample_order_df$Treatment),
]

ordered_samples <- sample_order_df$Sample

# reorder matrix
heatmap_matrix <- heatmap_matrix[, ordered_samples]


# =========================
# 4. COLUMN ANNOTATION
# =========================

annotation_col <- data.frame(
  Treatment = sample_order_df$Treatment
)

rownames(annotation_col) <- ordered_samples


# =========================
# 5. TREATMENT COLORS
# =========================

treatment_colors <- c(
  "DS25"  = "#006d2c",
  "DS50"  = "#31a354",
  "DS75"  = "#74c476",
  "DS100" = "#c7e9c0",
  "DT25"  = "#8c2d04",
  "DT50"  = "#d95f0e",
  "DT75"  = "#fe9929",
  "DT100" = "#feedde"
)

ann_colors <- list(
  Treatment = treatment_colors
)


# =========================
# 6. ROW ANNOTATION (LFC)
# =========================

row_annotation <- master_df_filtered %>%
  select(Feature, log2FoldChange)

# make rownames unique
rownames(row_annotation) <- make.unique(
  as.character(row_annotation$Feature)
)

row_annotation$Feature <- NULL

# match heatmap rownames
rownames(heatmap_matrix) <- make.unique(
  rownames(heatmap_matrix)
)

row_annotation <- row_annotation[
  rownames(heatmap_matrix),
  ,
  drop = FALSE
]


# =========================
# 7. COLOR SCALE OPTIMIZATION
# =========================

scaled_matrix <- t(scale(t(heatmap_matrix)))

range(scaled_matrix, na.rm = TRUE)
summary(as.vector(scaled_matrix))

# balanced visualization
break_limit <- 1.5

breaks <- seq(
  -break_limit,
  break_limit,
  length.out = 101
)

colors <- colorRampPalette(
  c("blue", "yellow", "red")
)(100)


# =========================
# 8. HEATMAP
# =========================

pheatmap(
  heatmap_matrix,
  
  scale = "row",
  
  # IMPORTANT:
  # keep treatment order fixed
  cluster_cols = FALSE,
  
  # cluster metabolites/features
  cluster_rows = TRUE,
  
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  
  annotation_colors = ann_colors,
  
  # remove labels beside annotations
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  
  show_rownames = FALSE,
  show_colnames = FALSE,
  
  border_color = NA,
  
  # heatmap colors
  color = colors,
  breaks = breaks
)


# =========================
# 9. SAVE HEATMAP MATRIX
# =========================

write.csv(
  as.data.frame(heatmap_matrix),
  file = "heatmap_matrix.csv",
  row.names = TRUE
)


# =========================
# 10. MERGE MATRIX + METADATA
# =========================

heatmap_df <- as.data.frame(heatmap_matrix)

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
# 11. SAVE MERGED TABLE
# =========================

write.csv(
  merged_heatmap,
  file = "heatmap_matrix_with_metadata.csv",
  row.names = FALSE
)

