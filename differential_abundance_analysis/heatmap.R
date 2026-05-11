library(dplyr)
library(pheatmap)

# =========================
# 1. FILTER FEATURES
# =========================

master_df_filtered <- master_df %>%
  filter(log2FoldChange >= 3)

# =========================
# 2. FILTER MATRIX
# =========================

common_features <- intersect(
  rownames(vsd),
  master_df_filtered$Feature
)

heatmap_matrix <- assay(vsd)[common_features, ]


# =========================
# 3. COLUMN ANNOTATION (ONLY TREATMENT)
# =========================

annotation_col <- data.frame(
  Treatment = gsub("(\\d+)_(DS|DT)_.*", "\\2\\1", colnames(heatmap_matrix))
)

rownames(annotation_col) <- colnames(heatmap_matrix)

annotation_col$Treatment <- factor(
  annotation_col$Treatment,
  levels = c(
    "DS25","DS50","DS75","DS100",
    "DT25","DT50","DT75","DT100"
  )
)


# =========================
# 4. TREATMENT COLORS
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
# 5. ROW ANNOTATION (LFC)
# =========================

row_annotation <- master_df_filtered %>%
  select(Feature, log2FoldChange)

rownames(row_annotation) <- make.unique(row_annotation$Feature)
row_annotation$Feature <- NULL

row_annotation <- row_annotation[
  rownames(heatmap_matrix),
  ,
  drop = FALSE
]


# =========================
# 6. LFC RED GRADIENT (IMPORTANT FIX)
# =========================

lfc_breaks <- seq(
  min(row_annotation$log2FoldChange, na.rm = TRUE),
  max(row_annotation$log2FoldChange, na.rm = TRUE),
  length.out = 101
)

lfc_colors <- colorRampPalette(c("white", "red"))(100)

ann_colors$log2FoldChange <- lfc_colors


# =========================
# 7. HEATMAP
# =========================

pheatmap(
  heatmap_matrix,
  
  scale = "row",
  
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  
  annotation_colors = ann_colors,
  
  # THIS REMOVES "LFC" TEXT FROM HEATMAP SIDE
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  
  show_rownames = FALSE,
  show_colnames = FALSE,
  
  border_color = NA
)

# =========================
# 8. SAVE HEATMAP MATRIX
# =========================

write.csv(
  as.data.frame(heatmap_matrix),
  file = "heatmap_matrix.csv",
  row.names = TRUE
)


# =========================
# MERGE HEATMAP MATRIX + FULL METADATA
# =========================

heatmap_df <- as.data.frame(heatmap_matrix)
heatmap_df$Feature <- rownames(heatmap_df)

merged_heatmap <- left_join(
  heatmap_df,
  master_df_filtered %>%
    select(Feature, Treatment, pvalue, log2FoldChange),
  by = "Feature"
)

# optional: reorder columns (clean structure)
merged_heatmap <- merged_heatmap %>%
  relocate(Feature, Treatment, log2FoldChange, pvalue)

# =========================
# MERGE HEATMAP MATRIX + FULL METADATA
# =========================

heatmap_df <- as.data.frame(heatmap_matrix)
heatmap_df$Feature <- rownames(heatmap_df)

merged_heatmap <- left_join(
  heatmap_df,
  master_df_filtered %>%
    select(Feature, Treatment, pvalue, log2FoldChange),
  by = "Feature"
)

# optional: reorder columns (clean structure)
merged_heatmap <- merged_heatmap %>%
  relocate(Feature, Treatment, log2FoldChange, pvalue)
write.csv(
  merged_heatmap,
  file = "heatmap_matrix_with_metadata.csv",
  row.names = FALSE
)



##new heatmap

breaks <- seq(0, 3, length.out = 101)
colors <- colorRampPalette(c("blue", "yellow", "red"))(100)

pheatmap(
  heatmap_matrix,
  
  scale = "row",
  
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  
  annotation_colors = ann_colors,
  
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  
  show_rownames = FALSE,
  show_colnames = FALSE,
  
  border_color = NA,
  
  # 🔥 FORCE SYMMETRIC COLOR DISTRIBUTION
  color = colors,
  breaks = breaks
)



###scale optimization

# =========================
# CUSTOM COLOR SCALE
# =========================

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
# HEATMAP
# =========================

pheatmap(
  heatmap_matrix,
  
  scale = "row",
  
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  
  annotation_col = annotation_col,
  annotation_row = row_annotation,
  
  annotation_colors = ann_colors,
  
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  
  show_rownames = FALSE,
  show_colnames = FALSE,
  
  border_color = NA,
  
  color = colors,
  breaks = breaks
)
