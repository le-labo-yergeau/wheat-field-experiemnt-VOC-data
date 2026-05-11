library(dplyr)

setwd("D:/D_drive/INRS_resources/Wheat metagenome/VOC_data/Deseq/VOC_an_2")

Field_exp_wheat_alignment0_5_sorted_nist_metabo <- read_excel("Field_exp_wheat_alignment0.5_sorted_nist_metabo.xlsx", 
                                                              sheet = "Field_exp_wheat_alignment0.5_so")
data_org_n<-Field_exp_wheat_alignment0_5_sorted_nist_metabo

meta<- read_excel("meta.xlsx", sheet = "meta")
# 1. Remove the first row if it contains metadata like "Time"
#data_org_n <- data_org_n[-1, ]  

# 2. Move the first column as row names
data_org_n <- data_org_n %>% column_to_rownames(var = colnames(data_org_n)[1])


# 3. Convert all columns to numeric
data_org_n <- data_org_n %>% mutate(across(everything(), as.numeric))

#Now data_org_n is a numeric matrix-like data frame
str(data_org_n)

# 2. Filter features with >20% missing values (keep ≤20% NA)
na_fraction <- rowMeans(is.na(data_org_n))
data_filtered <- data_org_n[na_fraction <= 0.90, ]


# 3. Half-minimum imputation (row-wise per VOC)
data_imputed <- t(apply(data_filtered, 1, function(x) {
  
  # If all values are NA, keep as NA (or remove earlier if preferred)
  if(all(is.na(x))) return(x)
  
  # Find minimum positive value
  min_val <- min(x[x > 0], na.rm = TRUE)
  
  # Replace NA with half-min
  x[is.na(x)] <- min_val / 2
  
  return(x)
}))

# 4. Convert back to data frame
data_org_n <- as.data.frame(data_imputed)


#data_org_0 <-data_org
# Subset only DS samples
ds_data <- data_org_n[, grepl("_DS_", colnames(data_org_n))]

# Subset only DT samples
dt_data <- data_org_n[, grepl("_DT_", colnames(data_org_n))]


# Subset DS samples
meta_DS <- meta %>% filter(Compartment == "DS")

# Subset DT samples
meta_DT <- meta %>% filter(Compartment == "DT")



# Ensure sample order matches count matrix
meta_DS <- meta_DS %>% filter(ID %in% colnames(ds_data))
ds_data <- ds_data[, meta_DS$ID]  # Reorder columns to match metadata order


# Convert Treatment to factor with 100 as reference level
meta_DS$Treatment <- as.character(meta_DS$Treatment)

# Then convert to factor with correct levels (keep all levels that exist!)
meta_DS$Treatment <- factor(meta_DS$Treatment, levels = c("100", "75", "50", "25"))

# Optional: set reference level (if you want 100 as baseline)
meta_DS$Treatment <- relevel(meta_DS$Treatment, ref = "100")

# Check result
table(meta_DS$Treatment)


#  Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(ds_data),  # counts must be integers
  colData   = meta_DS,
  design    = ~ Treatment
)

#  Filter low-count features (optional, improves power)
#dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

#  Get results (25 vs 100 comparison)
res <- results(dds, contrast = c("Treatment", "25", "100"))

# Order by adjusted p-value
res <- res[order(res$padj), ]

#  Save results
write.csv(as.data.frame(res), "DESeq2_DS_25_vs_100_results_n1.csv")

#  Quick summary
summary(res)

# Plot MA
plotMA(res, ylim = c(-5, 5))



# ---- Extract results for multiple comparisons ----
# Get results and order by adjusted p-value for all comparisons
res_25_vs_100 <- results(dds, contrast = c("Treatment", "25", "100"))
res_25_vs_100 <- res_25_vs_100[order(res_25_vs_100$padj), ]

res_50_vs_100 <- results(dds, contrast = c("Treatment", "50", "100"))
res_50_vs_100 <- res_50_vs_100[order(res_50_vs_100$padj), ]

res_75_vs_100 <- results(dds, contrast = c("Treatment", "75", "100"))
res_75_vs_100 <- res_75_vs_100[order(res_75_vs_100$padj), ]
# Optional: save results
write.csv(as.data.frame(res_25_vs_100), "DS_25_vs_100_an_1.csv")
write.csv(as.data.frame(res_50_vs_100), "DS_50_vs_100_an_1.csv")
write.csv(as.data.frame(res_75_vs_100), "DS_75_vs_100_an_1.csv")

# Function to summarize DESeq2 results
get_sig_summary <- function(res, dataset_name) {
  up <- sum(res$log2FoldChange > 0 & res$padj <= 0.05, na.rm = TRUE)
  down <- sum(res$log2FoldChange < 0 & res$padj <= 0.05, na.rm = TRUE)
  ns <- sum(res$padj > 0.05, na.rm = TRUE)
  
  data.frame(
    Comparison = dataset_name,
    Upregulated = up,
    Downregulated = down,
    Non_significant = ns
  )
}

# Create summary table for all three comparisons
summary_table <- do.call(rbind, list(
  get_sig_summary(res_25_vs_100, "25_vs_100"),
  get_sig_summary(res_50_vs_100, "50_vs_100"),
  get_sig_summary(res_75_vs_100, "75_vs_100")
))

# Print summary table
print(summary_table)

# Optionally save to CSV
write.csv(summary_table, "DS_VOC_summary_an_2.csv", row.names = FALSE)

# Optionally save as RDS
saveRDS(summary_table, "DS_VOC_summary.rds")




##########DT################
# Ensure sample order matches count matrix
meta_DT <- meta_DT %>% filter(ID %in% colnames(dt_data))
dt_data <- dt_data[, meta_DT$ID]  # Reorder columns to match metadata order


# Convert Treatment to factor with 100 as reference level
meta_DT$Treatment <- as.character(meta_DT$Treatment)

# Then convert to factor with correct levels (keep all levels that exist!)
meta_DT$Treatment <- factor(meta_DT$Treatment, levels = c("100", "75", "50", "25"))

# Optional: set reference level (if you want 100 as baseline)
meta_DT$Treatment <- relevel(meta_DT$Treatment, ref = "100")

# Check result
table(meta_DT$Treatment)


#  Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(dt_data),  # counts must be integers
  colData   = meta_DT,
  design    = ~ Treatment
)

#  Filter low-count features (optional, improves power)
#dds <- dds[rowSums(counts(dds)) > 10, ]

gm_mean <- function(x) {
  if(all(x == 0)) return(0)
  exp(mean(log(x[x > 0])))
}

geoMeans <- apply(counts(dds), 1, gm_mean)

dds <- estimateSizeFactors(dds, geoMeans = geoMeans)


# Run DESeq2
dds <- DESeq(dds)




# ---- Extract results for multiple comparisons ----
res_25_vs_100 <- results(dds, contrast = c("Treatment", "25", "100"))
res_50_vs_100 <- results(dds, contrast = c("Treatment", "50", "100"))
res_75_vs_100 <- results(dds, contrast = c("Treatment", "75", "100"))

# Optional: save results
write.csv(as.data.frame(res_25_vs_100), "DT_25_vs_100_an1.csv")
write.csv(as.data.frame(res_50_vs_100), "DT_50_vs_100_an1.csv")
write.csv(as.data.frame(res_75_vs_100), "DT_75_vs_100_an1.csv")






# Function to summarize DESeq2 results
get_sig_summary <- function(res, dataset_name) {
  up <- sum(res$log2FoldChange > 0 & res$padj <= 0.05, na.rm = TRUE)
  down <- sum(res$log2FoldChange < 0 & res$padj <= 0.05, na.rm = TRUE)
  ns <- sum(res$padj > 0.05, na.rm = TRUE)
  
  data.frame(
    Comparison = dataset_name,
    Upregulated = up,
    Downregulated = down,
    Non_significant = ns
  )
}

# Create summary table for all three comparisons
summary_table <- do.call(rbind, list(
  get_sig_summary(res_25_vs_100, "25_vs_100"),
  get_sig_summary(res_50_vs_100, "50_vs_100"),
  get_sig_summary(res_75_vs_100, "75_vs_100")
))

# Print summary table
print(summary_table)

# Optionally save to CSV
write.csv(summary_table, "DT_VOC_summary_an_1.csv", row.names = FALSE)

# Optionally save as RDS
saveRDS(summary_table, "DT_VOC_summary.rds")





###############DS_vs_DT####################

# Ensure sample order matches count matrix
#meta <- meta %>% filter(ID %in% colnames(data_filtered))
#data_filtered <- data_filtered[, data_filtered$ID]  # Reorder columns to match metadata order


# Convert Treatment to factor with 100 as reference level
meta$Compartment <- as.character(meta$Compartment)

# Then convert to factor with correct levels (keep all levels that exist!)
meta$Compartment <- factor(meta$Compartment, levels = c("DS", "DT"))

# Optional: set reference level (if you want 100 as baseline)
meta$Compartment <- relevel(meta$Compartment, ref = "DS")

# Check result
table(meta$Compartment)

ncol(data_filtered) == nrow(meta)
#  Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(data_filtered),  # counts must be integers
  colData   = meta,
  design    = ~ Compartment
)



#  Filter low-count features (optional, improves power)
#dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

#  Get results (25 vs 100 comparison)
res <- results(dds, contrast = c("Compartment", "DT", "DS"))

# Order by adjusted p-value
res <- res[order(res$padj), ]

#  Save results
write.csv(as.data.frame(res), "DS_vs_DT_results_n.csv")

#  Quick summary
summary(res)

# Plot MA
plotMA(res, ylim = c(-5, 5))




# Function to summarize DESeq2 results
get_sig_summary <- function(res, dataset_name) {
  up <- sum(res$log2FoldChange > 0 & res$padj <= 0.05, na.rm = TRUE)
  down <- sum(res$log2FoldChange < 0 & res$padj <= 0.05, na.rm = TRUE)
  ns <- sum(res$padj > 0.05, na.rm = TRUE)
  
  data.frame(
    Comparison = dataset_name,
    Upregulated = up,
    Downregulated = down,
    Non_significant = ns
  )
}

# Create summary table for all three comparisons
summary_table <- do.call(rbind, list(
  get_sig_summary(res, "DT_vs_DS")))

# Print summary table
print(summary_table)

# Optionally save to CSV
write.csv(summary_table, "DT_vs_DS.csv", row.names = FALSE)

# Optionally save as RDS
saveRDS(summary_table, "DS_VOC_summary.rds")



############# --- Replace 0 with half of the minimum positive value ---


data_filtered<-ds_data
# Find minimum positive value
min_positive <- min(data_filtered[data_filtered > 0], na.rm = TRUE)

# Compute replacement value (half of min positive)
replacement_value <- min_positive / 2

# Replace zeros
data_replaced <- data_filtered
data_replaced[data_replaced == 0] <- replacement_value

# Show confirmation
message("Replaced all zeros with: ", replacement_value)

ds_data<-data_replaced







################annotation extract

# Load libraries
library(dplyr)
library(stringr)
library(readr)
library(tools)

# -------------------------------
# Step 1: Import your VOC list
# -------------------------------
# Load libraries
library(dplyr)
library(stringr)
library(readr)
library(tools)

# -------------------------------
# Step 1: Import your VOC list
# -------------------------------
voc_table <- read_csv("voc_list.csv")   # columns: "Compound" and "Category"
voc_list <- voc_table$Compound

# -------------------------------
# Step 2: Automatically find all target files
# -------------------------------
files <- list.files(pattern = "_an_1\\.csv$")

# -------------------------------
# Step 3: Loop through each file
# -------------------------------
for (f in files) {
  
  # Read the CSV data correctly
  df <- read_csv(f)
  
  # Rename the first column to 'Feature' if needed
  colnames(df)[1] <- "Feature"
  
  # Step 4: Split the Feature into mz_rt and compound_info
  df <- df %>%
    mutate(
      mz_rt = str_extract(Feature, "^[^\\(]+") %>% str_trim(),
      compound_info = str_extract(Feature, "\\(([^\\)]+)\\)") %>%
        str_replace_all("[\\(\\)]", "") %>%
        str_trim()
    )
  
  # Step 5: Filter rows based on compounds in voc_list using literal matching
  df_voc <- df %>%
    filter(Reduce(`|`, lapply(voc_list, function(x) str_detect(compound_info, fixed(x))))) %>%
    # Step 5b: Filter padj <= 0.05
    filter(!is.na(padj) & padj <= 0.05)
  
  # Step 6: Add VOC category info
  df_voc <- df_voc %>%
    left_join(voc_table, by = c("compound_info" = "Compound"))
  
  # Step 7: Export as new CSV
  output_name <- paste0(file_path_sans_ext(f), "_VOC_filtered.csv")
  write_csv(df_voc, output_name)
  
  message("✅ Processed: ", f, " → ", output_name, " (padj ≤ 0.05 applied)")
}



Export padj-filtered original file

df_padj <- df %>%
  filter(!is.na(padj) & padj <= 0.05)

output_padj <- paste0(file_path_sans_ext(f), "_padj_filtered.csv")
write_csv(df_padj, output_padj)
