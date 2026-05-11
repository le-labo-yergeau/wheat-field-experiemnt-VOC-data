library(dplyr)
library(tibble)
library(missForest)

# 1. Load data

Field_exp_wheat_alignment0_5_sorted_nist_metabo <- read_excel("Field_exp_wheat_alignment0.5_sorted_nist_metabo.xlsx", 
                                                              sheet = "Field_exp_wheat_alignment0.5_so")

data_org_n <- Field_exp_wheat_alignment0_5_sorted_nist_metabo

# 2. Move first column (Treatment) to rownames
data_org_n <- data_org_n %>%
  column_to_rownames(var = colnames(data_org_n)[1])

# 3. Convert to numeric safely
count_table <- data_org_n %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))

str(count_table)

# 4. Define groups
group_names <- sub("_[0-9]+$", "", colnames(count_table))
unique_groups <- unique(group_names)

# 5. Replicate-based filtering (≥2 detections per group)
keep_features <- apply(count_table, 1, function(x) {
  any(sapply(unique_groups, function(g) {
    group_idx <- which(group_names == g)
    sum(x[group_idx] > 0, na.rm = TRUE) >= 3
  }))
})

data_filtered <- count_table[keep_features, , drop = FALSE]

# 6. IMPORTANT: transpose for missForest (it works column-wise better)
rf_input <- as.data.frame(t(data_filtered))

# 7. Random Forest imputation
set.seed(123)
rf_imputed <- missForest(rf_input, maxiter = 10, ntree = 100)

# 8. Extract imputed data
data_imputed <- t(rf_imputed$ximp)

# 9. Convert back to data frame
data_imputed <- as.data.frame(data_imputed)

# 10. Restore names
colnames(data_imputed) <- colnames(data_filtered)
rownames(data_imputed) <- rownames(data_filtered)

# 11. Check
dim(data_imputed)