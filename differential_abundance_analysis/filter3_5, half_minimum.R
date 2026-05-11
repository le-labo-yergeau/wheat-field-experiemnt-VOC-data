library(dplyr)

# 1. Remove non-count column if present
count_table <- data_org_n %>%
  select(-Treatment)

# 2. Convert to numeric (keep NA)
count_table <- count_table %>%
  mutate(across(everything(), as.numeric))

# 3. Define groups (e.g., 25_DS, 100_DT)
group_names <- sub("_[0-9]+$", "", colnames(count_table))
unique_groups <- unique(group_names)

# 4. FILTER FIRST (ignore NA)
keep_features <- apply(count_table, 1, function(x) {
  any(sapply(unique_groups, function(g) {
    group_idx <- which(group_names == g)
    sum(x[group_idx] > 0, na.rm = TRUE) >= 3
  }))
})

data_filtered <- count_table[keep_features, ]

# 5. Replace NA with half of the minimum (row-wise)
data_imputed <- t(apply(data_filtered, 1, function(x) {
  min_val <- min(x[x > 0], na.rm = TRUE)
  
  # If no positive values exist (rare after filtering), keep as is
  if (is.infinite(min_val)) {
    return(x)
  } else {
    x[is.na(x)] <- min_val / 2
    return(x)
  }
}))

# Convert back to data frame and restore column names
data_imputed <- as.data.frame(data_imputed)
colnames(data_imputed) <- colnames(data_filtered)

# Check
dim(data_imputed)