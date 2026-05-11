library(dplyr)

Field_exp_wheat_alignment0_5_sorted_nist_metabo <- read_excel("Field_exp_wheat_alignment0.5_sorted_nist_metabo.xlsx", 
                                                              sheet = "Field_exp_wheat_alignment0.5_so")
data_org_n<-Field_exp_wheat_alignment0_5_sorted_nist_metabo

# 1. Remove non-count column if still present
count_table <- data_org_n %>% dplyr::select(-Treatment)

# 2. Convert to numeric
count_table <- count_table %>%
  dplyr::mutate(across(everything(), as.numeric))

# 3. Replace NA with 0
count_table[is.na(count_table)] <- 0

# 4. Extract group names properly (e.g., "25_DS", "100_DT")
group_names <- sub("_[0-9]+$", "", colnames(count_table))
unique_groups <- unique(group_names)

# 5. Filter: present in ≥3 replicates within ≥1 group
keep_features <- apply(count_table, 1, function(x) {
  any(sapply(unique_groups, function(g) {
    group_idx <- which(group_names == g)
    sum(x[group_idx] > 0) >= 3
  }))
})

# 6. Apply filter
filtered_table <- count_table[keep_features, ]

# Check
dim(filtered_table)

