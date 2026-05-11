
#########for volatile data emission bargraph###

setwd("D:/INRS_resources/Wheat metagenome/VOC_data/Deseq/VOC_an_1_f")

Summary <- read_excel("Summary.xlsx", sheet = "graph")

ggplot(Summary, aes(x = Treatment, y = VOCs, fill = Status)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  
  # --- Add white labels inside bars ---
  geom_text(aes(label = VOCs),
            position = position_stack(vjust = 0.5),
            color = "white",
            fontface = "bold",
            size = 4) +
  # -------------------------------------

scale_fill_manual(
  values = c(
    "Emission decreased" = "#1E90FF",
    "Emission increased" = "#DC143C"
  )
) +
  labs(
    x = "Treatment",
    y = "Number of VOCs",
    fill = "Status",
    title = "Microbial VOC Emission"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    #panel.grid.major.x = element_blank()
  )

#############extract compond for ploting ##########

merged_compound_counts <- read_excel("merged_compound_counts.xlsx", 
                                     sheet = "imp_found_database")

merged <-merged_compound_counts

library(dplyr)
library(readr)
library(stringr)
library(purrr)


file_list <- c(
  "DT_75_vs_100_an_1_padj_filtered.csv",
  "DT_50_vs_100_an_1_padj_filtered.csv",
  "DT_25_vs_100_an_1_padj_filtered.csv",
  "DS_75_vs_100_an_1_padj_filtered.csv",
  "DS_50_vs_100_an_1_padj_filtered.csv",
  "DS_25_vs_100_an_1_padj_filtered.csv"
)

#Create a function to read each comparison file with treatment name
read_and_process <- function(file){
  
  df <- read_csv(file, show_col_types = FALSE) %>%
    mutate(
      Treatment = str_remove(file, "_an_1_padj_filtered.csv"),
      Status = ifelse(log2FoldChange > 0, "Increased", "Decreased")
    ) %>%
    select(compound_info, log2FoldChange, Status, Treatment)
  
  return(df)
}

#Read all files into one dataframe
all_results <- map_df(file_list, read_and_process)

#Match compounds to merged_compond_counts

matched_results <- merged %>%
  left_join(all_results, by = "compound_info")


VOC_summary <- matched_results %>%
  filter(!is.na(Status)) %>%
  group_by(Treatment, Status) %>%
  summarise(VOCs = n(), .groups = "drop")





summarize_file <- function(file) {
  
  df <- read_csv(file, show_col_types = FALSE)
  
  df %>%
    mutate(
      Treatment = str_remove(file, "_an_1_padj_filtered.csv"),
      Status = ifelse(log2FoldChange > 0, "Emission increased",
                      "Emission decreased")
    ) %>%
    group_by(Treatment, compound_info, Status) %>%
    summarise(Count = n(), .groups = "drop")
}

# Apply to all files & combine into master long file
Master_Long <- map_df(file_list, summarize_file)


write.csv(Master_Long, "Master_Long_compond_inf.csv")


###known compond mvoc###

library(dplyr)
library(readr)

# Read merged file (containing valid compound_info)
#merged <- read_csv("merged_compound_counts.csv", show_col_types = FALSE)

# Extract compound list
valid_compounds <- unique(merged$Compound_info)

# Filter Master_Long to keep only compounds also in merged file
Master_Filtered <- Master_Long %>%
  filter(compound_info %in% valid_compounds)

Master_Filtered



write.csv(Master_Filtered, "Master_Filtered.csv")


# Count number of unique compound_info entries
n_unique_compounds <- Master_Filtered %>%
  summarise(n = n_distinct(compound_info))

n_unique_compounds


library(dplyr)
library(ggplot2)
library(forcats)

# Step 1 — Compute average count per compound
avg_counts <- Master_Filtered %>%
  group_by(compound_info) %>%
  summarise(avg_count = mean(Count), .groups = "drop") %>%
  filter(avg_count >= 2)   # only compounds avg >= 2

# Step 2 — Filter original data to keep only these compounds
plot_data <- Master_Filtered %>%
  filter(compound_info %in% avg_counts$compound_info)

# Step 3 — Keep only first 6 treatments
top_treatments <- unique(plot_data$Treatment)[1:6]
plot_data <- plot_data %>% filter(Treatment %in% top_treatments)

# Step 4 — Reorder compound_info by total Count
plot_data <- plot_data %>%
  group_by(compound_info) %>%
  mutate(total_count = sum(Count)) %>%
  ungroup() %>%
  mutate(compound_info = fct_reorder(compound_info, total_count))

# Step 5 — Define status colors
status_colors <- c(
  "Emission increased" = "#DC143C",
  "Emission decreased" = "#1E90FF"
)

# Step 6 — Plot with facet_wrap by Status
ggplot(plot_data, aes(x = Treatment, y = compound_info, size = Count, color = Status)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  scale_size_continuous(range = c(3, 12), name = "Count") +
  facet_wrap(~Status, nrow = 1, scales = "free_x") +
  labs(
    x = "Treatment",
    y = "VOC",
    color = "Status",
    #title = "VOC Emission"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1, "lines")
  )



###########revised figures 1st draft#######################
setwd ("D:/D_drive/INRS_resources/Wheat metagenome/VOC_data/Deseq/VOC_an_1_f")

library(readxl)
Master_Long_compond_inf <- read_excel("Master_Long_compond_inf.xlsx")



### MVOC database matching VOC compunds ###

library(dplyr)
library(readr)

# Read merged file (containing valid compound_info)
#merged <- read_csv("merged_compound_counts.csv", show_col_types = FALSE)

# Extract compound list
mvoc_compounds <- unique(merged_compound_counts_f$Compound_info)

# Filter Master_Long to keep only compounds also in merged file
Master_Filtered <- Master_Long_compond_inf %>%
  filter(compound_info %in% mvoc_compounds)

Master_Filtered



write.csv(Master_Filtered, "Master_mvoc_compund.csv")


# Count number of unique compound_info entries
n_unique_compounds <- Master_Filtered %>%
  summarise(n = n_distinct(compound_info))

n_unique_compounds


library(dplyr)
library(ggplot2)
library(forcats)

# Step 1 — Compute average count per compound
avg_counts <- Master_Filtered %>%
  group_by(compound_info) %>%
  summarise(avg_count = mean(Count), .groups = "drop") %>%
  filter(avg_count >= 0)   # only compounds avg >= 2

# Step 2 — Filter original data to keep only these compounds
plot_data <- Master_Filtered %>%
  filter(compound_info %in% avg_counts$compound_info)

# Step 3 — Keep only first 6 treatments
top_treatments <- unique(plot_data$Treatment)[1:6]
plot_data <- plot_data %>% filter(Treatment %in% top_treatments)

# Step 4 — Reorder compound_info by total Count
plot_data <- plot_data %>%
  group_by(compound_info) %>%
  mutate(total_count = sum(Count)) %>%
  ungroup() %>%
  mutate(compound_info = fct_reorder(compound_info, total_count))

# Step 5 — Define status colors
status_colors <- c(
  "Emission increased" = "#DC143C",
  "Emission decreased" = "#1E90FF"
)

# Step 6 — Plot with facet_wrap by Status
ggplot(plot_data, aes(x = Treatment, y = compound_info, size = Count, color = Status)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  scale_size_continuous(range = c(3, 12), name = "Count") +
  facet_wrap(~Status, nrow = 1, scales = "free_x") +
  labs(
    x = "Treatment",
    y = "VOC",
    color = "Status",
    #title = "VOC Emission"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1, "lines")
  )

##############new plots only with emission increased########


library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)

# Step 1 — Filter only "Emission increased" compounds
plot_data <- Master_Filtered %>%
  filter(Status == "Emission increased") %>%
  # Step 2 — Assign facet based on treatment prefix
  mutate(Facet = case_when(
    str_starts(Treatment, "DS") ~ "DS",
    str_starts(Treatment, "DT") ~ "DT",
    TRUE ~ "Other"
  )) 

#Export plot data and clean if mass and charge ratio shown with componds
write.csv(plot_data, "plot_data_emission_increased.csv", row.names = FALSE)

#import clean plot data####
plot_data <- read.csv("D:/D_drive/INRS_resources/Wheat metagenome/VOC_data/Deseq/VOC_an_1_f/plot_data_emission_increased.csv")


# Step 3 — Reorder compounds by total Count
plot_data <- plot_data %>%
  group_by(compound_info) %>%
  mutate(total_count = sum(Count)) %>%
  ungroup() %>%
  mutate(compound_info = fct_reorder(compound_info, total_count))

# Step 4 — Define colors
status_colors <- c(
  "Emission increased" = "#DC143C"
)

# Step 5 — Plot
ggplot(plot_data, aes(x = Treatment, y = compound_info, size = Count, color = Status)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  scale_size_continuous(range = c(3, 12), name = "Count") +
  facet_wrap(~Facet, nrow = 1, scales = "free_x") +
  labs(
    x = "Treatment",
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1, "lines")
  )


##############barplot emission increased#########

library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

Summary <- read_excel("Summary.xlsx", sheet = "graph")

Summary_inc <- Summary %>%
  filter(Status == "Emission increased") %>%
  mutate(Facet = case_when(
    str_starts(Treatment, "DS") ~ "DS",
    str_starts(Treatment, "DT") ~ "DT"
  )) %>%
  filter(!is.na(Facet)) %>%  # keep only DS and DT
  mutate(Facet = factor(Facet, levels = c("DS", "DT")))

ggplot(Summary_inc, aes(x = Treatment, y = VOCs)) +
  geom_bar(stat = "identity", fill = "#DC143C", color = "black", linewidth = 0.2) +
  
  facet_wrap(~Facet, nrow = 1, scales = "free_x") +
  
  labs(
    x = "Treatment",
    y = "Number of DA VOCs",
    title = ""
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1))











# =========================
# FIGURE 7: Emission Increased (A = Bar, B = Dot)
# =========================

library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(readxl)
library(patchwork)
library(officer)
library(rvg)

# -------------------------
# Step 1 — Prepare bubble plot data
# -------------------------
plot_data <- plot_data %>%
  mutate(
    compound_info = str_trim(compound_info)  # removes leading/trailing spaces
  )


write.csv(plot_data, "plot_data_emission_increased.csv", row.names = FALSE)

plot_data <- plot_data %>%
  group_by(compound_info) %>%
  mutate(total_count = sum(Count)) %>%
  ungroup() %>%
  mutate(compound_info = fct_reorder(compound_info, total_count))

# -------------------------
# Step 2 — Dot plot (Figure 7B)
# -------------------------
status_colors <- c("Emission increased" = "#DC143C")

plot_dot <- ggplot(plot_data, aes(x = Treatment, y = compound_info, size = Count, color = Status)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = status_colors) +
  scale_size_continuous(range = c(3, 12), name = "Count") +
  facet_wrap(~Facet, nrow = 1, scales = "free_x") +
  labs(x = "Treatment", y = "") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.spacing = unit(1, "lines")
  )

# -------------------------
# Step 3 — Bar plot (Figure 7A)
# -------------------------
Summary <- read_excel("Summary.xlsx", sheet = "graph")

Summary_inc <- Summary %>%
  filter(Status == "Emission increased") %>%
  mutate(Facet = case_when(
    str_starts(Treatment, "DS") ~ "DS",
    str_starts(Treatment, "DT") ~ "DT"
  )) %>%
  filter(!is.na(Facet)) %>%
  mutate(Facet = factor(Facet, levels = c("DS", "DT")))

write.csv(Summary_inc, "barplot_DA_VOCs.csv", row.names = FALSE)

plot_bar <- ggplot(Summary_inc, aes(x = Treatment, y = VOCs)) +
  geom_bar(stat = "identity", fill = "#DC143C", color = "black", linewidth = 0.2) +
  facet_wrap(~Facet, nrow = 1, scales = "free_x") +
  labs(x = "Treatment", y = "Number of emitted VOCs") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------------
# Step 4 — Combined Figure (A on top, B below)
# -------------------------
combined_fig7 <- plot_bar / plot_dot +
  plot_layout(heights = c(1, 2.5)) +
  plot_annotation(tag_levels = "A")

# -------------------------
# Step 5 — Save as SVG
# -------------------------
ggsave("Figure7A_Bar.svg", plot_bar, width = 10, height = 5)
ggsave("Figure7B_Dot.svg", plot_dot, width = 10, height = 7)
ggsave("Figure7_Combined_AB.svg", combined_fig7, width = 10, height = 7.5)

# -------------------------
# Step 6 — Export to PowerPoint
# -------------------------
plot_list7 <- list(plot_bar, plot_dot, combined_fig7)

ppt <- read_pptx()

for(p in plot_list7){
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  
  if(identical(p, combined_fig7)){
    ppt <- ph_with(ppt, dml(ggobj = p),
                   location = ph_location(left = 0, top = 0, width = 10, height = 7.5))
  } else {
    ppt <- ph_with(ppt, dml(ggobj = p),
                   location = ph_location_fullsize())
  }
}

print(ppt, target = "Figure7_Emission_Increased_AB.pptx")





# =========================
# FIGURE 7 update: Combined (A = Bar, B = Line, C = Dot) - Fixed Colors and Layout
# =========================

library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(readxl)
library(patchwork)
library(officer)
library(rvg)

# -------------------------
# Step 1 — Prepare Dot Plot Data (C)
# -------------------------
Figure7_ABC_Data <- read_excel("Figure7_ABC_Data.xlsx", 
                               +     sheet = "C")

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
# Step 2 — Prepare Bar Plot (A)
# -------------------------
Summary_bar <- read_excel("Summary.xlsx", sheet = "graph") %>%
  mutate(
    Status = str_trim(Status),
    Status = case_when(
      Status %in% c("Emission increased", "Production") ~ "Production",
      Status %in% c("Emission decreased", "Consumption") ~ "Consumption",
      TRUE ~ Status
    ),
    Status = factor(Status, levels = c("Production", "Consumption"))
  )

plot_bar <- ggplot(Summary_bar, aes(x = Treatment, y = VOCs, fill = Status)) +
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
# Step 3 — Prepare Line Graph (B)
# -------------------------
Summary_line <- read_excel("Summary.xlsx", sheet = "line_graph") %>%
  mutate(Treatment_num = str_extract(Treatment, "\\d+"),
         Treatment_num = factor(Treatment_num, levels = c("25","50","75")))

plot_line <- ggplot(Summary_line, aes(x = Treatment_num, y = Percent_VOCs, group = Genotype, color = Genotype)) +
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
# Step 4 — Combine Figures
# -------------------------
# Top row: Bar (A) + Line (B)
top_row <- plot_bar + plot_line + plot_layout(widths = c(1, 1)) 

# Full Figure 7: Top row + Dot plot (C) below (increase C height)
figure7_combined <- top_row / plot_dot +
  plot_layout(heights = c(1, 1.1)) +  # increase height for Dot plot
  plot_annotation(tag_levels = "A")

# -------------------------
# Step 5 — Save SVG files
# -------------------------
ggsave("Figure7A_Bar.svg", plot_bar, width = 7, height = 5)
ggsave("Figure7B_Line.svg", plot_line, width = 7, height = 5)
ggsave("Figure7C_Dot.svg", plot_dot, width = 10, height = 7)
ggsave("Figure7_Combined_ABC_2.svg", figure7_combined, width = 10, height = 10)

# -------------------------
# Step 6 — Export to PowerPoint
# -------------------------
plot_list7 <- list(plot_bar, plot_line, plot_dot, figure7_combined)

ppt <- read_pptx()

for(p in plot_list7){
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  
  if(identical(p, figure7_combined)){
    ppt <- ph_with(ppt, dml(ggobj = p),
                   location = ph_location(left = 0, top = 0, width = 9, height = 9))
  } else {
    ppt <- ph_with(ppt, dml(ggobj = p),
                   location = ph_location_fullsize())
  }
}

print(ppt, target = "Figure7_ABC_fixed_n.pptx")


######save figure 7 data #########

library(readxl)
library(dplyr)
library(openxlsx)  # for writing multiple sheets

# =========================
# Step 1 — Prepare Bar Plot Data (A)
# =========================
Summary_bar <- read_excel("Summary.xlsx", sheet = "graph") %>%
  mutate(
    Status = str_trim(Status),
    Status = case_when(
      Status %in% c("Emission increased", "Production") ~ "Production",
      Status %in% c("Emission decreased", "Consumption") ~ "Consumption",
      TRUE ~ Status
    ),
    Status = factor(Status, levels = c("Production", "Consumption"))
  )

# =========================
# Step 2 — Prepare Line Graph Data (B)
# =========================
Summary_line <- read_excel("Summary.xlsx", sheet = "line_graph") %>%
  mutate(Treatment_num = str_extract(Treatment, "\\d+"),
         Treatment_num = factor(Treatment_num, levels = c("25","50","75")))

# =========================
# Step 3 — Prepare Dot Plot Data (C)
# =========================
plot_data_dot <- plot_data

# =========================
# Step 4 — Export All to ONE Excel file
# =========================
wb <- createWorkbook()

# Add Bar Plot Data
addWorksheet(wb, "Bar")
writeData(wb, "Bar", Summary_bar)

# Add Line Graph Data
addWorksheet(wb, "Line")
writeData(wb, "Line", Summary_line)

# Add Dot Plot Data
addWorksheet(wb, "Dot")
writeData(wb, "Dot", plot_data_dot)

# Save workbook
saveWorkbook(wb, "Figure7_ABC_Data.xlsx", overwrite = TRUE)
