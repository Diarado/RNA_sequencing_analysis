# ------------------------------------------------------------------------------
# Author: Jiahe Tian
# Date: 8/13/2024
# Purpose: It correlates individual gene counts across all samples with RIN and 
# PMI (if applicable) in NPH, NYBB, ROSMAP_polyA and ROSMAP_riboDepletion datasets
# ------------------------------------------------------------------------------
# Input:
#   - /rna_quality_paper/correlation_on_each_gene (supp 1&2)/variables.RData
# Output:
#   - /rna_quality_paper/correlation_on_each_gene (supp 1&2)/
# ------------------------------------------------------------------------------

library(psych)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(limma)

# ------------------------------------------------------------------------------
# Data Loading
# ------------------------------------------------------------------------------

# please specify your working directory
setwd("D:/Jiahe/Columbia")

# please specify whether to regress out sex and age
regress_age_sex = FALSE

load("rna_quality_paper/correlation_on_each_gene (supp 1&2)/processed_data.RData")

# ------------------------------------------------------------------------------
# NPH
# ------------------------------------------------------------------------------
# Filter out low-count genes
NPH_counts <- NPH_counts[!rowSums(NPH_counts < 5) >= 0.9 * ncol(NPH_counts), ]

NPH_counts <- NPH_counts[, !colnames(NPH_counts) %in% c("001.1", "002.1")]
# Convert count data to matrix
countData <- as.matrix(NPH_counts)
sampleNames <- colnames(countData)

# Create colData with dummy variable since you need it for DESeq2
colData <- data.frame(row.names = sampleNames, dummy = rep(1, length(sampleNames)))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 1)

# Apply variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Extract transformed data
vst_counts <- assay(vsd)

if (regress_age_sex) {
  # Convert Gender to numeric
  ID_age_sex <- NPH_md %>%
    select(SampleID, Age, Gender) %>%
    mutate(Gender = ifelse(Gender == "M", 1, 0))
  
  # Ensure the order of samples in vst_counts matches the order in ID_age_sex
  ID_age_sex <- ID_age_sex[match(colnames(vst_counts), ID_age_sex$SampleID), ]
  
  # Remove batch effect using Age and Gender as covariates
  batch_effects <- removeBatchEffect(vst_counts, batch = NULL, covariates = as.matrix(ID_age_sex[, c("Age", "Gender")]))
  # batch_effects now contains the data with Age and Gender regressed out
  normalized_NPH_counts <- batch_effects
} else {
  normalized_NPH_counts <- vst_counts
}

# Transpose normalized counts matrix
normalized_NPH_counts_t <- t(normalized_NPH_counts)
normalized_NPH_counts_t <- normalized_NPH_counts_t[!duplicated(rownames(normalized_NPH_counts_t)), ]

# Extract RIN values
first_001_index <- which(NPH_md$SampleID == "001")[1]
NPH_md <- NPH_md[-first_001_index, ]

# Find the first occurrence of '002' and change it to '2'
first_002_index <- which(NPH_md$SampleID == "002")[1]
NPH_md <- NPH_md[-first_002_index, ]

NPH_RIN <- data.frame(RIN = NPH_md$RIN)
row.names(NPH_RIN) <- NPH_md$SampleID
NPH_RIN <- NPH_RIN %>%
  mutate(RIN = ifelse(is.na(RIN), mean(RIN, na.rm = TRUE), RIN))


df1 <- data.frame(normalized_NPH_counts_t)
df1$SampleID <- rownames(normalized_NPH_counts_t)

df2 <- data.frame(NPH_RIN)
df2$SampleID <- rownames(NPH_RIN)
NPH_ready <- left_join(df2, df1, by = "SampleID")
NPH_ready <- NPH_ready |>
  select(-SampleID)

cor_results <- corr.test(
  x = NPH_ready[, c("RIN")],
  y = NPH_ready[, !(names(NPH_ready) %in% c("RIN"))],
  method = "spearman",
  adjust = "fdr"
)

# Extract correlation matrix
cor_matrix <- cor_results$r

# Extract adjusted p-values
p_values <- cor_results$p
padj_values <- cor_results$p.adj

# Combine the matrices into one data frame
combined_df <- data.frame(
  RIN = cor_matrix[1, ],
  RIN_p = p_values[1, ],
  RIN_p_adj = padj_values[1, ]
)

# Split the data frame into two separate tables
RIN_table <- combined_df[, c("RIN", "RIN_p", "RIN_p_adj")]
colnames(RIN_table) <- c("RIN", "p", "p.adj")
write.csv(RIN_table, "rna_quality_paper/correlation_on_each_gene (supp 1&2)/NPH_all_gene_RIN_corr.csv", row.names = TRUE)

# ------------------------------------------------------------------------------
# NYBB
# ------------------------------------------------------------------------------

BB_counts[is.na(BB_counts)] <- 0
BB_counts[!sapply(BB_counts, is.numeric)] <- 0
BB_counts <- BB_counts[!rowSums(BB_counts < 5) >= 0.9 * ncol(BB_counts), ]

age_sex_info_BB <- read.csv("rna_quality_paper/data/Inventory3212023.csv") 
BB_md <- BB_md |>
  mutate(Client.Sample = str_replace(Client.Sample, "4839-2", "4839"), Client.Sample = str_replace(Client.Sample, "4306-2", "4306"),
         Client.Sample = str_replace_all(Client.Sample, "[A-Za-z]", ""))

BB_md$Client.Sample <- as.numeric(BB_md$Client.Sample)

age_sex_info_BB <- age_sex_info_BB %>%
  select(Brain.ID, AgeYear, Gender) |>
  filter(Brain.ID %in% BB_md$Client.Sample)

# VST (including normalization)
countData <- as.matrix(BB_counts)  
sampleNames <- colnames(countData)
colData <- data.frame(row.names = sampleNames, dummy = rep(1, length(sampleNames)))

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ 1)

# Apply variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Extract transformed data
vst_counts <- assay(vsd)

if (regress_age_sex) {
  # Convert Gender to numeric
  age_sex_info_BB <- age_sex_info_BB %>%
    select(Brain.ID, AgeYear, Gender) %>%
    mutate(Gender = ifelse(Gender == "Man", 1, 0)) |>
    drop_na()
  
  # Identify common samples
  common_samples <- intersect(colnames(vst_counts), age_sex_info_BB$Brain.ID)
  
  # Subset vst_counts and age_sex_info_BB to include only the common samples
  vst_counts <- vst_counts[, common_samples]
  age_sex_info_BB <- age_sex_info_BB[age_sex_info_BB$Brain.ID %in% common_samples, ]
  
  # Align the order of samples
  age_sex_info_BB <- age_sex_info_BB[match(colnames(vst_counts), age_sex_info_BB$Brain.ID), ]
  
  # Remove batch effect using Age and Gender as covariates
  batch_effects <- removeBatchEffect(vst_counts, batch = NULL, covariates = as.matrix(age_sex_info_BB[, c("AgeYear", "Gender")]))
  
  
  # batch_effects now contains the data with Age and Gender regressed out
  normalized_BB_counts <- batch_effects
} else {
  normalized_BB_counts <- vst_counts
}

# Transpose normalized counts matrix
normalized_BB_counts_t <- t(normalized_BB_counts)
normalized_BB_counts_t <- normalized_BB_counts_t[!duplicated(rownames(normalized_BB_counts_t)), ]

# Extract RIN values
BB_md <- BB_md |>
  distinct(Client.Sample, .keep_all = TRUE)

BB_RIN_PMI <- data.frame(RIN = BB_md$BA.RIN)
row.names(BB_RIN_PMI) <- BB_md$Client.Sample

matching_rows <- row.names(BB_RIN_PMI) %in% BB_pmi$T

# Extract PMIRTFzn values based on the matching rows
new_row_values <- BB_pmi$PMIRTFzn[match(row.names(BB_RIN_PMI), BB_pmi$T)]

# Create a new row name
new_row_name <- "PMIRTFzn"

# Add the new row to BB_RIN_PMI
BB_RIN_PMI[new_row_name] <- new_row_values
BB_RIN_PMI <- BB_RIN_PMI[!is.nan(BB_RIN_PMI$RIN) & !is.nan(BB_RIN_PMI$PMIRTFzn) & !is.na(BB_RIN_PMI$RIN) & !is.na(BB_RIN_PMI$PMIRTFzn), ]


df1 <- data.frame(normalized_BB_counts_t)
df1$SampleID <- rownames(normalized_BB_counts_t)

df2 <- data.frame(BB_RIN_PMI)
df2$SampleID <- rownames(BB_RIN_PMI)
BB_ready <- left_join(df2, df1, by = "SampleID")
BB_ready <- BB_ready |>
  select(-SampleID)

# print(names(BB_ready))
cor_results <- corr.test(
  x = BB_ready[, c("RIN", "PMIRTFzn")],
  y = BB_ready[, !(names(BB_ready) %in% c("RIN", "PMIRTFzn"))],
  method = "spearman",
  adjust = "fdr"
)


# Extract correlation matrix
cor_matrix <- cor_results$r

# Extract adjusted p-values
p_values <- cor_results$p
padj_values <- cor_results$p.adj

# Combine the matrices into one data frame
combined_df <- data.frame(
  PMI = cor_matrix["PMIRTFzn", ],
  PMI_p = p_values["PMIRTFzn", ],
  PMI_p_adj = padj_values["PMIRTFzn", ],
  RIN = cor_matrix["RIN", ],
  RIN_p = p_values["RIN", ],
  RIN_p_adj = padj_values["RIN", ]
)

# Split the data frame into two separate tables
PMI_table <- combined_df[, c("PMI", "PMI_p", "PMI_p_adj")]
colnames(PMI_table) <- c("PMI", "p", "p.adj")
write.csv(PMI_table, "rna_quality_paper/correlation_on_each_gene (supp 1&2)/BB_all_gene_PMI_corr.csv", row.names = TRUE)
RIN_table <- combined_df[, c("RIN", "RIN_p", "RIN_p_adj")]
colnames(RIN_table) <- c("RIN", "p", "p.adj")
write.csv(RIN_table, "rna_quality_paper/correlation_on_each_gene (supp 1&2)/BB_all_gene_RIN_corr.csv", row.names = TRUE)

# ------------------------------------------------------------------------------
# ROSMAP_oligo and ROSMAP_riboDepletion
# ------------------------------------------------------------------------------

ROSMAP_counts[is.na(ROSMAP_counts)] <- 0
ROSMAP_counts[!sapply(ROSMAP_counts, is.numeric)] <- 0

# filter to find oligo batch
projid_keeps <- ROSMAP_md |>
  pull(projid)

ROSMAP_counts_oligo <- ROSMAP_counts[, colnames(ROSMAP_counts) %in% projid_keeps]
ROSMAP_counts_oligo <- ROSMAP_counts_oligo[!rowSums(ROSMAP_counts_oligo < 5) >= 0.9 * ncol(ROSMAP_counts_oligo), ]

ROSMAP_counts_riboDep <- ROSMAP_counts[, !(colnames(ROSMAP_counts) %in% projid_keeps)]
ROSMAP_counts_riboDep <- ROSMAP_counts_riboDep[!rowSums(ROSMAP_counts_riboDep < 5) >= 0.9 * ncol(ROSMAP_counts_riboDep), ]

sex_info <- read.csv("rna_quality_paper/data/dataset_1267_basic_10-02-2023.csv")

generate_tables <- function(md, counts, name) {
  ROSMAP_counts <- counts
  ROSMAP_md <- md
  age_death_data <- ROSMAP_md %>%
    select(projid, age_death)
  
  ID_age_sex_ROSMAP <- age_death_data %>%
    inner_join(sex_info %>% select(projid, msex), by = "projid")
  
  
  # VST (including normalization)
  countData <- as.matrix(ROSMAP_counts)  
  sampleNames <- colnames(countData)
  colData <- data.frame(row.names = sampleNames, dummy = rep(1, length(sampleNames)))
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ 1)
  
  # Apply variance stabilizing transformation
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  # Extract transformed data
  vst_counts <- assay(vsd)
  
  if (regress_age_sex) {
    ID_age_sex_ROSMAP <- ID_age_sex_ROSMAP[match(colnames(vst_counts), ID_age_sex_ROSMAP$projid), ]
    
    # Remove batch effect using Age and Gender as covariates
    batch_effects <- removeBatchEffect(vst_counts, batch = NULL, covariates = as.matrix(ID_age_sex_ROSMAP[, c("age_death", "msex")]))
    
    # batch_effects now contains the data with Age and Gender regressed out
    normalized_ROSMAP_counts <- batch_effects
  } else {
    normalized_ROSMAP_counts <- vst_counts
  }
  
  # Transpose normalized counts matrix
  normalized_ROSMAP_counts_t <- t(normalized_ROSMAP_counts)
  normalized_ROSMAP_counts_t <- normalized_ROSMAP_counts_t[!duplicated(rownames(normalized_ROSMAP_counts_t)), ]
  
  # Extract RIN values
  ROSMAP_RIN_PMI <- data.frame(RIN = ROSMAP_md$RIN, PMI = ROSMAP_md$pmi)
  row.names(ROSMAP_RIN_PMI) <- ROSMAP_md$projid
  ROSMAP_RIN_PMI <- ROSMAP_RIN_PMI[!is.nan(ROSMAP_RIN_PMI$RIN) & !is.nan(ROSMAP_RIN_PMI$PMI) & !is.na(ROSMAP_RIN_PMI$RIN) & !is.na(ROSMAP_RIN_PMI$PMI), ]
  
  
  df1 <- data.frame(normalized_ROSMAP_counts_t)
  df1$SampleID <- rownames(normalized_ROSMAP_counts_t)
  
  df2 <- data.frame(ROSMAP_RIN_PMI)
  df2$SampleID <- rownames(ROSMAP_RIN_PMI)
  ROSMAP_ready <- left_join(df2, df1, by = "SampleID")
  ROSMAP_ready <- ROSMAP_ready |>
    select(-SampleID)
  
  cor_results <- corr.test(
    x = ROSMAP_ready[, c("RIN", "PMI")],
    y = ROSMAP_ready[, !(names(ROSMAP_ready) %in% c("RIN", "PMI"))],
    method = "spearman",
    adjust = "fdr"
  )
  
  # Extract correlation matrix
  cor_matrix <- cor_results$r
  
  # Extract adjusted p-values
  p_values <- cor_results$p
  padj_values <- cor_results$p.adj
  
  # Combine the matrices into one data frame
  combined_df <- data.frame(
    PMI = cor_matrix["PMI", ],
    PMI_p = p_values["PMI", ],
    PMI_p_adj = padj_values["PMI", ],
    RIN = cor_matrix["RIN", ],
    RIN_p = p_values["RIN", ],
    RIN_p_adj = padj_values["RIN", ]
  )
  
  
  # Split the data frame into two separate tables
  PMI_table <- combined_df[, c("PMI", "PMI_p", "PMI_p_adj")]
  colnames(PMI_table) <- c("PMI", "p", "p.adj")
  file_path <- paste0("rna_quality_paper/correlation_on_each_gene (supp 1&2)/", name, "_gene_PMI_corr.csv")
  write.csv(PMI_table, file_path, row.names = TRUE)
  
  RIN_table <- combined_df[, c("RIN", "RIN_p", "RIN_p_adj")]
  colnames(RIN_table) <- c("RIN", "p", "p.adj")
  file_path <- paste0("rna_quality_paper/correlation_on_each_gene (supp 1&2)/", name, "_gene_RIN_corr.csv")
  write.csv(RIN_table, file_path, row.names = TRUE)
}

generate_tables(ROSMAP_md, ROSMAP_counts_oligo, "ROSMAP_oligo")
generate_tables(ROSMAP_md_other, ROSMAP_counts_riboDep, "ROSMAP_riboDep")
