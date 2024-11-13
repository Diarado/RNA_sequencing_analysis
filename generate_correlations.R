# new table for correlation among ROSMAP * 2 and NYBB in RIN, PMI, PTT 

# ------------------------------------------------------------------------------
# Author: Jiahe Tian, based on Sophie Ross's script
# Date: 7/11/2024
# Purpose: It generates Correlation Matrices of RIN, PMI, and other RNA quality measures on 
# NPH, NYBB, ROSMAP_polyA and ROSMAP_riboDepletion datasets
# ------------------------------------------------------------------------------
# Input:
#   - /rna_quality_paper/data
# Output:
#   - /rna_quality_paper/correlation_matrices (table 1&2)
# ------------------------------------------------------------------------------

# Load required packages
library(dplyr)
library(stringr)
library(tidyverse)
library(readxl)
library(psych)
library(openxlsx)
library(sva)
library(purrr)

# ------------------------------------------------------------------------------
# Data Loading
# ------------------------------------------------------------------------------

# please specify your working directory
setwd("C:/Users/jt3586/drs/")

#loading in the data
BB_md1 <- read.csv("rna_quality_paper/data/p1-A_NYBB_Teich(in).csv")

sheet_names <- excel_sheets("rna_quality_paper/data/Teich NYBB Sample Tracking.xlsx")
sheets <- lapply(sheet_names, function(sheet) read_excel("rna_quality_paper/data/Teich NYBB Sample Tracking.xlsx", sheet = sheet))
names(sheets) <- sheet_names

# Access individual sheets
sheet1 <- sheets[["QC_TEI_14989"]]
sheet2 <- sheets[["QC-DEJ_14869"]]
sheet3 <- sheets[["Initial QC_TEI_14060"]]
sheet4 <- sheets[["QC-TEI_13969_B01_EXS_RNA-92-567"]]
sheet5 <- sheets[["QC-TEI_14202 Testing"]]
sheet6 <- sheets[["QC-TEI_14266 KAPA Seq Metrics"]]


combined_sheets <- bind_rows(
  sheet1 %>% drop_na(DV200) %>% select(`NYGC Sample Name`, DV200) %>% mutate(DV200 = as.numeric(DV200)),
  sheet2 %>% drop_na(DV200) %>% select(`NYGC Sample Name`, DV200) %>% mutate(DV200 = as.numeric(DV200)),
  sheet3 %>% drop_na(DV200) %>% select(`NYGC Sample Name`, DV200) %>% mutate(DV200 = as.numeric(DV200)),
  sheet4 %>% drop_na(DV200) %>% select(`NYGC Sample Name`, DV200) %>% mutate(DV200 = as.numeric(DV200)),
  sheet5 %>% drop_na(DV200) %>% select(`NYGC Sample Name`, DV200) %>% mutate(DV200 = as.numeric(DV200)),
  sheet6 %>% drop_na(DV200) %>% select(`NYGC Sample Name`, DV200) %>% mutate(DV200 = as.numeric(DV200))
)

# Join combined_sheets to BB_md1
BB_md1 <- BB_md1 %>%
  left_join(combined_sheets, 
            by = c("NYGC.Sample.Name" = "NYGC Sample Name"))

# DV200 old and new that didn't match, if any
discrepancies <- list()

# Function to merge columns and check for discrepancies
merge_columns <- function(ba, dv, row_index) {
  if (!is.na(ba) & !is.na(dv)) {
    if (ba != dv) {
      discrepancies[[length(discrepancies) + 1]] <<- list(row = row_index, BA.DV200 = ba, DV200 = dv)
    }
    return(ba)
  } else if (!is.na(ba)) {
    return(ba)
  } else if (!is.na(dv)) {
    return(dv)
  } else {
    return(NA)
  }
}

# Apply the function to each row
BB_md1$BA.DV200 <- mapply(merge_columns, BB_md1$BA.DV200, BB_md1$DV200, seq_len(nrow(BB_md1)))
BB_md <- BB_md1
BB_counts <- read.csv("rna_quality_paper/data/p1_raw_counts_geneID.csv", row.names=1, check.names = F)
BB_pmi <- read.csv("rna_quality_paper/data/All_P-1_Demos_with_Rs.csv", check.names = F)

NPH_md <- read.csv("rna_quality_paper/data/NPH_metadata.csv")
NPH_medians <- read.csv("rna_quality_paper/data/median_nph_picard_500_t1000.csv") # use this for 
#NPH_medians <- read.csv("rna_quality_paper/data/median_nph_uncheck.csv") # use this for Picard

NPH_medians <- NPH_medians |>
  rename(SampleID = SAMPLE)
NPH_counts <- read.csv("rna_quality_paper/data/NPH_raw_counts.csv", row.names=1, check.names = F)

ROSMAP_md <-  read.csv("rna_quality_paper/data/ROSMAP_metadata_full_899_final.csv")
ROSMAP_medians <- read.csv("rna_quality_paper/data/median_rosmap_picard_500_t1000.csv") # use this for 
#ROSMAP_medians <- read.csv("rna_quality_paper/data/median_rosmap_match.csv") # use this for Picard

ROSMAP_medians <- ROSMAP_medians |>
  rename(specimenID = SAMPLE)
ROSMAP_counts <-  read.csv("rna_quality_paper/data/ROSMAP_DLPFC_final_raw_59653x899_projid_counts 1.csv", row.names=1, check.names = F)


################################################################################
# getting the correct formatting 
# so that the column names of the counts match the row names of the metadata

#BB
BB_counts <- BB_counts[, colnames(BB_counts)[match(BB_md$Client.Sample, colnames(BB_counts))]]
BB_pmi <- BB_pmi[BB_pmi$T %in% BB_md$Client.Sample,]
BB_pmi <- BB_pmi[match(BB_md$Client.Sample, BB_pmi$T),]

#NPH
NPH_medians$SampleID <- sapply(str_extract_all(NPH_medians$SampleID, "[0-9]+"), paste, collapse = "")
print(duplicated(NPH_medians$SampleID))
NPH_md$SampleID <- sapply(str_extract_all(NPH_md$SampleID, "[0-9]+"), paste, collapse = "")

NPH_tmp <- inner_join(NPH_medians, NPH_md, by = "SampleID")
NPH_md <- NPH_tmp

colnames(NPH_counts) <- gsub("[^0-9]", "", colnames(NPH_counts))
NPH_counts <- NPH_counts[, colnames(NPH_counts)[match(NPH_md$SampleID, colnames(NPH_counts))]]

#ROSMAP
ROSMAP_md$specimenID <- sub("^Sample_", "", ROSMAP_md$specimenID)
# this extracts the sample id that belongs to ROSMAP oligodT (aka. batch 1)

# Merge before splitting, ensuring all samples in ROSMAP_md are retained
ROSMAP_tmp <- merge(ROSMAP_medians, ROSMAP_md, by = "specimenID", all.y = TRUE)
ROSMAP_md <- ROSMAP_tmp

# Then split the merged data
samples_keep <- readLines("rna_quality_paper/data/batch1_samples.txt")
ROSMAP_md_1 <- ROSMAP_md[ROSMAP_md$specimenID %in% samples_keep, ]
ROSMAP_md_other <- ROSMAP_md[!(ROSMAP_md$specimenID %in% samples_keep), ]



###############################################################################

#getting the fraction of counts from the top 10 genes over the total expression
top_10_ratio <- function(gene_counts){
  # making a vector to hold the top 10 gene ratios
  top_10_ratios <- numeric(ncol(gene_counts))
  
  for (i in 1:ncol(gene_counts)) {
    # Get the top 10 genes for each sample
    top_10_genes <- head(sort(gene_counts[, i], decreasing = TRUE), 10)
    
    # Sum the counts of the top 10 genes
    top_10_sum <- sum(top_10_genes)
    
    # Calculate the total gene expression for the sample
    total_expression <- sum(gene_counts[, i])
    
    # Calculate the ratio of the top 10 genes to total gene expression
    top_10_ratio <- top_10_sum / total_expression
    
    # Store the ratio in the vector
    top_10_ratios[i] <- top_10_ratio
  }
  
  # Print or use the results as needed
  return(top_10_ratios)
}

BB_top_10 <- top_10_ratio(BB_counts)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05574 0.07767 0.09497 0.11496 0.13160 0.40524


NPH_top_10 <- top_10_ratio(NPH_counts)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1275  0.2356  0.2971  0.3085  0.3731  0.5601 

# filter_counts_ROSMAP:
# Filter ROSMAP_counts to keep columns where colnames are in ROSMAP_md$projid
# oligo dT
ROSMAP_counts_filtered <- ROSMAP_counts[, colnames(ROSMAP_counts) %in% ROSMAP_md_1$projid]
ROSMAP_top_10 <- top_10_ratio(ROSMAP_counts_filtered)

ROSMAP_counts_other <- ROSMAP_counts[, colnames(ROSMAP_counts) %in% ROSMAP_md_other$projid]
ROSMAP_top_10_2 <- top_10_ratio(ROSMAP_counts_other)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05712 0.09961 0.20803 0.23143 0.32084 0.84062 


# making dataframes for each individual dataset
#BB:
# convert time from HH:MM to mins
BB_pmi <- BB_pmi |>
  mutate(PMIRTFzn = str_split(PMIRTFzn, ":")) |>
  mutate(PMIRTFzn = map_dbl(PMIRTFzn, ~ {
    hours <- as.numeric(.x[1])
    minutes <- as.numeric(.x[2])
    seconds <- ifelse(length(.x) >= 3, as.numeric(.x[3]), 0)
    hours + (minutes / 60) + (seconds / 3600)
  }))

# generate data frames ready for the correlation analysis
BB_relevant <- data.frame(RIN = BB_md$BA.FQN, PMI = BB_pmi$PMIRTFzn,
                          PTT = BB_top_10, dv200 = BB_md$BA.DV200)

NPH_relevant <- data.frame(RIN = NPH_md$RIN, PTT = NPH_top_10,
                           median_5_bias = NPH_md$MEDIAN_5PRIME_BIAS, median_3_bias = NPH_md$MEDIAN_3PRIME_BIAS, 
                           median_5_3_bias = NPH_md$MEDIAN_5PRIME_TO_3PRIME_BIAS, Batch = NPH_md$Batch)


ROSMAP_oligo_relevant <- data.frame(RIN = ROSMAP_md_1$RIN, PMI = ROSMAP_md_1$pmi, PTT = ROSMAP_top_10,
                                    median_5_bias = ROSMAP_md_1$MEDIAN_5PRIME_BIAS, 
                                    median_3_bias = ROSMAP_md_1$MEDIAN_3PRIME_BIAS, median_5_3_bias = ROSMAP_md_1$MEDIAN_5PRIME_TO_3PRIME_BIAS
                                    , Batch = ROSMAP_md_1$sequencing_batch)

ROSMAP_riboDep_relevant <- data.frame(RIN = ROSMAP_md_other$RIN, PMI = ROSMAP_md_other$pmi, PTT = ROSMAP_top_10_2)

autopsy_relevant <- rbind(
  BB_relevant[, c("RIN", "PMI", "PTT")],
  ROSMAP_oligo_relevant[, c("RIN", "PMI", "PTT")],
  ROSMAP_riboDep_relevant[, c("RIN", "PMI", "PTT")]
)

# ------------------------------------------------------------------------------
# Generate Correlation Matrices
# ------------------------------------------------------------------------------


# Calculate the minimum RIN values from the NPH_relevant dataset
min_rin <- min(NPH_relevant$RIN, na.rm = TRUE)
print(min_rin)
# Filter the ROSMAP_relevant dataset to include only rows with RIN within the range [min_rin, inf]

filtered_BB_relevant <- BB_relevant[BB_relevant$RIN >= min_rin, ]
filtered_ROSMAP_oligo_relevant <- ROSMAP_oligo_relevant[ROSMAP_oligo_relevant$RIN >= min_rin, ]
filtered_ROSMAP_riboDep_relevant <- ROSMAP_riboDep_relevant[ROSMAP_riboDep_relevant$RIN >= min_rin, ]

filtered_autopsy_relevant <- autopsy_relevant[autopsy_relevant$RIN >= min_rin,]
# Replace any non-numeric values in PMI with NA

filtered_BB_relevant$PMI <- as.numeric(as.character(filtered_BB_relevant$PMI))
filtered_BB_relevant$PMI[is.nan(filtered_BB_relevant$PMI)] <- NA

filtered_ROSMAP_oligo_relevant$PMI <- as.numeric(as.character(filtered_ROSMAP_oligo_relevant$PMI))
filtered_ROSMAP_oligo_relevant$PMI[is.nan(filtered_ROSMAP_oligo_relevant$PMI)] <- NA

filtered_ROSMAP_riboDep_relevant$PMI <- as.numeric(as.character(filtered_ROSMAP_riboDep_relevant$PMI))
filtered_ROSMAP_riboDep_relevant$PMI[is.nan(filtered_ROSMAP_riboDep_relevant$PMI)] <- NA

filtered_autopsy_relevant$PMI <- as.numeric(as.character(filtered_autopsy_relevant$PMI))
filtered_autopsy_relevant$PMI[is.nan(filtered_autopsy_relevant$PMI)] <- NA

generate_correlation_tables <- function(data, name) {
  # Select all numeric columns for correlation
  numeric_data <- data[, sapply(data, is.numeric)]
  
  # Perform Spearman's correlations
  cor_results <- corr.test(numeric_data, numeric_data, method = "spearman", adjust = "fdr")
  
  # Adjust p-values
  p_adj <- cor_results$p.adj
  cor_matrix <- cor_results$r
  
  rownames(cor_matrix) <- colnames(numeric_data)
  colnames(cor_matrix) <- colnames(numeric_data)
  
  # Combine correlation coefficients and adjusted p-values into one table
  combined_matrix <- cor_matrix
  for (i in 1:nrow(p_adj)) {
    for (j in 1:ncol(p_adj)) {
      combined_matrix[i, j] <- ifelse(i > j, paste0(round(cor_matrix[i, j], 3), " (p = ", round(p_adj[i, j], 3), ")"), "")
    }
  }
  
  # Convert to data frame
  combined_df <- as.data.frame(combined_matrix)
  
  return(list(name = name, data = combined_df))
}

# Create a list to hold all the sheets
sheet_list <- list()

# Generate all correlation tables and add them to the sheet list
sheet_list[[1]] <- generate_correlation_tables(NPH_relevant, "NPH")
sheet_list[[2]] <- generate_correlation_tables(BB_relevant, "NYBB")
sheet_list[[3]] <- generate_correlation_tables(filtered_BB_relevant, "NYBB_after_filter")
sheet_list[[4]] <- generate_correlation_tables(ROSMAP_oligo_relevant, "ROSMAP_oligo")
sheet_list[[5]] <- generate_correlation_tables(filtered_ROSMAP_oligo_relevant, "ROSMAP_oligo_after_filter")
sheet_list[[6]] <- generate_correlation_tables(ROSMAP_riboDep_relevant, "ROSMAP_riboDep")
sheet_list[[7]] <- generate_correlation_tables(filtered_ROSMAP_riboDep_relevant, "ROSMAP_riboDep_after_filter")
sheet_list[[8]] <- generate_correlation_tables(autopsy_relevant, "autopsy_combined")
sheet_list[[9]] <- generate_correlation_tables(filtered_autopsy_relevant, "autopsy_combined_after_filter")

# a new Excel workbook
wb <- createWorkbook()
for (sheet in sheet_list) {
  addWorksheet(wb, sheet$name)
  writeData(wb, sheet$name, sheet$data, rowNames = TRUE)
}

saveWorkbook(wb, "rna_quality_paper/correlation_matrices (table 1&2)/Correlations_median_500_t1000.xlsx", overwrite = TRUE)

# save the varibles created for individual gene correlation analysis
save.image(file = "rna_quality_paper/correlation_on_each_gene (supp 1&2)/processed_data.RData")
