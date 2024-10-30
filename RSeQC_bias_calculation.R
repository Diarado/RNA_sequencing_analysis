#!/usr/bin/env Rscript
library(dplyr)
library(stringr)

# Function to read and process each file
process_file <- function(filepath) {
  
  lines <- readLines(filepath)
  
  # Check that there are at least two lines (header and data)
  if (length(lines) < 2) {
    message(paste("Skipping file", filepath, "- does not contain enough lines."))
    return(NULL)
  }
  
  # Extract the header line
  header_line <- lines[1]
  header_values <- unlist(strsplit(header_line, "\t"))
  
  # Confirm that the header starts with "Percentile"
  if (header_values[1] != "Percentile") {
    message(paste("Skipping file", filepath, "- does not have 'Percentile' as the first column in the header."))
    return(NULL)
  }
  
  # Extract bin numbers from the header (should be 1 to 100)
  bin_numbers <- as.numeric(header_values[-1])
  
  # Check that bin_numbers are from 1 to 100
  if (!all(bin_numbers == 1:100)) {
    message(paste("Skipping file", filepath, "- bins are not numbered from 1 to 100."))
    return(NULL)
  }
  
  # Extract the data line (second line)
  data_line <- lines[2]
  
  # Split the data line into values
  data_values <- unlist(strsplit(data_line, "\t"))
  
  # Ensure that the number of data values matches the header (sample name + 100 bins)
  if (length(data_values) != length(header_values)) {
    message(paste("Skipping file", filepath, "- mismatch between header and data columns."))
    return(NULL)
  }
  
  # Extract the sample name using the desired method
  sample_name <- sub("\\..*", "", basename(filepath))
  
  # Extract coverage values and convert to numeric
  coverage_values <- as.numeric(data_values[-1])
  
  # Check for NA values in coverage (non-numeric entries)
  if (any(is.na(coverage_values))) {
    message(paste("Skipping file", filepath, "- contains non-numeric coverage values."))
    return(NULL)
  }
  
  # Define regions (adjust the percentage as needed)
  percent_coverage <- 0.1  # Represents 10% of the gene body
  
  num_bins <- length(coverage_values)  # Should be 100
  print(num_bins)
  # Calculate indices for 5' and 3' regions
  bins_5_prime <- coverage_values[1:floor(num_bins * percent_coverage)]  # First 10%
  print(bins_5_prime)
  bins_3_prime <- coverage_values[(num_bins - floor(num_bins * percent_coverage) + 1):num_bins]  # Last 10%
  print(bins_3_prime)
  # Compute median coverage in each region and overall
  median_5_prime <- median(bins_5_prime, na.rm = TRUE)
  median_3_prime <- median(bins_3_prime, na.rm = TRUE)
  median_overall <- median(coverage_values, na.rm = TRUE)
  
  # Compute biases as ratios of median coverages
  median_5_prime_bias <- median_5_prime / median_overall
  median_3_prime_bias <- median_3_prime / median_overall
  median_5_prime_to_3_prime_bias <- median_5_prime / median_3_prime
  
  # Create a data frame with the required values
  tb <- data.frame(
    SAMPLE = sample_name,
    MEDIAN_5PRIME_BIAS = median_5_prime_bias,
    MEDIAN_3PRIME_BIAS = median_3_prime_bias,
    MEDIAN_5PRIME_TO_3PRIME_BIAS = median_5_prime_to_3_prime_bias,
    stringsAsFactors = FALSE
  )
  
  return(tb)
}

# Directories containing the files to process
#directories <- c("/mnt/data/jiahe/drs/NPH2/")
directories <- c("C:/Users/jt3586/drs/rna_quality/test")

# Initialize an empty data frame to store results
tbs <- data.frame(
  SAMPLE = character(),
  MEDIAN_5PRIME_BIAS = numeric(),
  MEDIAN_3PRIME_BIAS = numeric(),
  MEDIAN_5PRIME_TO_3PRIME_BIAS = numeric(),
  stringsAsFactors = FALSE
)

# Process each file in the directories
for (dir in directories) {
  files <- list.files(dir, pattern = "*.txt", full.names = TRUE)
  
  for (file in files) {
    tb <- process_file(file)
    if (!is.null(tb)) {
      tbs <- bind_rows(tbs, tb)
      message(paste("Processed file", file))
    }
    # If tb is NULL, the file was skipped
  }
}

# Finalize the results
final_tb <- tbs %>%
  distinct() %>%
  arrange(SAMPLE)

# Write the results to a CSV file
#write.csv(final_tb, file = "/mnt/data/jiahe/drs/median_nph_RSeQC.csv", row.names = FALSE)
write.csv(final_tb, file = "C:/Users/jt3586/drs/rna_quality/test/median_nph_RSeQC.csv", row.names = FALSE)

