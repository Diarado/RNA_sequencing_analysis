# Load the necessary packages
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(dplyr)

# -----------------------------------------------------------------------------------
# Define a Sanitization Function for Directory Names
# -----------------------------------------------------------------------------------
sanitize_dir_name <- function(name) {
  # Define a pattern of invalid characters for Windows
  invalid_chars <- "[<>:\"/\\\\|?*()]+"
  # Replace invalid characters with an underscore
  clean_name <- gsub(invalid_chars, "_", name)
  # Optionally, replace spaces with underscores for consistency
  clean_name <- gsub(" ", "_", clean_name)
  return(clean_name)
}

# -----------------------------------------------------------------------------------
# Define Variable Label Mapping
# -----------------------------------------------------------------------------------
var_label_map <- c(
  "median_5_bias" = "5' Bias",
  "median_3_bias" = "3' Bias",
  "median_5_3_bias" = "5' to 3' Ratio",
  "dv200" = "DV200"
)

# -----------------------------------------------------------------------------------
# 3. Define a Function to Create Scatter Plots with Conditional Batch Coloring
# -----------------------------------------------------------------------------------
create_scatter_plot <- function(data, x_var, y_var, plot_title, x_label = x_var, y_label = y_var) {
  # Ensure the variables are numeric
  data[[x_var]] <- as.numeric(data[[x_var]])
  data[[y_var]] <- as.numeric(data[[y_var]])
  
  # Check if 'Batch' exists in the dataset
  has_batch <- "Batch" %in% colnames(data)
  
  if (has_batch) {
    # Remove rows with NA in x_var, y_var, or Batch
    plot_data <- data %>% select(all_of(c(x_var, y_var, "Batch"))) %>% na.omit()
    
    # Check if there are enough data points to perform correlation
    if (nrow(plot_data) < 3) {
      warning(paste("Not enough data points for", plot_title))
      return(NULL)
    }
    
    # Calculate Spearman correlation (excluding Batch)
    cor_test <- suppressWarnings(cor.test(plot_data[[x_var]], plot_data[[y_var]], method = "spearman"))
    r_value <- round(cor_test$estimate, 2)
    p_value <- signif(cor_test$p.value, 2)
    
    # Create the scatter plot with Batch coloring
    p <- ggplot(plot_data, aes_string(x = x_var, y = y_var, color = "Batch")) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      theme_minimal() +
      labs(title = plot_title,
           x = x_label,
           y = y_label) +
      annotate("text", 
               x = Inf, y = -Inf, 
               label = paste0("r = ", r_value, "\np = ", p_value),
               hjust = 1.1, vjust = -0.1,
               size = 4, color = "black")
    
  } else {
    # Remove rows with NA in x_var and y_var
    plot_data <- data %>% select(all_of(c(x_var, y_var))) %>% na.omit()
    
    # Check if there are enough data points to perform correlation
    if (nrow(plot_data) < 3) {
      warning(paste("Not enough data points for", plot_title))
      return(NULL)
    }
    
    # Calculate Spearman correlation
    cor_test <- suppressWarnings(cor.test(plot_data[[x_var]], plot_data[[y_var]], method = "spearman"))
    r_value <- round(cor_test$estimate, 2)
    p_value <- signif(cor_test$p.value, 2)
    
    # Create the scatter plot without Batch coloring
    p <- ggplot(plot_data, aes_string(x = x_var, y = y_var)) +
      geom_point(color = "blue", alpha = 0.6) +
      geom_smooth(method = "lm", color = "red", se = FALSE) +
      theme_minimal() +
      labs(title = plot_title,
           x = x_label,
           y = y_label) +
      annotate("text", 
               x = Inf, y = -Inf, 
               label = paste0("r = ", r_value, "\np = ", p_value),
               hjust = 1.1, vjust = -0.1,
               size = 4, color = "black")
  }
  
  return(p)
}

# -----------------------------------------------------------------------------------
# 4. Define a Function to Generate and Save Scatter Plots for All Variable Pairs
# -----------------------------------------------------------------------------------
generate_and_save_plots <- function(data, dataset_name, output_dir = "Scatter_Plots") {
  # Sanitize the dataset name for directory creation
  clean_dataset_name <- sanitize_dir_name(dataset_name)
  
  # Create a subdirectory for the dataset
  dataset_dir <- file.path(output_dir, clean_dataset_name)
  if (!dir.exists(dataset_dir)) {
    dir.create(dataset_dir, recursive = TRUE)
  }
  
  
  # Select only numeric columns (excluding 'Batch' if present)
  numeric_vars <- data %>% 
    select_if(is.numeric) %>% 
    colnames()
  
  # Generate all unique pairs of numeric variables
  var_pairs <- combn(numeric_vars, 2, simplify = FALSE)
  
  # Initialize a counter for progress tracking
  total_plots <- length(var_pairs)
  plot_counter <- 1
  
  # Iterate through each pair and generate plots
  for (pair in var_pairs) {
    x_var <- pair[1]
    y_var <- pair[2]
    
    # Map variable names to labels if they exist in var_label_map
    x_label <- ifelse(x_var %in% names(var_label_map), var_label_map[x_var], x_var)
    y_label <- ifelse(y_var %in% names(var_label_map), var_label_map[y_var], y_var)
    
    # Construct the plot title using the mapped labels
    plot_title <- paste(dataset_name, ":", x_label, "vs", y_label)
    
    # Create the scatter plot with custom labels
    plot <- create_scatter_plot(data, x_var, y_var, plot_title, x_label, y_label)
    
    # Check if plot was created successfully
    if (!is.null(plot)) {
      # Define the filename (use original variable names for filenames)
      filename <- paste0(x_var, "_vs_", y_var, ".png")
      filepath <- file.path(dataset_dir, filename)
      
      # Save the plot as PNG
      ggsave(filename = filepath, plot = plot, width = 6, height = 4, dpi = 300)
      
      cat(sprintf("Saved plot %d/%d: %s\n", plot_counter, total_plots, filepath))
    } else {
      cat(sprintf("Skipped plot %d/%d: %s vs %s (Insufficient data)\n", 
                  plot_counter, total_plots, x_var, y_var))
    }
    
    plot_counter <- plot_counter + 1
  }
}

# -----------------------------------------------------------------------------------
# 5. Define the List of Datasets (Including Filtered and Unfiltered)
# -----------------------------------------------------------------------------------

# Define the list of datasets with appropriate names
datasets <- list(
  `NYBB Autopsy` = BB_relevant,
  `NYBB Autopsy (RIN > 6)` = filtered_BB_relevant,
  `NPH Biopsy` = NPH_relevant,
  `ROSMAP poly-A Selection` = ROSMAP_oligo_relevant,
  `ROSMAP poly-A Selection (RIN > 6)` = filtered_ROSMAP_oligo_relevant,
  `ROSMAP Ribosomal Depletion` = ROSMAP_riboDep_relevant,
  `ROSMAP Ribosomal Depletion (RIN > 6)` = filtered_ROSMAP_riboDep_relevant
)

# -----------------------------------------------------------------------------------
# 6. Generate and Save Scatter Plots for Each Dataset
# -----------------------------------------------------------------------------------

# Define the main output directory for scatter plots
main_output_dir <- "C:/Users/jt3586/drs/rna_quality_paper/supp_figures/Scatter_Plots"

# Iterate through each dataset and generate plots
for (dataset_name in names(datasets)) {
  cat(sprintf("Processing dataset: %s\n", dataset_name))
  generate_and_save_plots(datasets[[dataset_name]], dataset_name, main_output_dir)
  cat(sprintf("Completed dataset: %s\n\n", dataset_name))
}
# Load necessary libraries
# Prepare the data
# Prepare the data
# Prepare the data
data <- data.frame(PMI = autopsy_relevant$PMI)

# Create a layout for two plots side by side and increase the overall size of the plot
par(mfrow = c(1, 2), mar = c(5, 4, 2, 2) + 0.1, oma = c(0, 0, 4, 0), cex.main = 0.9)

# Draw the boxplot on the left with a color
boxplot(data$PMI, 
        col = "lightgreen", 
        ylab = "PMI (hour)")

# Draw the histogram on the right with more breaks and no individual title
hist(data$PMI, 
     xlab = "PMI (hour)", 
     main = "",  # Ensures no title is added to the histogram
     col = "lightblue", 
     border = "black", 
     breaks = 20)

# Add a centered, two-line grand title for both plots
mtext("Distribution of PMI Across Autopsy Samples (n = 1412)\n(NYBB, ROSMAP poly-A Selection, and ROSMAP Ribosomal Depletion)", 
      outer = TRUE, cex = 1.2, line = 1)

# Reset plotting layout
par(mfrow = c(1, 1))
