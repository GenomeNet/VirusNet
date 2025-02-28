#!/usr/bin/env Rscript

# Suppress TF messages
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)

# Load libraries
suppressWarnings(suppressPackageStartupMessages({
  library(deepG)
  library(magrittr)
  library(microseq)
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(zoo)
  library(keras)
  library(reticulate)
  library(futile.logger)
}))

# Get the base directory of the Conda environment
conda_prefix <- Sys.getenv("CONDA_PREFIX")

# Source the R scripts
invisible(source(file.path(conda_prefix, "bin", "utils.r")))
invisible(source(file.path(conda_prefix, "bin", "setup_logger.r")))

# Function to safely write FASTA if not empty
safe_write_fasta <- function(fasta_subset, file_path, prediction_df = NULL) {
  if (!is.null(fasta_subset) && nrow(fasta_subset) > 0) {
    # If prediction_df is provided, add prediction values to headers
    if (!is.null(prediction_df)) {
      # Create a mapping of contig names to prediction values
      pred_map <- setNames(prediction_df$virus, prediction_df$contig_name)
      
      # Update headers with prediction values
      for (i in 1:nrow(fasta_subset)) {
        header <- fasta_subset$Header[i]
        # Find the prediction value for this contig
        if (header %in% names(pred_map)) {
          pred_value <- pred_map[header]
          # Format with 2 decimal places
          pred_formatted <- sprintf("%.2f", pred_value)
          # Update the header
          fasta_subset$Header[i] <- paste0(header, " | virus_prop=", pred_formatted)
        }
      }
    }
    
    writeFasta(fasta_subset, file_path)
    custom_log("INFO", paste("FASTA data written to:", file_path), "3/8")
  } else {
    custom_log("WARN", "No contigs found", "3/8")
  }
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
command_args <- list()
for (i in seq(1, length(args), 2)) {
  command_args[[gsub("--", "", args[i])]] <- args[i + 1]
}

# Set default virus threshold if not provided
if (is.null(command_args$virus_threshold)) {
  command_args$virus_threshold <- 0.5
} else {
  command_args$virus_threshold <- as.numeric(command_args$virus_threshold)
}

# Check if step_size is not the default value (1000) and warn that it's ignored in metagenome mode
if (!is.null(command_args$step_size) && as.numeric(command_args$step_size) != 1000) {
  custom_log("WARN", paste("Step size (", command_args$step_size, ") is ignored in metagenome mode. Metagenome mode always uses one prediction per contig."), "1/8")
}

# Check if the output directory exists, if not, create it
output_directory <- command_args$output
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

custom_log("INFO", "Checking input file", "1/8")
fasta_data <- readFasta(command_args$input)
num_fasta_entries <- nrow(fasta_data)

custom_log("INFO", paste("Number of FASTA entries in the file:", num_fasta_entries), "1/8")

if (num_fasta_entries == 0) {
  custom_log("ERROR", "Input FASTA file is empty. Please provide a non-empty FASTA file.", "1/8")
  stop("Empty input file.")
}

custom_log("INFO", "Loading binary model", "2/8")

suppressWarnings({
  model_binary <- keras::load_model_hdf5(command_args$model_binary, custom_objects = custom_objects)
})

custom_log("INFO", "Performing predictions", "3/8")
temp_file_binary <- tempfile(fileext = ".h5")

custom_log("INFO", "Using metagenomic mode", "3/8")
result <- tryCatch({
  sink("/dev/null")
  pred <- predict_model(
    output_format = "one_pred_per_entry",
    model = model_binary,
    layer_name = "dense_26",
    path_input = command_args$input,
    round_digits = 3,
    step = as.numeric(command_args$step_size),
    batch_size = as.numeric(command_args$batch_size),
    verbose = TRUE,
    return_states = FALSE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = temp_file_binary,
    return_int = length(model_binary$input_shape) == 2
  )
  sink()
  success <<- TRUE
  custom_log("INFO", "Successful prediction", "3/8")
}, error = function(e) {
  custom_log("WARN", paste("Error during prediction: ", e$message), "3/8")
})

if (!success) {
  custom_log("ERROR", "Failed to predict. Please check your input data and model.", "3/8")
  stop("Prediction failed.")
}

# Interpret the output file
prediction_df <- load_prediction(h5_path = temp_file_binary, get_sample_position = FALSE, verbose = FALSE)
prediction_df <- as.data.frame(prediction_df$states)

colnames(prediction_df) <- c("non_virus", "virus")
prediction_df$contig_name <- fasta_data$Header

# Make sure virus column has proper values
# Output summary using the logging system
virus_threshold <- as.numeric(command_args$virus_threshold)
viral_count <- length(which(prediction_df$virus >= virus_threshold))
non_viral_count <- length(which(prediction_df$virus < virus_threshold))

# Ensure virus column values are correct and not empty
cat("Debugging virus column values:\n")
cat(sprintf("Range: min=%f, max=%f\n", min(prediction_df$virus), max(prediction_df$virus)))
cat(sprintf("Using virus threshold: %f\n", virus_threshold))

custom_log("INFO", paste("Number of contigs classified as viral:", viral_count), "3/8")
custom_log("INFO", paste("Number of contigs classified as non-viral:", non_viral_count), "3/8")

# Add classification column
prediction_df$classification <- ifelse(prediction_df$virus >= virus_threshold, "virus", "non-virus")

# Make sure contig names don't contain commas which could break CSV parsing
prediction_df$contig_name <- gsub(",", ";", prediction_df$contig_name)

# Print debug information about contig names
cat("Debug - First few contig names in prediction_df:\n")
print(head(prediction_df$contig_name))

# Write CSV with strict formatting to ensure proper reading by pandas
non_summarized_output_path <- file.path(output_directory, "binary_results.csv")
write.csv(prediction_df, non_summarized_output_path, row.names = FALSE, quote = TRUE)

# Generate summarized results with standard deviation
# Group by contig_name and calculate statistics
if (nrow(prediction_df) > 0) {
  # For metagenome mode, we don't have multiple entries per contig
  # Create a simplified summarized file with only relevant fields
  contig_summary <- prediction_df %>%
    mutate(
      virus_probability = virus,
      is_virus = virus >= virus_threshold
    ) %>%
    select(contig_name, virus_probability, is_virus)
  
  # Write summarized results
  summarized_output_path <- file.path(output_directory, "binary_results_summarized.csv")
  write.csv(contig_summary, summarized_output_path, row.names = FALSE, quote = TRUE)
  
  custom_log("INFO", paste("Summarized results written to:", summarized_output_path), "3/8")
  custom_log("INFO", "Note: Metagenome mode uses one prediction per contig, so no standard deviation is calculated", "3/8")
}

# Print the first few rows of the CSV to verify format
cat("Verifying CSV output format:\n")
write.table(head(prediction_df, 5), quote = TRUE, row.names = FALSE, sep = ",")

# Subset and write FASTA data for viral and non-viral contigs
viral_contigs <- prediction_df$contig_name[prediction_df$virus >= virus_threshold]
non_viral_contigs <- prediction_df$contig_name[prediction_df$virus < virus_threshold]

# Debug viral contigs identification
cat("Debug - Number of viral contigs identified:", length(viral_contigs), "\n")
if (length(viral_contigs) > 0) {
  cat("Debug - First few viral contigs:\n")
  print(head(viral_contigs))
  
  # Extract contig IDs for debugging
  viral_contig_ids <- sapply(viral_contigs, function(x) strsplit(x, " ")[[1]][1])
  cat("Debug - First few viral contig IDs (first part before space):\n")
  print(head(viral_contig_ids))
}

# Debug FASTA headers
cat("Debug - First few FASTA headers:\n")
print(head(fasta_data$Header))

# Check if headers match exactly
if (length(viral_contigs) > 0) {
  # Replace semicolons back to commas for matching
  matching_headers <- gsub(";", ",", viral_contigs)
  
  # Check for exact matches
  exact_matches <- sum(matching_headers %in% fasta_data$Header)
  cat("Debug - Exact matches between viral contigs and FASTA headers:", exact_matches, "\n")
  
  # If no exact matches, try partial matching
  if (exact_matches == 0) {
    cat("Debug - Attempting partial matching...\n")
    # Extract just the first part of the contig name (before any spaces)
    short_viral_contigs <- gsub(" .*$", "", matching_headers)
    short_fasta_headers <- gsub(" .*$", "", fasta_data$Header)
    
    partial_matches <- sum(short_viral_contigs %in% short_fasta_headers)
    cat("Debug - Partial matches using first part of name:", partial_matches, "\n")
    
    # Use partial matching if needed
    if (partial_matches > 0) {
      viral_fasta_indices <- which(short_fasta_headers %in% short_viral_contigs)
      viral_fasta <- fasta_data[viral_fasta_indices, ]
      cat("Debug - Using partial matching for FASTA extraction\n")
    } else {
      viral_fasta <- NULL
    }
  } else {
    # Use exact matching
    viral_fasta <- fasta_data[fasta_data$Header %in% matching_headers, ]
  }
} else {
  viral_fasta <- NULL
}

# Same for non-viral contigs
if (length(non_viral_contigs) > 0) {
  # Replace semicolons back to commas for matching
  matching_headers <- gsub(";", ",", non_viral_contigs)
  
  # Check for exact matches
  exact_matches <- sum(matching_headers %in% fasta_data$Header)
  
  # If no exact matches, try partial matching
  if (exact_matches == 0) {
    # Extract just the first part of the contig name (before any spaces)
    short_non_viral_contigs <- gsub(" .*$", "", matching_headers)
    short_fasta_headers <- gsub(" .*$", "", fasta_data$Header)
    
    # Use partial matching if needed
    if (sum(short_non_viral_contigs %in% short_fasta_headers) > 0) {
      non_viral_fasta_indices <- which(short_fasta_headers %in% short_non_viral_contigs)
      non_viral_fasta <- fasta_data[non_viral_fasta_indices, ]
    } else {
      non_viral_fasta <- NULL
    }
  } else {
    # Use exact matching
    non_viral_fasta <- fasta_data[fasta_data$Header %in% matching_headers, ]
  }
} else {
  non_viral_fasta <- NULL
}

# Check if we found any viral sequences
if (!is.null(viral_fasta) && nrow(viral_fasta) > 0) {
  cat("Debug - Found", nrow(viral_fasta), "viral sequences to write to FASTA\n")
} else {
  cat("Debug - No viral sequences found to write to FASTA\n")
}

# Pass prediction_df to safe_write_fasta to add prediction values to headers
safe_write_fasta(viral_fasta, file.path(output_directory, "viral_contigs.fasta"), prediction_df)
safe_write_fasta(non_viral_fasta, file.path(output_directory, "non_viral_contigs.fasta"), prediction_df)

on.exit({
  unlink(temp_file_binary)
}, add = TRUE)