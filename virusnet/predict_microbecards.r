#!/usr/bin/env Rscript

# Suppress TensorFlow messages
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)

# Load required libraries
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
  library(progress)
}))

source("utils.r")

# Define custom objects for the model (adjust as needed for your model)
custom_objects <- list(
  "layer_pos_embedding" = layer_pos_embedding,
  "layer_pos_sinusoid" = layer_pos_sinusoid,
  "layer_transformer_block" = layer_transformer_block,
  "layer_aggregate_time_dist" = layer_aggregate_time_dist
)

# Utility function for logging
log_info <- function(message) {
  cat("[INFO]", message, "\n")
}

process_fasta <- function(fasta_file, model, output_csv, window_size, step) {
  library(dplyr)
  
  # Read all contigs from the FASTA file
  fasta_data_all <- readFasta(fasta_file)
  if (nrow(fasta_data_all) == 0) stop("FASTA file is empty.")
  
  # Add an index column to track original positions
  fasta_data_all$index <- seq_len(nrow(fasta_data_all))
  fasta_data_all$Length <- nchar(fasta_data_all$Sequence)
  fasta_data_all$SanitizedHeader <- make.names(fasta_data_all$Header, unique = TRUE)
  
  # Filter contigs with Length >= window_size
  fasta_data <- fasta_data_all %>% filter(Length >= window_size)
  
  # Define output CSV file path
  base_name <- tools::file_path_sans_ext(basename(fasta_file))
  csv_file <- file.path(output_csv, paste0(base_name, ".csv"))
  
  # Create a temporary directory for H5 files
  temp_output_dir <- file.path(tempdir(), base_name)
  dir.create(temp_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Perform prediction for the entire FASTA file
  predict_model(
    output_format = "by_entry",
    model = model,
    layer_name = "dense_26",
    path_input = fasta_file,
    step = step,
    batch_size = 1,
    verbose = FALSE,
    padding = "self",
    mode = "label",
    format = "fasta",
    output_dir = temp_output_dir,
    return_int = TRUE
  )
  
  # Initialize list for summary statistics
  summary_list <- list()
  
  # Initialize counters for virus and non-virus regions
  virus_regions <- 0
  non_virus_regions <- 0
  total_regions <- 0
  
  # Initialize counters for nucleotides
  viral_nt <- 0
  non_viral_nt <- 0
  
  # Process only the filtered contigs
  for (i in 1:nrow(fasta_data)) {
    # Get the original index of the contig
    idx <- fasta_data$index[i]
    contig <- fasta_data$SanitizedHeader[i]
    contig_length <- fasta_data$Length[i]
    
    # Define the H5 file path
    h5_file_path <- file.path(temp_output_dir, paste0("states_nr_", idx, ".h5"))
    if (!file.exists(h5_file_path)) {
      warning(paste("H5 file not found for contig", contig, "at index", idx))
      next
    }
    
    # Load predictions for this contig
    preds_contig <- load_prediction(h5_path = h5_file_path, get_sample_position = TRUE, verbose = FALSE)
    
    # Extract states and sample_end_position
    states <- as.data.frame(preds_contig$states)
    colnames(states) <- c("non_virus", "virus")
    virus_probs <- states$virus
    num_preds <- length(preds_contig$sample_end_position)
    
    # Calculate summary statistics
    mean_prob <- mean(virus_probs, na.rm = TRUE)
    median_prob <- median(virus_probs, na.rm = TRUE)
    sd_prob <- sd(virus_probs, na.rm = TRUE)
    # Calculate IQR
    iqr_prob <- IQR(virus_probs, na.rm = TRUE)
    # Determine classification
    classification <- ifelse(mean_prob > 0.5, "virus", "non-virus")
    
    # Update virus/non-virus counters
    if (mean_prob > 0.5) {
      virus_regions <- virus_regions + 1
      viral_nt <- viral_nt + contig_length
    } else {
      non_virus_regions <- non_virus_regions + 1
      non_viral_nt <- non_viral_nt + contig_length
    }
    total_regions <- total_regions + 1
    
    # Store summary stats
    summary_list[[i]] <- data.frame(
      contig = contig,
      mean_prediction = mean_prob,
      median_prediction = median_prob,
      sd_prediction = sd_prob,
      iqr_prediction = iqr_prob,
      classification = classification,
      num_predictions = num_preds,
      contig_length = contig_length
    )
  }
  
  # Combine results and write to CSV with semicolon delimiter and no quotes
  summary_df <- bind_rows(summary_list)
  write.table(summary_df, csv_file, sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Delete temporary H5 files
  unlink(temp_output_dir, recursive = TRUE)
  
  message("Processed:", fasta_file)
  
  # Return summary stats for the whole file
  return(list(
    file_name = base_name,
    virus_percentage = ifelse(total_regions > 0, virus_regions / total_regions * 100, 0),
    non_virus_percentage = ifelse(total_regions > 0, non_virus_regions / total_regions * 100, 0),
    viral_nt = viral_nt,
    non_viral_nt = non_viral_nt
  ))
}

# Main function
main <- function() {
  # Hardcoded parameters (modify as needed)
  input_folder <- "../test"  # Replace with your input folder path
  model_binary <- "~/.virusnet/transfer_learning_virus_bert_5_85.h5"  # Replace with your model path
  output_csv <- "output_csv5"
  window_size <- 5000
  step <- 5000  # Step size equals window size for coarse predictions
  
  # Create output directory for CSV files
  dir.create(output_csv, showWarnings = FALSE, recursive = TRUE)
  
  # List all FASTA files in the input folder
  fasta_files <- list.files(input_folder, pattern = "\\.fasta$", full.names = TRUE)
  if (length(fasta_files) == 0) stop("No FASTA files found in the input folder.")
  
  # Load the model
  model <- keras::load_model_hdf5(model_binary, custom_objects = custom_objects)
  
  # Initialize progress bar
  pb <- progress_bar$new(total = length(fasta_files), format = "[:bar] :percent :eta")
  
  # Initialize list for summary output
  summary_output_list <- list()
  
  # Process each FASTA file
  for (i in 1:length(fasta_files)) {
    fasta_file <- fasta_files[i]
    summary_stats <- process_fasta(fasta_file, model, output_csv, window_size, step)
    
    # Add to summary list
    summary_output_list[[i]] <- data.frame(
      file_name = summary_stats$file_name,
      virus_percentage = summary_stats$virus_percentage,
      non_virus_percentage = summary_stats$non_virus_percentage,
      viral_nt = summary_stats$viral_nt,
      non_viral_nt = summary_stats$non_viral_nt
    )
    
    pb$tick()
  }
  
  # Combine and write summary output
  summary_output_df <- bind_rows(summary_output_list)
  summary_output_path <- file.path(output_csv, "summary_output.csv")
  write.table(summary_output_df, summary_output_path, sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  log_info("Processing complete.")
}

# Run the main function
main()