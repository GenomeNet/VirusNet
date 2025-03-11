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
  library(keras)
  library(reticulate)
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

# Function to process a single FASTA file
process_fasta <- function(fasta_file, model, model_genus, genus_labels, output_csv, window_size, step) {
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
  
  # Create a temporary directory for binary prediction H5 files
  temp_output_dir <- file.path(tempdir(), base_name)
  dir.create(temp_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Perform binary prediction for the entire FASTA file
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
    contig <- fasta_data$Header[i]
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
    iqr_prob <- IQR(virus_probs, na.rm = TRUE)
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
    
    # Store summary stats with original index
    summary_list[[i]] <- data.frame(
      original_index = fasta_data$index[i],
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
  
  # Combine results into a data frame
  summary_df <- bind_rows(summary_list)
  
  # Write viral contigs to a new FASTA file
  viral_indices <- summary_df %>% filter(classification == "virus") %>% pull(original_index)
  viral_contigs <- fasta_data_all %>% filter(index %in% viral_indices)
  viral_fasta_file <- file.path(output_csv, paste0(base_name, "_viral.fasta"))
  
  # Process genus prediction if there are viral contigs
  if (nrow(viral_contigs) > 0) {
    # Write viral contigs to FASTA file using original headers
    writeFasta(viral_contigs[, c("Header", "Sequence")], viral_fasta_file)
    
    # Create a temporary directory for genus prediction H5 files
    temp_genus_output_dir <- file.path(tempdir(), paste0(base_name, "_genus"))
    dir.create(temp_genus_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Run genus prediction on the viral FASTA file
    predict_model(
      output_format = "by_entry",
      model = model_genus,
      layer_name = "dense_3",
      path_input = viral_fasta_file,
      step = 500,
      batch_size = 10,
      verbose = FALSE,
      padding = "self",
      mode = "label",
      format = "fasta",
      output_dir = temp_genus_output_dir,
      return_int = length(model_genus$input_shape) == 2
    )
    
    # Initialize list for genus predictions
    genus_list <- list()
    
    # Process each viral contig for genus prediction
    for (j in 1:nrow(viral_contigs)) {
      h5_file_path <- file.path(temp_genus_output_dir, paste0("states_nr_", j, ".h5"))
      if (!file.exists(h5_file_path)) {
        warning(paste("Genus H5 file not found for viral contig", j))
        next
      }
      
      # Load genus predictions
      preds_genus <- load_prediction(h5_path = h5_file_path, get_sample_position = FALSE, verbose = FALSE)
      states <- preds_genus$states
      if (is.null(states)) {
        warning(paste("No states found in H5 file for viral contig", j))
        next
      }
      
      # Compute average probabilities across all positions for this contig
      avg_probs <- colMeans(states)
      max_idx <- which.max(avg_probs)
      predicted_genus <- genus_labels[max_idx]
      genus_probability <- avg_probs[max_idx]
      
      # Map back to original index
      genus_list[[j]] <- data.frame(
        original_index = viral_contigs$index[j],
        genus = predicted_genus,
        genus_probability = genus_probability
      )
    }
    
    # Combine genus predictions
    genus_df <- bind_rows(genus_list)
    
    # Clean up temporary genus H5 files
    unlink(temp_genus_output_dir, recursive = TRUE)
  } else {
    # If no viral contigs, create an empty genus data frame
    genus_df <- data.frame(original_index = integer(0), genus = character(0), genus_probability = numeric(0))
  }
  
  # Add genus predictions to summary_df (NA for non-viral contigs)
  summary_df <- summary_df %>%
    left_join(genus_df, by = "original_index") %>%
    select(contig, mean_prediction, median_prediction, sd_prediction, iqr_prediction, classification, num_predictions, contig_length, genus, genus_probability)
  
  # Write detailed results to CSV
  write.table(summary_df, csv_file, sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Compute unique genera with probability > 0.5
  unique_genera <- genus_df %>% 
    filter(genus_probability > 0.5) %>% 
    pull(genus) %>% 
    unique() %>% 
    paste(collapse = ", ")
  
  # Delete temporary binary H5 files
  unlink(temp_output_dir, recursive = TRUE)
  
  message("Processed:", fasta_file)
  
  # Return summary stats for the whole file
  return(list(
    file_name = base_name,
    virus_percentage = ifelse(total_regions > 0, virus_regions / total_regions * 100, 0),
    non_virus_percentage = ifelse(total_regions > 0, non_virus_regions / total_regions * 100, 0),
    viral_nt = viral_nt,
    non_viral_nt = non_viral_nt,
    unique_genera = unique_genera
  ))
}

# Parse command-line arguments
option_list <- list(
  make_option(c("-f", "--fasta_file"), type = "character", default = NULL, help = "Input FASTA file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, help = "Output directory", metavar = "character"),
  make_option(c("-b", "--model_binary"), type = "character", default = NULL, help = "Binary model path", metavar = "character"),
  make_option(c("-g", "--model_genus"), type = "character", default = NULL, help = "Genus model path", metavar = "character"),
  make_option(c("-l", "--genus_labels"), type = "character", default = NULL, help = "Genus labels RDS path", metavar = "character"),
  make_option(c("-w", "--window_size"), type = "integer", default = 1000, help = "Window size [default=1000]", metavar = "integer"),
  make_option(c("-s", "--step"), type = "integer", default = 5000, help = "Step size [default=5000]", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$fasta_file) || is.null(opt$output_dir) || is.null(opt$model_binary) || 
    is.null(opt$model_genus) || is.null(opt$genus_labels)) {
  stop("Missing required arguments. Use --help for usage information.")
}

# Create output directory if it doesn't exist
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# Load models and labels
model <- keras::load_model_hdf5(opt$model_binary, custom_objects = custom_objects)
model_genus <- keras::load_model_hdf5(opt$model_genus, compile = FALSE)
genus_labels <- readRDS(opt$genus_labels)

# Process the single FASTA file
summary_stats <- process_fasta(opt$fasta_file, model, model_genus, genus_labels, opt$output_dir, opt$window_size, opt$step)

# Write summary to a separate CSV file
base_name <- tools::file_path_sans_ext(basename(opt$fasta_file))
summary_csv_file <- file.path(opt$output_dir, paste0(base_name, "_summary.csv"))
summary_df <- data.frame(
  file_name = summary_stats$file_name,
  virus_percentage = summary_stats$virus_percentage,
  non_virus_percentage = summary_stats$non_virus_percentage,
  viral_nt = summary_stats$viral_nt,
  non_viral_nt = summary_stats$non_viral_nt,
  unique_genera = summary_stats$unique_genera
)
write.table(summary_df, summary_csv_file, sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)

log_info(paste("Processing complete for", opt$fasta_file))