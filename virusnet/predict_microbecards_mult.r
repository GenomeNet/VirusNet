#!/usr/bin/env Rscript

# Suppress TensorFlow messages
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)

# Load required libraries in the master process
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
  library(parallel)
}))

source("utils.r")

# Define the process_fasta function for individual FASTA files
process_fasta <- function(fasta_file) {
  # Access variables from the worker's global environment
  model <- .GlobalEnv[["model"]]
  model_genus <- .GlobalEnv[["model_genus"]]
  genus_labels <- .GlobalEnv[["genus_labels"]]
  output_csv <- .GlobalEnv[["output_csv"]]
  window_size <- .GlobalEnv[["window_size"]]
  step <- .GlobalEnv[["step"]]
  
  # Read all contigs from the FASTA file
  fasta_data_all <- readFasta(fasta_file)
  if (nrow(fasta_data_all) == 0) stop("FASTA file is empty.")
  
  # Add an index column and compute lengths
  fasta_data_all$index <- seq_len(nrow(fasta_data_all))
  fasta_data_all$Length <- nchar(fasta_data_all$Sequence)
  fasta_data_all$SanitizedHeader <- make.names(fasta_data_all$Header, unique = TRUE)
  
  # Filter contigs with Length >= window_size
  fasta_data <- fasta_data_all %>% filter(Length >= window_size)
  
  # Define output CSV file path
  base_name <- tools::file_path_sans_ext(basename(fasta_file))
  csv_file <- file.path(output_csv, paste0(base_name, ".csv"))
  
  # Create a unique temporary directory for this worker
  temp_output_dir <- file.path(tempdir(), paste0(base_name, "_", Sys.getpid()))
  dir.create(temp_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Perform binary prediction
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
  
  # Initialize counters
  virus_regions <- 0
  non_virus_regions <- 0
  total_regions <- 0
  viral_nt <- 0
  non_viral_nt <- 0
  
  # Process filtered contigs
  for (i in 1:nrow(fasta_data)) {
    idx <- fasta_data$index[i]
    contig <- fasta_data$Header[i]
    contig_length <- fasta_data$Length[i]
    
    h5_file_path <- file.path(temp_output_dir, paste0("states_nr_", idx, ".h5"))
    if (!file.exists(h5_file_path)) {
      warning(paste("H5 file not found for contig", contig, "at index", idx))
      next
    }
    
    preds_contig <- load_prediction(h5_path = h5_file_path, get_sample_position = TRUE, verbose = FALSE)
    states <- as.data.frame(preds_contig$states)
    colnames(states) <- c("non_virus", "virus")
    virus_probs <- states$virus
    num_preds <- length(preds_contig$sample_end_position)
    
    mean_prob <- mean(virus_probs, na.rm = TRUE)
    median_prob <- median(virus_probs, na.rm = TRUE)
    sd_prob <- sd(virus_probs, na.rm = TRUE)
    iqr_prob <- IQR(virus_probs, na.rm = TRUE)
    classification <- ifelse(mean_prob > 0.5, "virus", "non-virus")
    
    if (mean_prob > 0.5) {
      virus_regions <- virus_regions + 1
      viral_nt <- viral_nt + contig_length
    } else {
      non_virus_regions <- non_virus_regions + 1
      non_viral_nt <- non_viral_nt + contig_length
    }
    total_regions <- total_regions + 1
    
    summary_list[[i]] <- data.frame(
      original_index = idx,
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
  
  # Combine summary results
  summary_df <- bind_rows(summary_list)
  
  # Write viral contigs to a new FASTA file
  viral_indices <- summary_df %>% filter(classification == "virus") %>% pull(original_index)
  viral_contigs <- fasta_data_all %>% filter(index %in% viral_indices)
  viral_fasta_file <- file.path(output_csv, paste0(base_name, "_viral.fasta"))
  
  if (nrow(viral_contigs) > 0) {
    writeFasta(viral_contigs[, c("Header", "Sequence")], viral_fasta_file)
    
    temp_genus_output_dir <- file.path(tempdir(), paste0(base_name, "_genus_", Sys.getpid()))
    dir.create(temp_genus_output_dir, showWarnings = FALSE, recursive = TRUE)
    
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
    
    genus_list <- list()
    for (j in 1:nrow(viral_contigs)) {
      h5_file_path <- file.path(temp_genus_output_dir, paste0("states_nr_", j, ".h5"))
      if (!file.exists(h5_file_path)) {
        warning(paste("Genus H5 file not found for viral contig", j))
        next
      }
      
      preds_genus <- load_prediction(h5_path = h5_file_path, get_sample_position = FALSE, verbose = FALSE)
      states <- preds_genus$states
      if (is.null(states)) {
        warning(paste("No states found in H5 file for viral contig", j))
        next
      }
      
      avg_probs <- colMeans(states)
      max_idx <- which.max(avg_probs)
      predicted_genus <- genus_labels[max_idx]
      genus_probability <- avg_probs[max_idx]
      
      genus_list[[j]] <- data.frame(
        original_index = viral_contigs$index[j],
        genus = predicted_genus,
        genus_probability = genus_probability
      )
    }
    
    genus_df <- bind_rows(genus_list)
    unlink(temp_genus_output_dir, recursive = TRUE)
  } else {
    genus_df <- data.frame(original_index = integer(0), genus = character(0), genus_probability = numeric(0))
  }
  
  # Merge genus predictions with summary
  summary_df <- summary_df %>%
    left_join(genus_df, by = "original_index") %>%
    select(contig, mean_prediction, median_prediction, sd_prediction, iqr_prediction, classification, num_predictions, contig_length, genus, genus_probability)
  
  # Write results to CSV
  write.table(summary_df, csv_file, sep = ";", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Compute unique genera
  unique_genera <- genus_df %>% 
    filter(genus_probability > 0.5) %>% 
    pull(genus) %>% 
    unique() %>% 
    paste(collapse = ", ")
  
  # Clean up temporary files
  unlink(temp_output_dir, recursive = TRUE)
  
  message("Processed:", fasta_file)
  
  # Return summary stats
  return(list(
    file_name = base_name,
    virus_percentage = ifelse(total_regions > 0, virus_regions / total_regions * 100, 0),
    non_virus_percentage = ifelse(total_regions > 0, non_virus_regions / total_regions * 100, 0),
    viral_nt = viral_nt,
    non_viral_nt = non_viral_nt,
    unique_genera = unique_genera
  ))
}

# Main function for parallel processing
main <- function() {
  # Hardcoded parameters (modify as needed)
  input_folder <- "../../../merged"
  model_binary <- "~/.virusnet/transfer_learning_virus_bert_5_85.h5"
  model_genus_path <- "~/.virusnet/virus_genus_2023-01-23.hdf5"
  genus_labels_path <- "~/.virusnet/genus_labels.rds"
  output_csv <- "output_merged_mult"
  window_size <- 1000
  step <- 5000
  num_gpus <- 2  # Adjust based on available GPUs
  
  # Create output directory
  dir.create(output_csv, showWarnings = FALSE, recursive = TRUE)
  
  # List FASTA files
  fasta_files <- list.files(input_folder, pattern = "\\.fasta$", full.names = TRUE)
  if (length(fasta_files) == 0) stop("No FASTA files found in the input folder.")
  
  # Set up parallel cluster
  cl <- makeCluster(num_gpus)
  on.exit(stopCluster(cl))
  
  # Assign GPUs to workers
  clusterApply(cl, 0:(num_gpus - 1), function(gpu_id) {
    Sys.setenv(CUDA_VISIBLE_DEVICES = gpu_id)
  })
  
  # Export variables to workers
  clusterExport(cl, c("output_csv", "window_size", "step", "model_binary", 
                     "model_genus_path", "genus_labels_path"), envir = environment())
  
  # Initialize workers with libraries and models
  clusterEvalQ(cl, {
    Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)
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
    }))
    source("utils.r")
    
    custom_objects <- list(
      "layer_pos_embedding" = layer_pos_embedding,
      "layer_pos_sinusoid" = layer_pos_sinusoid,
      "layer_transformer_block" = layer_transformer_block,
      "layer_aggregate_time_dist" = layer_aggregate_time_dist
    )
    
    model <- keras::load_model_hdf5(model_binary, custom_objects = custom_objects)
    model_genus <- keras::load_model_hdf5(model_genus_path, compile = FALSE)
    genus_labels <- readRDS(genus_labels_path)
  })
  
  # Process files in parallel
  summary_output_list <- parLapply(cl, fasta_files, process_fasta)
  
  # Combine and write summary
  summary_output_df <- bind_rows(summary_output_list)
  summary_output_path <- file.path(output_csv, "summary_output.csv")
  write.table(summary_output_df, summary_output_path, sep = ";", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  cat("Processing complete.\n")
}

# Run main function if not interactive
if (!interactive()) {
  main()
}