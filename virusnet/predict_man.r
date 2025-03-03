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
}))

custom_objects <- list(
  "layer_pos_embedding" = layer_pos_embedding,
  "layer_pos_sinusoid" = layer_pos_sinusoid,
  "layer_transformer_block" = layer_transformer_block,
  "layer_aggregate_time_dist" = layer_aggregate_time_dist
)

# Define utility functions
safe_write_fasta <- function(fasta_subset, file_path) {
  if (!is.null(fasta_subset) && nrow(fasta_subset) > 0 && all(nchar(fasta_subset$Sequence) > 0, na.rm = TRUE)) {
    writeFasta(fasta_subset, file_path)
    return(TRUE)
  }
  return(FALSE)
}

# Function to merge overlapping intervals
merge_intervals <- function(starts, ends) {
  if (length(starts) == 0) return(data.frame(start = numeric(0), end = numeric(0)))
  intervals <- data.frame(start = starts, end = ends) %>% arrange(start)
  merged <- intervals[1, ]
  for (i in 2:nrow(intervals)) {
    if (intervals$start[i] <= merged$end[nrow(merged)] + 1) {  # Merge if overlapping or adjacent
      merged$end[nrow(merged)] <- max(merged$end[nrow(merged)], intervals$end[i])
    } else {
      merged <- rbind(merged, intervals[i, ])
    }
  }
  merged
}

# Logging function
log_info <- function(message) {
  cat("[INFO]", message, "\n")
}

# Command-line arguments (hardcoded for demonstration)
args <- list(
  input = "../../../merged_1contig/GCF_034479145.1_ASM3447914v1_genomic.fasta",  # Replace with your input FASTA path
  model_binary = "~/.virusnet/transfer_learning_virus_bert_5_85.h5",              # Replace with your model path
  output = "prophage_detection_output",
  batch_size = 4,
  verbose = TRUE
)

# Parameters
window_size <- 1000           # Model window size (base pairs)
coarse_step <- 10000          # Coarse step size for initial scan (bp)
fine_step <- 500              # Fine step size for detailed scan (bp)
boundary_step <- 50           # Ultra-fine step size for boundary refinement (bp)
prob_threshold <- 0.7         # Probability threshold for identifying potential prophage regions
boundary_threshold <- 0.5     # Threshold to define prophage boundaries
initial_window <- 5000        # Initial window size around high-probability positions (bp)
min_prophage_size <- 7000     # Minimum prophage size (bp)

# Create output directory
if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

# Load FASTA data
fasta_data <- readFasta(args$input)
if (nrow(fasta_data) == 0) stop("Input FASTA file is empty.")
fasta_data$Length <- nchar(fasta_data$Sequence)
fasta_data <- fasta_data %>% filter(Length >= window_size)
fasta_data$SanitizedHeader <- make.names(fasta_data$Header, unique = TRUE)

# Load the model
model <- keras::load_model_hdf5(args$model_binary, custom_objects = custom_objects)

# --- Step 1: Initial Coarse Sampling ---
log_info("Starting coarse sampling with 10k step size.")
temp_file_coarse <- tempfile(fileext = ".h5")
pred_coarse <- predict_model(
  output_format = "by_entry_one_file",
  model = model,
  layer_name = "dense_26",
  path_input = args$input,
  step = coarse_step,
  batch_size = as.numeric(args$batch_size),
  verbose = FALSE,
  padding = "standard",
  mode = "label",
  format = "fasta",
  filename = temp_file_coarse,
  return_int = length(model$input_shape) == 2
)
coarse_preds <- load_prediction(h5_path = temp_file_coarse, get_sample_position = TRUE, verbose = FALSE)

# Process coarse predictions
coarse_df <- lapply(1:nrow(fasta_data), function(i) {
  states <- as.data.frame(coarse_preds[[i]]$states)
  positions <- coarse_preds[[i]]$sample_end_position
  colnames(states) <- c("non_virus", "virus")
  data.frame(
    contig = fasta_data$SanitizedHeader[i],
    position = positions,
    virus_prob = states$virus,
    type = "coarse"
  )
}) %>% bind_rows()

# Filter positions where virus_prob > 0.7
high_prob_positions <- coarse_df %>% filter(virus_prob > prob_threshold)
log_info(paste("Found", nrow(high_prob_positions), "high-probability positions after coarse sampling."))

# --- Step 2: Merge Initial Windows and Perform Fine-Grained Prediction ---
log_info("Merging initial windows and performing fine-grained prediction.")
all_fine_preds <- list()
for (current_contig in unique(high_prob_positions$contig)) {
  contig_high_prob <- high_prob_positions %>% filter(contig == current_contig)
  contig_length <- fasta_data$Length[fasta_data$SanitizedHeader == current_contig]
  
  # Create initial windows
  starts <- pmax(1, contig_high_prob$position - initial_window)
  ends <- pmin(contig_length, contig_high_prob$position + initial_window)
  
  # Merge overlapping windows
  merged_windows <- merge_intervals(starts, ends)
  
  # Perform fine-grained prediction on each merged window
  for (j in 1:nrow(merged_windows)) {
    window_start <- merged_windows$start[j]
    window_end <- merged_windows$end[j]
    
    # Extract subsequence
    fasta_subset <- fasta_data %>%
      filter(SanitizedHeader == current_contig) %>%
      mutate(
        Sequence = substr(Sequence, window_start, window_end),
        Header = paste0(SanitizedHeader, "_fine_", window_start),
        start_offset = window_start
      )
    
    temp_fasta <- tempfile(fileext = ".fasta")
    safe_write_fasta(fasta_subset %>% select(Header, Sequence), temp_fasta)
    
    # Perform fine-grained prediction
    log_info(paste("Running fine-grained prediction for contig:", current_contig, "window:", window_start, "-", window_end))
    temp_file_fine <- tempfile(fileext = ".h5")
    pred_fine <- predict_model(
      output_format = "by_entry_one_file",
      model = model,
      layer_name = "dense_26",
      path_input = temp_fasta,
      step = fine_step,
      batch_size = as.numeric(args$batch_size),
      verbose = FALSE,
      padding = "standard",
      mode = "label",
      format = "fasta",
      filename = temp_file_fine,
      return_int = length(model$input_shape) == 2
    )
    fine_preds <- load_prediction(h5_path = temp_file_fine, get_sample_position = TRUE, verbose = FALSE)
    
    fine_df <- data.frame(
      contig = current_contig,
      position = fine_preds[[1]]$sample_end_position + fasta_subset$start_offset - 1,
      virus_prob = fine_preds[[1]]$states[, 2],
      type = "fine"
    )
    
    all_fine_preds[[length(all_fine_preds) + 1]] <- fine_df
  }
}
fine_df_combined <- bind_rows(all_fine_preds)

# --- Step 3: Identify High-Probability Regions ---
log_info("Identifying high-probability regions from fine-grained predictions.")
high_prob_regions <- list()
for (current_contig in unique(fine_df_combined$contig)) {
  contig_fine <- fine_df_combined %>% filter(contig == current_contig) %>% arrange(position)
  is_high <- contig_fine$virus_prob >= boundary_threshold
  runs <- rle(is_high)
  run_starts <- cumsum(c(1, runs$lengths[-length(runs$lengths)]))
  run_ends <- cumsum(runs$lengths)
  high_runs <- which(runs$values)
  for (run_idx in high_runs) {
    start_idx <- run_starts[run_idx]
    end_idx <- run_ends[run_idx]
    start_pos <- contig_fine$position[start_idx]
    end_pos <- contig_fine$position[end_idx]
    high_prob_regions[[length(high_prob_regions) + 1]] <- data.frame(
      contig = current_contig,
      approximate_start = start_pos,
      approximate_end = end_pos
    )
  }
}
high_prob_regions_df <- bind_rows(high_prob_regions)

# --- Step 4: Refine Boundaries with Ultra-Fine Sampling ---
log_info("Refining prophage boundaries with ultra-fine sampling.")
all_boundary_preds <- list()  # Initialize list to collect boundary predictions
refined_regions <- list()
for (i in 1:nrow(high_prob_regions_df)) {
  current_contig <- high_prob_regions_df$contig[i]
  approximate_start <- high_prob_regions_df$approximate_start[i]
  approximate_end <- high_prob_regions_df$approximate_end[i]
  contig_length <- fasta_data$Length[fasta_data$SanitizedHeader == current_contig]
  
  # Refine left boundary
  left_refine_start <- max(1, approximate_start - fine_step)
  left_refine_end <- approximate_start
  extract_start <- max(1, left_refine_start - window_size + 1)
  extract_end <- min(contig_length, left_refine_end + window_size - 1)
  fasta_subset <- fasta_data %>%
    filter(SanitizedHeader == current_contig) %>%
    mutate(
      Sequence = substr(Sequence, extract_start, extract_end),
      Header = paste0(SanitizedHeader, "_refine_left_", extract_start),
      start_offset = extract_start
    )
  temp_fasta <- tempfile(fileext = ".fasta")
  safe_write_fasta(fasta_subset %>% select(Header, Sequence), temp_fasta)
  
  temp_file_refine <- tempfile(fileext = ".h5")
  pred_refine <- predict_model(
    output_format = "by_entry_one_file",
    model = model,
    layer_name = "dense_26",
    path_input = temp_fasta,
    step = boundary_step,
    batch_size = as.numeric(args$batch_size),
    verbose = FALSE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = temp_file_refine,
    return_int = length(model$input_shape) == 2
  )
  refine_preds <- load_prediction(h5_path = temp_file_refine, get_sample_position = TRUE, verbose = FALSE)
  
  refine_df <- data.frame(
    position = refine_preds[[1]]$sample_end_position + fasta_subset$start_offset - 1,
    virus_prob = refine_preds[[1]]$states[, 2]
  )
  # Save boundary predictions
  boundary_left_df <- data.frame(
    contig = current_contig,
    position = refine_preds[[1]]$sample_end_position + fasta_subset$start_offset - 1,
    virus_prob = refine_preds[[1]]$states[, 2],
    type = "boundary"
  )
  all_boundary_preds[[length(all_boundary_preds) + 1]] <- boundary_left_df
  
  refined_start_candidates <- refine_df %>% filter(virus_prob >= boundary_threshold)
  refined_start <- if (nrow(refined_start_candidates) > 0) min(refined_start_candidates$position) else approximate_start
  
  # Refine right boundary
  right_refine_start <- approximate_end
  right_refine_end <- min(contig_length, approximate_end + fine_step)
  extract_start <- max(1, right_refine_start - window_size + 1)
  extract_end <- min(contig_length, right_refine_end + window_size - 1)
  fasta_subset <- fasta_data %>%
    filter(SanitizedHeader == current_contig) %>%
    mutate(
      Sequence = substr(Sequence, extract_start, extract_end),
      Header = paste0(SanitizedHeader, "_refine_right_", extract_start),
      start_offset = extract_start
    )
  temp_fasta <- tempfile(fileext = ".fasta")
  safe_write_fasta(fasta_subset %>% select(Header, Sequence), temp_fasta)
  
  temp_file_refine <- tempfile(fileext = ".h5")
  pred_refine <- predict_model(
    output_format = "by_entry_one_file",
    model = model,
    layer_name = "dense_26",
    path_input = temp_fasta,
    step = boundary_step,
    batch_size = as.numeric(args$batch_size),
    verbose = FALSE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = temp_file_refine,
    return_int = length(model$input_shape) == 2
  )
  refine_preds <- load_prediction(h5_path = temp_file_refine, get_sample_position = TRUE, verbose = FALSE)
  
  refine_df <- data.frame(
    position = refine_preds[[1]]$sample_end_position + fasta_subset$start_offset - 1,
    virus_prob = refine_preds[[1]]$states[, 2]
  )
  # Save boundary predictions
  boundary_right_df <- data.frame(
    contig = current_contig,
    position = refine_preds[[1]]$sample_end_position + fasta_subset$start_offset - 1,
    virus_prob = refine_preds[[1]]$states[, 2],
    type = "boundary"
  )
  all_boundary_preds[[length(all_boundary_preds) + 1]] <- boundary_right_df
  
  refined_end_candidates <- refine_df %>% filter(virus_prob >= boundary_threshold)
  refined_end <- if (nrow(refined_end_candidates) > 0) max(refined_end_candidates$position) else approximate_end
  
  # Store refined region if it meets the minimum size
  if (!is.na(refined_start) && !is.na(refined_end) && refined_end - refined_start >= min_prophage_size) {
    refined_regions[[i]] <- data.frame(
      contig = current_contig,
      start = refined_start,
      end = refined_end
    )
  }
}
refined_regions_df <- bind_rows(refined_regions)

# Combine all predictions (coarse, fine, and boundary) for plotting
boundary_df_combined <- bind_rows(all_boundary_preds)
all_preds <- bind_rows(coarse_df, fine_df_combined, boundary_df_combined)

# --- Step 5: Output Results ---
log_info("Completed analysis, writing results.")
write.csv(coarse_df, file.path(args$output, "coarse_predictions.csv"), row.names = FALSE)
write.csv(fine_df_combined, file.path(args$output, "fine_predictions.csv"), row.names = FALSE)
write.csv(refined_regions_df, file.path(args$output, "prophage_regions.csv"), row.names = FALSE)

# Export prophage sequences to FASTA
for (i in 1:nrow(refined_regions_df)) {
  current_contig <- refined_regions_df$contig[i]
  start <- refined_regions_df$start[i]
  end <- refined_regions_df$end[i]
  fasta_subset <- fasta_data %>%
    filter(SanitizedHeader == current_contig) %>%
    mutate(
      Sequence = substr(Sequence, start, end),
      Header = paste0(SanitizedHeader, "_prophage_", start, "_", end)
    )
  safe_write_fasta(
    fasta_subset %>% select(Header, Sequence),
    file.path(args$output, paste0("prophage_", i, ".fasta"))
  )
}

# --- Step 6: Plotting Full Contigs ---
log_info("Generating plots for each contig.")
for (current_contig in unique(all_preds$contig)) {
  # Subset predictions for the current contig
  contig_preds <- all_preds %>% filter(contig == current_contig)
  
  # Subset refined regions for the current contig
  contig_regions <- refined_regions_df %>% filter(contig == current_contig)
  
  # Get contig length
  contig_length <- fasta_data$Length[fasta_data$SanitizedHeader == current_contig]
  
  # Create background_df
  if (nrow(contig_regions) > 0) {
    contig_regions <- contig_regions %>% arrange(start)
    non_virus_intervals <- list()
    if (contig_regions$start[1] > 1) {
      non_virus_intervals[[1]] <- data.frame(
        start = 1,
        end = contig_regions$start[1] - 1,
        classification = "non-virus"
      )
    }
    for (k in 1:(nrow(contig_regions) - 1)) {
      if (contig_regions$end[k] < contig_regions$start[k + 1] - 1) {
        non_virus_intervals[[length(non_virus_intervals) + 1]] <- data.frame(
          start = contig_regions$end[k] + 1,
          end = contig_regions$start[k + 1] - 1,
          classification = "non-virus"
        )
      }
    }
    if (contig_regions$end[nrow(contig_regions)] < contig_length) {
      non_virus_intervals[[length(non_virus_intervals) + 1]] <- data.frame(
        start = contig_regions$end[nrow(contig_regions)] + 1,
        end = contig_length,
        classification = "non-virus"
      )
    }
    non_virus_df <- bind_rows(non_virus_intervals)
    virus_df <- contig_regions %>% mutate(classification = "virus")
    background_df <- bind_rows(non_virus_df, virus_df) %>% arrange(start)
  } else {
    background_df <- data.frame(
      start = 1,
      end = contig_length,
      classification = "non-virus"
    )
  }
  
  # Create the plot
  p <- ggplot() +
    geom_rect(data = background_df, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = classification), alpha = 0.3) +
    scale_fill_manual(values = c("virus" = "red", "non-virus" = "green")) +
    geom_point(data = contig_preds, aes(x = position, y = virus_prob, color = type), alpha = 0.5) +
    scale_color_manual(values = c("coarse" = "blue", "fine" = "purple", "boundary" = "orange")) +  # Updated to include boundary
    geom_line(data = contig_preds %>% arrange(position), aes(x = position, y = virus_prob), color = "black") +
    geom_text(data = contig_regions %>% mutate(prophage_number = row_number()), 
              aes(x = (start + end) / 2, y = 1.05, label = paste("Prophage", prophage_number)), 
              size = 3, color = "black") +
    labs(title = paste("Contig:", current_contig),
         subtitle = "Virus Probability and Prophage Regions",
         x = "Position (nt)",
         y = "Virus Probability",
         fill = "Classification",
         color = "Prediction Type") +
    theme_minimal() +
    ylim(0, 1.1)
  
  # Save the plot
  ggsave(file.path(args$output, paste0(current_contig, "_prophage_plot.pdf")), p, width = 12, height = 6)
}

# --- Step 7: Plotting Individual Prophage Regions ---
log_info("Generating individual PDF plots for each prophage region with ±500 nt.")
for (i in 1:nrow(refined_regions_df)) {
  current_contig <- refined_regions_df$contig[i]
  prophage_start <- refined_regions_df$start[i]
  prophage_end <- refined_regions_df$end[i]
  
  # Define the extended range (±500 nt)
  contig_length <- fasta_data$Length[fasta_data$SanitizedHeader == current_contig]
  plot_start <- max(1, prophage_start - 500)
  plot_end <- min(contig_length, prophage_end + 500)
  
  # Subset predictions for this range
  prophage_preds <- all_preds %>% 
    filter(contig == current_contig, position >= plot_start, position <= plot_end)
  
  # Create background data frame for shading
  background_df <- data.frame(
    start = c(plot_start, prophage_start),
    end = c(prophage_start - 1, prophage_end),
    classification = c("non-virus", "virus")
  )
  if (prophage_end < plot_end) {
    background_df <- bind_rows(
      background_df,
      data.frame(start = prophage_end + 1, end = plot_end, classification = "non-virus")
    )
  }
  
  # Create the plot
  p <- ggplot() +
    geom_rect(data = background_df, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = classification), alpha = 0.3) +
    scale_fill_manual(values = c("virus" = "red", "non-virus" = "green")) +
    geom_point(data = prophage_preds, aes(x = position, y = virus_prob, color = type), alpha = 0.5) +
    scale_color_manual(values = c("coarse" = "blue", "fine" = "purple", "boundary" = "orange")) +  # Updated to include boundary
    geom_line(data = prophage_preds %>% arrange(position), aes(x = position, y = virus_prob), color = "black") +
    labs(title = paste("Prophage", i, "in Contig:", current_contig),
         subtitle = paste("Region:", prophage_start, "-", prophage_end, "with ±500 nt"),
         x = "Position (nt)",
         y = "Virus Probability",
         fill = "Classification",
         color = "Prediction Type") +
    theme_minimal() +
    ylim(0, 1.1)
  
  # Save the plot
  ggsave(file.path(args$output, paste0("prophage_", i, "_plot.pdf")), p, width = 8, height = 4)
}

cat("Processing complete. Check the output directory:", args$output, "\n")