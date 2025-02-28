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
  library(HiddenMarkov)
  library(HMM)
}))

# Source utility scripts
conda_prefix <- Sys.getenv("CONDA_PREFIX")
invisible(source(file.path(conda_prefix, "bin", "utils.r")))
invisible(source(file.path(conda_prefix, "bin", "setup_logger.r")))
source("utils.r")
source("setup_logger.r")

# Command arguments
command_args <- list()
command_args$output <- "output_r3"
command_args$input <- "/mnt/md0/pmuench/virusnet/merged_1contig/GCA_900637335.1_50618_H02_genomic.fasta"
command_args$model_binary <- "~/.virusnet/transfer_learning_virus_bert_5_85.h5"
command_args$batch_size <- 4
command_args$verbose <- "TRUE"

# Set minimum prophage size threshold (in base pairs)
prophage_min_size <- 5000  # Minimum size to consider as a prophage region

# Function to safely write FASTA
safe_write_fasta <- function(fasta_subset, file_path) {
  if (!is.null(fasta_subset) && nrow(fasta_subset) > 0 && all(nchar(fasta_subset$Sequence) > 0, na.rm = TRUE)) {
    writeFasta(fasta_subset, file_path)
    custom_log("INFO", paste("FASTA data written to:", file_path), "3/8")
    return(TRUE)
  } else {
    custom_log("WARN", "No valid contigs to write", "3/8")
    return(FALSE)
  }
}

# Create output directory
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

# Validate and clean FASTA data
fasta_data$Length <- nchar(fasta_data$Sequence)
fasta_data <- fasta_data %>% filter(Length > 0 & !is.na(Sequence) & !is.na(Header))
total_nt <- sum(fasta_data$Length)
custom_log("INFO", paste("Total nucleotides after validation:", total_nt), "1/8")
custom_log("INFO", paste("Contig lengths:", paste(fasta_data$Length, collapse = ", ")), "1/8")

# Create sanitized versions of contig names for data frame operations
fasta_data$SanitizedHeader <- make.names(fasta_data$Header, unique = TRUE)
header_mapping <- setNames(fasta_data$Header, fasta_data$SanitizedHeader)
reverse_header_mapping <- setNames(fasta_data$SanitizedHeader, fasta_data$Header)

# Add diagnostic logging for header sanitization
custom_log("INFO", "=== DIAGNOSTIC: Header Sanitization ===", "1/8")
for (i in 1:nrow(fasta_data)) {
  custom_log("INFO", paste("Original:", fasta_data$Header[i], "-> Sanitized:", fasta_data$SanitizedHeader[i]), "1/8")
}

custom_log("INFO", "Loading binary model", "2/8")
suppressWarnings({
  model_binary <- keras::load_model_hdf5(command_args$model_binary, custom_objects = custom_objects)
})

# --- Step 1: Initial Coarse Prediction with Step Size 10000 ---
custom_log("INFO", paste("Performing initial coarse predictions with step size 10000, total nt:", total_nt), "3/8")
temp_file_coarse <- tempfile(fileext = ".h5")
pred_coarse <- predict_model(
  output_format = "by_entry_one_file",
  model = model_binary,
  layer_name = "dense_26",
  path_input = command_args$input,
  round_digits = 3,
  step = 10000,
  batch_size = as.numeric(command_args$batch_size),
  verbose = FALSE,
  return_states = FALSE,
  padding = "standard",
  mode = "label",
  format = "fasta",
  filename = temp_file_coarse,
  return_int = length(model_binary$input_shape) == 2
)

custom_log("INFO", "Loading coarse predictions", "3/8")
prediction_df_coarse <- load_prediction(h5_path = temp_file_coarse, get_sample_position = TRUE, verbose = FALSE)

custom_log("INFO", "Processing coarse predictions", "3/8")
coarse_results <- lapply(1:num_fasta_entries, function(i) {
  states <- as.data.frame(prediction_df_coarse[[i]]$states)
  positions <- prediction_df_coarse[[i]]$sample_end_position
  colnames(states) <- c("is_non_virus", "is_virus")
  starts <- pmax(1, positions - 10000 + 1)
  data.frame(
    contig_name = fasta_data$SanitizedHeader[i],
    original_name = fasta_data$Header[i],
    start = starts,
    end = positions,
    width = 10000,
    is_virus = states$is_virus
  )
})
coarse_df <- do.call(rbind, coarse_results)

# Initial classification with adjusted thresholds
coarse_summary <- coarse_df %>%
  group_by(contig_name) %>%
  summarise(
    original_name = first(original_name),
    mean_is_virus = mean(is_virus, na.rm = TRUE),
    min_is_virus = min(is_virus, na.rm = TRUE),
    max_is_virus = max(is_virus, na.rm = TRUE),
    sd_is_virus = sd(is_virus, na.rm = TRUE),
    max_diff = if (n() > 1) max(abs(diff(is_virus)), na.rm = TRUE) else 0,
    case = case_when(
      mean_is_virus < 0.1 & max_diff < 0.2 ~ "non-virus",
      mean_is_virus > 0.9 & max_diff < 0.2 ~ "virus",
      TRUE ~ "mixed"
    )
  ) %>%
  ungroup()

custom_log("INFO", "Step 1 prediction stats:", "3/8")
custom_log("INFO", paste("Mean is_virus distribution:", paste(round(coarse_summary$mean_is_virus, 3), collapse = ", ")), "3/8")
custom_log("INFO", paste("Min is_virus distribution:", paste(round(coarse_summary$min_is_virus, 3), collapse = ", ")), "3/8")
custom_log("INFO", paste("Max is_virus distribution:", paste(round(coarse_summary$max_is_virus, 3), collapse = ", ")), "3/8")
custom_log("INFO", paste("SD is_virus distribution:", paste(round(coarse_summary$sd_is_virus, 3), collapse = ", ")), "3/8")
custom_log("INFO", paste("Max diff distribution:", paste(round(coarse_summary$max_diff, 3), collapse = ", ")), "3/8")
custom_log("INFO", paste("Case counts - non-virus:", sum(coarse_summary$case == "non-virus"), "virus:", sum(coarse_summary$case == "virus"), "mixed:", sum(coarse_summary$case == "mixed")), "3/8")

refined_contigs <- coarse_summary$contig_name[coarse_summary$case == "mixed" & !is.na(coarse_summary$contig_name)]
if (any(is.na(coarse_summary$contig_name))) {
  custom_log("WARN", "NA contig names detected in coarse summary", "3/8")
}

# --- Step 2: Targeted Medium Refinement with Step Size 5000 ---
medium_df <- data.frame()
if (length(refined_contigs) > 0) {
  refined_regions <- list()
  refined_nt <- 0
  for (contig_sanitized in refined_contigs) {
    # Get original contig name for display and sequence lookup
    original_contig_name <- coarse_summary$original_name[coarse_summary$contig_name == contig_sanitized]
    
    contig_preds <- coarse_df %>% 
      filter(contig_name == contig_sanitized) %>% 
      arrange(start)
      
    consecutive_viral_regions <- contig_preds %>%
      mutate(
        is_viral = is_virus >= 0.5,
        run_id = cumsum(c(1, diff(is_viral != lag(is_viral, default = FALSE))))
      ) %>%
      filter(is_viral) %>%
      group_by(run_id) %>%
      summarise(start = min(start), end = max(end), .groups = "drop")
    
    if (nrow(consecutive_viral_regions) > 0) {
      # Use original name for display but sanitized name for data operations
      contig_seq <- fasta_data$Sequence[fasta_data$SanitizedHeader == contig_sanitized]
      seq_length <- nchar(contig_seq)
      custom_log("INFO", paste("Contig", original_contig_name, "has", nrow(consecutive_viral_regions), "consecutive viral regions"), "3/8")
      merged_start <- min(consecutive_viral_regions$start)
      merged_end <- max(consecutive_viral_regions$end)
      start_pos <- max(1, merged_start - 5000)
      end_pos <- min(merged_end + 5000, seq_length)
      subseq <- substr(contig_seq, start_pos, end_pos)
      region_nt <- end_pos - start_pos + 1
      refined_regions <- c(refined_regions, list(data.frame(
        Header = paste0(contig_sanitized, "_region"),
        Sequence = subseq,
        Start = start_pos,
        End = end_pos,
        OriginalName = original_contig_name
      )))
      refined_nt <- refined_nt + region_nt
      custom_log("INFO", paste("Region for", original_contig_name, ":", start_pos, "-", end_pos, "nt:", region_nt), "3/8")
    }
  }
  
  if (length(refined_regions) > 0) {
    refined_fasta <- do.call(rbind, refined_regions)
    custom_log("INFO", paste("Performing targeted medium refinement for", nrow(refined_fasta), "regions in", length(refined_contigs), "mixed contigs with step size 5000, total nt:", refined_nt), "3/8")
    
    temp_refine_file_5000 <- tempfile(fileext = ".fasta")
    # Write FASTA with sanitized header to avoid issues
    refined_fasta_to_write <- refined_fasta %>% select(Header, Sequence)
    if (safe_write_fasta(refined_fasta_to_write, temp_refine_file_5000)) {
      custom_log("INFO", "Running predict_model for medium refinement", "3/8")
      temp_file_5000 <- tempfile(fileext = ".h5")
      pred_5000 <- predict_model(
        output_format = "by_entry_one_file",
        model = model_binary,
        layer_name = "dense_26",
        path_input = temp_refine_file_5000,
        round_digits = 3,
        step = 5000,
        batch_size = as.numeric(command_args$batch_size),
        verbose = FALSE,
        return_states = FALSE,
        padding = "standard",
        mode = "label",
        format = "fasta",
        filename = temp_file_5000,
        return_int = length(model_binary$input_shape) == 2
      )
      
      custom_log("INFO", "Loading medium predictions", "3/8")
      prediction_df_5000 <- load_prediction(h5_path = temp_file_5000, get_sample_position = TRUE, verbose = FALSE)
      if (length(prediction_df_5000) != nrow(refined_fasta)) {
        custom_log("ERROR", "Mismatch between number of FASTA entries and predictions in Step 2", "3/8")
        stop("Prediction length does not match FASTA input.")
      }
      
      custom_log("INFO", "Processing medium predictions", "3/8")
      prediction_df_5000_named <- setNames(prediction_df_5000, refined_fasta$Header)
      
      medium_results <- lapply(1:nrow(refined_fasta), function(i) {
        states <- as.data.frame(prediction_df_5000[[i]]$states)
        if (nrow(states) == 0) {
          custom_log("WARN", paste("Empty states for region:", refined_fasta$Header[i]), "3/8")
          return(NULL)
        }
        positions <- prediction_df_5000[[i]]$sample_end_position
        colnames(states) <- c("is_non_virus", "is_virus")
        starts <- pmax(1, positions - 5000 + 1) + refined_fasta$Start[i] - 1
        ends <- refined_fasta$Start[i] - 1 + positions
        
        # Extract the sanitized contig name from the region header
        sanitized_contig_name <- sub("_region$", "", refined_fasta$Header[i])
        
        data.frame(
          contig_name = sanitized_contig_name,
          original_name = refined_fasta$OriginalName[i],
          start = starts,
          end = ends,
          width = 5000,
          is_virus = states$is_virus
        )
      })
      
      medium_results <- Filter(Negate(is.null), medium_results)
      if (length(medium_results) > 0) {
        medium_df <- do.call(rbind, medium_results)
      } else {
        custom_log("WARN", "No valid medium predictions generated", "3/8")
      }
    }
  } else {
    custom_log("INFO", "No viral regions identified for medium refinement", "3/8")
  }
} else {
  custom_log("INFO", "No contigs require medium refinement", "3/8")
}

# --- Categorize Contigs After Step 2 ---
combined_df <- bind_rows(
  coarse_df %>% filter(!(contig_name %in% refined_contigs)),
  medium_df
)

# Add diagnostic logging for combined_df
custom_log("INFO", "=== DIAGNOSTIC: Combined DF ===", "3/8")
custom_log("INFO", paste("Combined DF rows:", nrow(combined_df)), "3/8")
custom_log("INFO", paste("Combined DF unique contig_names:", paste(unique(combined_df$contig_name), collapse=", ")), "3/8")

contig_classification <- combined_df %>%
  group_by(contig_name) %>%
  summarise(
    original_name = first(original_name),
    mean_is_virus = mean(is_virus, na.rm = TRUE),
    sd_is_virus = sd(is_virus, na.rm = TRUE),
    case = case_when(
      mean_is_virus < 0.2 ~ "non-virus",
      mean_is_virus > 0.8 ~ "virus",
      TRUE ~ "prophage"
    )
  ) %>%
  ungroup()

# Add diagnostic logging for contig_classification
custom_log("INFO", "=== DIAGNOSTIC: Contig Classification ===", "3/8")
custom_log("INFO", paste("Classification rows:", nrow(contig_classification)), "3/8")
custom_log("INFO", paste("Classification unique contig_names:", paste(unique(contig_classification$contig_name), collapse=", ")), "3/8")

# Ensure all contigs are classified
custom_log("INFO", "=== DIAGNOSTIC: Checking for missing contigs ===", "3/8")
custom_log("INFO", paste("fasta_data$SanitizedHeader:", paste(fasta_data$SanitizedHeader, collapse=", ")), "3/8")
custom_log("INFO", paste("contig_classification$contig_name:", paste(contig_classification$contig_name, collapse=", ")), "3/8")

# Use tryCatch to catch any errors in the setdiff operation
tryCatch({
  missing_contigs_sanitized <- setdiff(fasta_data$SanitizedHeader, contig_classification$contig_name)
  custom_log("INFO", "=== DIAGNOSTIC: Missing Contigs ===", "3/8")
  custom_log("INFO", paste("Missing sanitized contigs:", paste(missing_contigs_sanitized, collapse=", ")), "3/8")
  custom_log("INFO", paste("All sanitized headers:", paste(fasta_data$SanitizedHeader, collapse=", ")), "3/8")
  custom_log("INFO", paste("All classified contigs:", paste(contig_classification$contig_name, collapse=", ")), "3/8")
  
  if (length(missing_contigs_sanitized) > 0) {
    # Create a mapping of sanitized to original headers for the missing contigs
    custom_log("INFO", "Creating mapping for missing contigs", "3/8")
    for (i in seq_along(missing_contigs_sanitized)) {
      sanitized <- missing_contigs_sanitized[i]
      original <- fasta_data$Header[fasta_data$SanitizedHeader == sanitized]
      custom_log("INFO", paste("Mapping", i, ":", sanitized, "->", original), "3/8")
    }
    
    missing_original_headers <- fasta_data$Header[match(missing_contigs_sanitized, fasta_data$SanitizedHeader)]
    custom_log("INFO", paste("Missing original headers:", paste(missing_original_headers, collapse=", ")), "3/8")
    
    missing_df <- data.frame(
      contig_name = missing_contigs_sanitized,
      original_name = missing_original_headers,
      mean_is_virus = 0,
      sd_is_virus = 0,
      case = "non-virus"
    )
    custom_log("INFO", "Created missing_df", "3/8")
    custom_log("INFO", paste("missing_df rows:", nrow(missing_df)), "3/8")
    
    contig_classification <- bind_rows(contig_classification, missing_df)
    custom_log("INFO", "Added missing contigs to contig_classification", "3/8")
    custom_log("INFO", paste("contig_classification rows after adding missing:", nrow(contig_classification)), "3/8")
    custom_log("WARN", paste("Added", length(missing_contigs_sanitized), "missing contigs as non-virus"), "3/8")
  } else {
    custom_log("INFO", "No missing contigs found", "3/8")
  }
}, error = function(e) {
  custom_log("ERROR", paste("Error in missing contigs processing:", e$message), "3/8")
  print(traceback())
})

# --- Step 3: Fine-Grained Boundary Refinement for Prophage Case (Batched) ---
fine_df <- data.frame(
  contig_name = character(),
  original_name = character(),
  start = integer(),
  end = integer(),
  width = integer(),
  is_virus = numeric(),
  stringsAsFactors = FALSE
)

prophage_contigs <- contig_classification$contig_name[contig_classification$case == "prophage"]
if (length(prophage_contigs) > 0) {
  custom_log("INFO", "=== DIAGNOSTIC: Prophage Contigs ===", "3/8")
  custom_log("INFO", paste("Number of prophage contigs:", length(prophage_contigs)), "3/8")
  custom_log("INFO", paste("Prophage contigs:", paste(prophage_contigs, collapse=", ")), "3/8")
  
  custom_log("INFO", paste("Processing", length(prophage_contigs), "prophage contigs"), "3/8")
  
  # Add tryCatch to identify where the error occurs
  tryCatch({
    custom_log("INFO", "Filtering prophage_df", "3/8")
    prophage_df <- combined_df %>% filter(contig_name %in% prophage_contigs)
    custom_log("INFO", paste("prophage_df rows:", nrow(prophage_df)), "3/8")
    
    if (nrow(prophage_df) > 0) {
      custom_log("INFO", "Creating virus_regions", "3/8")
      virus_regions <- prophage_df %>%
        group_by(contig_name) %>%
        mutate(
          is_virus_binary = is_virus >= 0.5,
          transition = is_virus_binary != lag(is_virus_binary, default = FALSE) | 
                      is_virus_binary != lead(is_virus_binary, default = FALSE)
        ) %>%
        filter(transition | is_virus_binary) %>%
        summarise(
          boundaries = list(unique(c(start[transition], end[transition], start[is_virus_binary], end[is_virus_binary]))),
          original_name = first(original_name),
          .groups = "drop"
        )
      custom_log("INFO", paste("virus_regions rows:", nrow(virus_regions)), "3/8")
      
      fine_nt <- 0
      for (i in 1:nrow(virus_regions)) {
        tryCatch({
          sanitized_contig_name <- virus_regions$contig_name[i]
          original_contig_name <- virus_regions$original_name[i]
          custom_log("INFO", paste("Processing virus region for contig:", sanitized_contig_name, "->", original_contig_name), "3/8")
          
          boundaries <- sort(unique(virus_regions$boundaries[[i]]))
          custom_log("INFO", paste("Number of boundaries:", length(boundaries)), "3/8")
          
          # Use sanitized header for lookup in fasta_data
          custom_log("INFO", paste("Looking up sequence for sanitized header:", sanitized_contig_name), "3/8")
          contig_seq <- fasta_data$Sequence[fasta_data$SanitizedHeader == sanitized_contig_name]
          seq_length <- nchar(contig_seq)
          custom_log("INFO", paste("Sequence length:", seq_length), "3/8")
          
          custom_log("INFO", paste("Refining", length(boundaries), "boundaries for contig", original_contig_name), "3/8")
          filtered_boundaries <- c()
          last_boundary <- -Inf
          for (b in boundaries) {
            if (b - last_boundary > 500) {
              filtered_boundaries <- c(filtered_boundaries, b)
              last_boundary <- b
            }
          }
          
          boundary_fasta <- data.frame(
            Header = character(), 
            Sequence = character(), 
            Start = integer(), 
            OriginalName = character(), 
            stringsAsFactors = FALSE
          )
          
          for (boundary in filtered_boundaries) {
            start_pos <- max(1, boundary - 200)
            end_pos <- min(boundary + 200, seq_length)
            region_length <- end_pos - start_pos + 1
            if (region_length < 100) {
              custom_log("WARN", paste("Boundary region too short for contig", original_contig_name, "at", boundary, ":", region_length, "nt"), "3/8")
              next
            }
            subseq <- substr(contig_seq, start_pos, end_pos)
            boundary_fasta <- rbind(boundary_fasta, data.frame(
              Header = paste0(sanitized_contig_name, "_boundary_", boundary),
              Sequence = subseq,
              Start = start_pos,
              OriginalName = original_contig_name
            ))
            fine_nt <- fine_nt + region_length
          }
          
          if (nrow(boundary_fasta) > 0) {
            temp_fasta_file <- tempfile(fileext = ".fasta")
            custom_log("INFO", paste("Writing batched boundary FASTA for", original_contig_name, "with", nrow(boundary_fasta), "sequences, total nt:", fine_nt), "3/8")
            boundary_fasta_to_write <- boundary_fasta %>% select(Header, Sequence)
            if (!safe_write_fasta(boundary_fasta_to_write, temp_fasta_file)) next
            
            temp_file_fine <- tempfile(fileext = ".h5")
            custom_log("INFO", paste("Running batched predict_model for", original_contig_name), "3/8")
            pred_fine <- predict_model(
              output_format = "by_entry_one_file",
              model = model_binary,
              layer_name = "dense_26",
              path_input = temp_fasta_file,
              round_digits = 3,
              step = 100,
              batch_size = as.numeric(command_args$batch_size),
              verbose = FALSE,
              return_states = FALSE,
              padding = "standard",
              mode = "label",
              format = "fasta",
              filename = temp_file_fine,
              return_int = length(model_binary$input_shape) == 2
            )
            
            custom_log("INFO", paste("Loading batched fine predictions for", original_contig_name), "3/8")
            prediction_df_fine <- load_prediction(h5_path = temp_file_fine, get_sample_position = TRUE, verbose = FALSE)
            if (length(prediction_df_fine) != nrow(boundary_fasta)) {
              custom_log("ERROR", paste("Mismatch in batched predictions for", original_contig_name, ":", length(prediction_df_fine), "vs", nrow(boundary_fasta)), "3/8")
              next
            }
            
            prediction_df_fine_named <- setNames(prediction_df_fine, boundary_fasta$Header)
            for (j in 1:nrow(boundary_fasta)) {
              states <- as.data.frame(prediction_df_fine[[j]]$states)
              if (nrow(states) == 0) {
                custom_log("WARN", paste("Empty states for boundary", boundary_fasta$Header[j]), "3/8")
                next
              }
              positions <- prediction_df_fine[[j]]$sample_end_position
              colnames(states) <- c("is_non_virus", "is_virus")
              fine_starts <- boundary_fasta$Start[j] + positions - 100 + 1
              fine_ends <- boundary_fasta$Start[j] + positions
              fine_df <- rbind(fine_df, data.frame(
                contig_name = sanitized_contig_name,  # Ensure we use sanitized name
                original_name = boundary_fasta$OriginalName[j],
                start = fine_starts,
                end = fine_ends,
                width = 100,
                is_virus = states$is_virus
              ))
            }
          }
        }, error = function(e) {
          custom_log("ERROR", paste("Error processing boundary for contig", i, ":", e$message), "3/8")
        })
      }
      custom_log("INFO", paste("Refined boundaries for prophage contigs with step size 100, total nt:", fine_nt), "3/8")
    } else {
      custom_log("INFO", "No rows in prophage_df, skipping refinement", "3/8")
    }
  }, error = function(e) {
    custom_log("ERROR", paste("Error in prophage processing:", e$message), "3/8")
    print(traceback())
  })
}

# --- Generate Final Output with Merged Regions and Probabilities ---
# First, deduplicate and combine overlapping predictions
tryCatch({
  custom_log("INFO", "Binding all predictions", "3/8")
  
  # Check if fine_df is empty and handle accordingly
  if (nrow(fine_df) == 0) {
    custom_log("INFO", "No fine-grained predictions to include", "3/8")
    all_preds <- bind_rows(
      coarse_df %>% select(contig_name, original_name, start, end, is_virus),
      medium_df %>% select(contig_name, original_name, start, end, is_virus)
    )
  } else {
    all_preds <- bind_rows(
      coarse_df %>% select(contig_name, original_name, start, end, is_virus),
      medium_df %>% select(contig_name, original_name, start, end, is_virus),
      fine_df %>% select(contig_name, original_name, start, end, is_virus)
    )
  }
  
  # Arrange and remove duplicates only if we have data
  if (nrow(all_preds) > 0) {
    all_preds <- all_preds %>%
      arrange(contig_name, start, end) %>%
      distinct()
  }
  
  # Add diagnostic logging for all_preds
  custom_log("INFO", "=== DIAGNOSTIC: All Predictions ===", "3/8")
  custom_log("INFO", paste("All predictions rows:", nrow(all_preds)), "3/8")
  custom_log("INFO", paste("All predictions unique contig_names:", paste(unique(all_preds$contig_name), collapse=", ")), "3/8")
}, error = function(e) {
  custom_log("ERROR", paste("Error in binding all predictions:", e$message), "3/8")
  # Create a fallback empty dataframe with the correct structure
  all_preds <<- data.frame(
    contig_name = character(),
    original_name = character(),
    start = integer(),
    end = integer(),
    is_virus = numeric(),
    stringsAsFactors = FALSE
  )
  custom_log("WARN", "Created empty fallback prediction dataframe", "3/8")
  print(traceback())
})

final_output <- data.frame(
  contig_name = character(),
  original_name = character(),
  start = integer(),
  end = integer(),
  length = integer(),
  prediction = character(),
  probability_virus = numeric(),
  mean_prob = numeric(),
  sd_prob = numeric(),
  n_merged = integer(),
  stringsAsFactors = FALSE
)

custom_log("INFO", "Generating final output", "3/8")
for (sanitized_header in unique(fasta_data$SanitizedHeader)) {
  tryCatch({
    original_header <- fasta_data$Header[fasta_data$SanitizedHeader == sanitized_header]
    custom_log("INFO", paste("Processing contig:", sanitized_header, "->", original_header), "3/8")
    
    contig_case <- contig_classification$case[contig_classification$contig_name == sanitized_header]
    contig_mean <- contig_classification$mean_is_virus[contig_classification$contig_name == sanitized_header]
    contig_len <- fasta_data$Length[fasta_data$SanitizedHeader == sanitized_header]
    
    custom_log("INFO", paste("  Case:", contig_case, "Mean:", contig_mean, "Length:", contig_len), "3/8")
    
    if (is.na(contig_case)) {
      contig_case <- "non-virus"
      contig_mean <- 0
      custom_log("WARN", paste("Contig", original_header, "missing classification, defaulting to non-virus"), "3/8")
    }
    
    if (contig_case %in% c("non-virus", "virus")) {
      final_output <- rbind(final_output, data.frame(
        contig_name = sanitized_header,
        original_name = original_header,
        start = 1,
        end = contig_len,
        length = contig_len,
        prediction = contig_case,
        probability_virus = round(contig_mean, 3),
        mean_prob = round(contig_mean, 3),
        sd_prob = 0,
        n_merged = 1
      ))
    } else if (contig_case == "prophage") {
      custom_log("INFO", paste("  Filtering predictions for prophage contig:", sanitized_header), "3/8")
      contig_preds <- all_preds %>% 
        filter(contig_name == sanitized_header) %>%
        mutate(prediction = ifelse(is_virus >= 0.5, "virus", "non-virus"))
      
      custom_log("INFO", paste("  Found", nrow(contig_preds), "predictions for prophage contig"), "3/8")
      
      # If no predictions for this contig, skip
      if (nrow(contig_preds) == 0) {
        custom_log("WARN", paste("No predictions for prophage contig", sanitized_header, "- treating as non-virus"), "3/8")
        final_output <- rbind(final_output, data.frame(
          contig_name = sanitized_header,
          original_name = original_header,
          start = 1,
          end = contig_len,
          length = contig_len,
          prediction = "non-virus",
          probability_virus = 0,
          mean_prob = 0,
          sd_prob = 0,
          n_merged = 1
        ))
        next
      }
      
      # Identify all positions where prediction changes
      positions <- sort(unique(c(
        1,
        contig_preds$start,
        contig_preds$end + 1,
        contig_len + 1
      )))
      
      # Create segments for analysis
      segments <- data.frame(
        start = positions[-length(positions)],
        end = positions[-1] - 1
      ) %>%
        filter(end >= start) # Remove invalid segments
      
      regions <- data.frame()
      last_pred <- NULL
      last_prob <- NULL
      last_sd <- NULL
      last_n <- NULL
      region_start <- NULL
      
      # For each segment, determine predominant prediction
      for (i in 1:nrow(segments)) {
        seg_start <- segments$start[i]
        seg_end <- segments$end[i]
        
        # Find overlapping predictions and calculate weighted average
        overlaps <- contig_preds %>%
          filter(start <= seg_end & end >= seg_start) %>%
          mutate(
            overlap_start = pmax(start, seg_start),
            overlap_end = pmin(end, seg_end),
            overlap_length = overlap_end - overlap_start + 1
          )
        
        if (nrow(overlaps) == 0) {
          # Use nearest neighbor or default to non-virus if no data
          if (!is.null(last_pred)) {
            curr_pred <- last_pred
            curr_prob <- last_prob
            curr_sd <- 0
            curr_n <- 1
          } else {
            # Default to non-virus for segments with no data
            curr_pred <- "non-virus"
            curr_prob <- 0
            curr_sd <- 0
            curr_n <- 0
          }
        } else {
          # Calculate weighted average for this segment
          weighted_sum <- sum(overlaps$is_virus * overlaps$overlap_length)
          total_length <- sum(overlaps$overlap_length)
          curr_prob <- weighted_sum / total_length
          
          # Calculate weighted standard deviation if we have multiple predictions
          if (nrow(overlaps) > 1) {
            weighted_var <- sum(overlaps$overlap_length * (overlaps$is_virus - curr_prob)^2) / total_length
            curr_sd <- sqrt(weighted_var)
          } else {
            curr_sd <- 0
          }
          
          curr_n <- nrow(overlaps)
          curr_pred <- ifelse(curr_prob >= 0.5, "virus", "non-virus")
        }
        
        # Start a new region or extend current one
        if (is.null(region_start) || curr_pred != last_pred) {
          # Save previous region if it exists
          if (!is.null(region_start)) {
            regions <- rbind(regions, data.frame(
              contig_name = sanitized_header,
              original_name = original_header,
              start = region_start,
              end = seg_start - 1,
              length = seg_start - region_start,
              prediction = last_pred,
              probability_virus = round(last_prob, 3),
              mean_prob = round(last_prob, 3),
              sd_prob = round(last_sd, 3),
              n_merged = last_n
            ))
          }
          # Start new region
          region_start <- seg_start
        } else {
          # For continuing regions, update the standard deviation and count
          if (!is.null(last_sd) && !is.null(last_n) && !is.null(last_prob)) {
            # Combine the statistics for the merged region
            total_len <- (seg_start - region_start) + (seg_end - seg_start + 1)
            weight1 <- (seg_start - region_start) / total_len
            weight2 <- (seg_end - seg_start + 1) / total_len
            
            # Calculate combined mean
            combined_mean <- weight1 * last_prob + weight2 * curr_prob
            
            # Calculate combined variance using parallel axis theorem
            combined_var <- weight1 * (last_sd^2 + (last_prob - combined_mean)^2) + 
                           weight2 * (curr_sd^2 + (curr_prob - combined_mean)^2)
            
            curr_sd <- sqrt(combined_var)
            curr_n <- last_n + curr_n
            curr_prob <- combined_mean
          }
        }
        
        # Update for next iteration
        last_pred <- curr_pred
        last_prob <- curr_prob
        last_sd <- curr_sd
        last_n <- curr_n
      }
      
      # Add the final region
      if (!is.null(region_start)) {
        regions <- rbind(regions, data.frame(
          contig_name = sanitized_header,
          original_name = original_header,
          start = region_start,
          end = seg_end,
          length = seg_end - region_start + 1,
          prediction = last_pred,
          probability_virus = round(last_prob, 3),
          mean_prob = round(last_prob, 3),
          sd_prob = round(last_sd, 3),
          n_merged = last_n
        ))
      }
      
      # Apply prophage size filter - reclassify small viral regions in bacterial background as non-virus
      if (nrow(regions) > 0) {
        # Calculate region sizes
        regions$size <- regions$end - regions$start + 1
        
        # Log the number of regions before filtering
        custom_log("INFO", paste("  Found", sum(regions$prediction == "virus"), "viral regions before size filtering"), "3/8")
        
        # Find bacterial background regions (those classified as non-virus)
        has_bacterial_background <- any(regions$prediction == "non-virus")
        
        # Apply size filter only when there's a bacterial background
        if (has_bacterial_background) {
          # Reclassify viral regions smaller than threshold
          small_viral_regions <- regions$prediction == "virus" & regions$size < prophage_min_size
          
          if (any(small_viral_regions)) {
            custom_log("INFO", paste("  Reclassifying", sum(small_viral_regions), "small viral regions (< ", prophage_min_size, "bp) to non-virus"), "3/8")
            regions$prediction[small_viral_regions] <- "non-virus"
            # Keep probability but set below threshold
            regions$probability_virus[small_viral_regions] <- pmin(regions$probability_virus[small_viral_regions], 0.49)
          }
          
          # Try to merge adjacent non-virus regions
          if (nrow(regions) > 1) {
            merged_regions <- data.frame()
            current_region <- regions[1,]
            
            for (i in 2:nrow(regions)) {
              if (regions$prediction[i] == current_region$prediction && 
                  regions$prediction[i] == "non-virus" &&
                  regions$start[i] == current_region$end + 1) {
                # Merge with current region
                curr_size <- current_region$size
                new_size <- regions$size[i]
                total_size <- curr_size + new_size
                
                # Update region end and size
                current_region$end <- regions$end[i]
                current_region$size <- total_size
                current_region$length <- total_size
                
                # Compute new weighted average probability
                current_region$probability_virus <- 
                  (current_region$probability_virus * curr_size + 
                     regions$probability_virus[i] * new_size) / total_size
                
                # Update mean_prob to match probability_virus
                current_region$mean_prob <- current_region$probability_virus
                
                # Update sd_prob using parallel axis theorem for combining variances
                w1 <- curr_size / total_size
                w2 <- new_size / total_size
                combined_mean <- current_region$probability_virus
                
                # Compute combined variance 
                var1 <- current_region$sd_prob^2
                var2 <- regions$sd_prob[i]^2
                mean1 <- current_region$mean_prob
                mean2 <- regions$mean_prob[i]
                
                combined_var <- w1 * (var1 + (mean1 - combined_mean)^2) + 
                               w2 * (var2 + (mean2 - combined_mean)^2)
                
                current_region$sd_prob <- sqrt(combined_var)
                
                # Update the count of merged regions
                current_region$n_merged <- current_region$n_merged + regions$n_merged[i]
              } else {
                # Add current region to results and start a new one
                merged_regions <- rbind(merged_regions, current_region)
                current_region <- regions[i,]
              }
            }
            # Add the last region
            merged_regions <- rbind(merged_regions, current_region)
            regions <- merged_regions
          }
        }
        
        # Log the number of regions after filtering
        custom_log("INFO", paste("  After filtering:", sum(regions$prediction == "virus"), "viral regions and", 
                              sum(regions$prediction == "non-virus"), "non-viral regions"), "3/8")
        
        # Remove the size column before adding to final output
        regions <- regions %>% select(-size)
      }
      
      # Add to final output
      if (nrow(regions) > 0) {
        final_output <- rbind(final_output, regions)
      } else {
        # Fallback if no regions were created
        final_output <- rbind(final_output, data.frame(
          contig_name = sanitized_header,
          original_name = original_header,
          start = 1,
          end = contig_len,
          length = contig_len,
          prediction = "non-virus",
          probability_virus = 0,
          mean_prob = 0,
          sd_prob = 0,
          n_merged = 1
        ))
      }
    }
  }, error = function(e) {
    custom_log("ERROR", paste("Error processing contig", sanitized_header, ":", e$message), "3/8")
    print(traceback())
  })
}

# For output CSV, map back to original contig names
output_csv <- final_output %>%
  select(-contig_name) %>%
  rename(contig_name = original_name) %>%
  select(contig_name, start, end, length, prediction, probability_virus, mean_prob, sd_prob, n_merged)

# Write output files
output_path <- file.path(output_directory, "results.csv")
write.csv(output_csv, output_path, row.names = FALSE, quote = FALSE)
custom_log("INFO", paste("Results written to:", output_path), "3/8")

# Generate plots with probabilities
custom_log("INFO", "Generating visualization plots", "3/8")
pdf_path <- paste0(output_directory, "/results.pdf")
pdf(pdf_path, width = 20, height = 10)

# First identify which contigs have interesting patterns to plot
contigs_to_plot <- final_output %>%
  group_by(contig_name) %>%
  summarize(
    original_name = first(original_name),
    has_multiple_regions = n() > 1,
    has_virus = any(prediction == "virus"),
    max_prob = max(probability_virus),
    .groups = "drop"
  ) %>%
  filter(has_multiple_regions | has_virus | max_prob > 0.3) %>%
  pull(contig_name)

if (length(contigs_to_plot) == 0) {
  # Create a simple message plot if nothing interesting to show
  plot(1, type="n", axes=FALSE, xlab="", ylab="")
  text(1, 1, "No viral regions detected in any contig", cex=2)
  custom_log("INFO", "No contigs with viral regions to plot", "3/8")
} else {
  custom_log("INFO", paste("Plotting", length(contigs_to_plot), "contigs with interesting patterns"), "3/8")
  
  for (i in seq_along(contigs_to_plot)) {
    sanitized_contig <- contigs_to_plot[i]
    original_contig <- final_output$original_name[final_output$contig_name == sanitized_contig][1]
    
    custom_log("INFO", paste("Plotting contig", i, "of", length(contigs_to_plot), ":", original_contig), "3/8")
    
    contig_data <- final_output %>% 
      filter(contig_name == sanitized_contig) %>%
      arrange(start)
    
    # Create the plot using the original name for display
    p <- ggplot(contig_data, aes(xmin = start, xmax = end, ymin = 0, ymax = probability_virus, fill = prediction)) +
      geom_rect() +
      scale_fill_manual(values = c("non-virus" = "blue", "virus" = "red")) +
      labs(title = paste("Contig:", original_contig),
           subtitle = paste("Virus prediction by region (", i, "of", length(contigs_to_plot), ")", sep=""),
           x = "Position (bp)",
           y = "Probability of Virus",
           fill = "Prediction") +
      theme_minimal() +
      ylim(0, 1)
    
    print(p)
  }
}

# Properly close the PDF device
dev.off()
custom_log("INFO", paste("Plots written to:", pdf_path), "3/8")

# Classify and write FASTA subsets
viral_contigs_sanitized <- final_output %>% filter(prediction == "virus") %>% pull(contig_name) %>% unique()
viral_contigs <- fasta_data$Header[fasta_data$SanitizedHeader %in% viral_contigs_sanitized]

non_viral_contigs_sanitized <- final_output %>% 
  filter(prediction == "non-virus" & !contig_name %in% viral_contigs_sanitized) %>% 
  pull(contig_name) %>% 
  unique()
non_viral_contigs <- fasta_data$Header[fasta_data$SanitizedHeader %in% non_viral_contigs_sanitized]

custom_log("INFO", paste("Number of contigs with viral regions:", length(viral_contigs)), "3/8")
custom_log("INFO", paste("Number of fully non-viral contigs:", length(non_viral_contigs)), "3/8")

# Use Header column for lookup since we're writing the original FASTA format
viral_fasta <- if (length(viral_contigs) > 0) fasta_data[fasta_data$Header %in% viral_contigs, ] else NULL
safe_write_fasta(viral_fasta, file.path(output_directory, "viral_contigs.fasta"))

non_viral_fasta <- if (length(non_viral_contigs) > 0) fasta_data[fasta_data$Header %in% non_viral_contigs, ] else NULL
safe_write_fasta(non_viral_fasta, file.path(output_directory, "non_viral_contigs.fasta"))

# Clean up temporary files
on.exit({
  unlink(temp_file_coarse)
  if (exists("temp_refine_file_5000")) unlink(temp_refine_file_5000)
  if (exists("temp_file_5000")) unlink(temp_file_5000)
  if (exists("temp_fasta_file")) unlink(temp_fasta_file)
  if (exists("temp_file_fine")) unlink(temp_file_fine)
}, add = TRUE)