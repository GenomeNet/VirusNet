#!/usr/bin/env Rscript

# Suppress TF messages
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)

# Set library paths to the Conda environment's R library path
#.libPaths(c(paste0(Sys.getenv("PREFIX"), "/lib/R/library")))

# Load libraries
suppressWarnings(suppressPackageStartupMessages({
  library(deepG)
  library(magrittr)
  library(microseq)
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(ggplot2)
  library(keras)
  library(reticulate)
  library(futile.logger)
}))

source("virusnet/utils.r")
source("virusnet/setup_logger.r")

#use_python(paste0(Sys.getenv("PREFIX"), "/bin/python"), required = TRUE)
# use_python("/home/pmuench/micromamba/envs/build-env2/conda-bld/virusnet-0_1717763173756/_test_e/bin/python")

args <- commandArgs(trailingOnly = TRUE)

# Parse command line arguments
arg_list <- list()
for (i in seq(1, length(args), 2)) {
  arg_list[[gsub("--", "", args[i])]] <- args[i + 1]
}

# Check if the output directory exists, if not, create it
output_dir <- arg_list$output

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

custom_log("INFO", "Checking input file", "1/8")
# Read the FASTA file using microseq
fasta_data <- readFasta(arg_list$input)

# Count the number of FASTA entries
num_fasta_entries <- nrow(fasta_data)

# Print the number of FASTA entries
message(paste("Number of FASTA entries in the file:", num_fasta_entries))

custom_log("INFO", "Loading binary model", "2/8")

#  arg_list$model_binary <- "~/.virusnet/pmuench_transfer_learning_virus_bert_5_85.hdf5"
model_binary <- keras::load_model_hdf5(arg_list$model_binary, custom_objects = custom_objects)

custom_log("INFO", "Perfomring predictions", "3/8")

temp_file_genus <- tempfile(fileext = ".h5")
temp_file_binary <- tempfile(fileext = ".h5")

batch_sizes <- c(100000, 10000, 5000, 2500, 1000, 512)  # Define potential batch sizes
success <- FALSE

for (batch_size in batch_sizes) {
  message(paste("Trying batch size:", batch_size))  # Debugging output
  result <- tryCatch({
    sink("/dev/null")
    pred <- predict_model(
      output_format = "by_entry_one_file",
      model = model_binary,
      layer_name = "dense_26",
      path_input = arg_list$input,
      round_digits = 3,
      step = 1000,
      batch_size = batch_size,
      verbose = TRUE,  # Set verbose to TRUE for more detailed output
      return_states = FALSE,
      padding = "standard",
      mode = "label",
      format = "fasta",
      filename = temp_file_binary,
      return_int = length(model_binary$input_shape) == 2
    )
    # Restore normal output
    sink()
    success <- TRUE
    message(paste("Successful prediction with batch size", batch_size))
    # No need to return, just set success to TRUE
  }, error = function(e) {
    message(paste("Error with batch size", batch_size, ": ", e$message))
    # No need to return, just continue to the next iteration
  })
  
  if (success) {
    break  # Exit the loop if prediction succeeds
  }
}

if (!success) {
  message("Failed to predict with any of the batch sizes")
}

# Interprete the output file
df <- load_prediction(h5_path = temp_file_binary, get_sample_position = TRUE, verbose = FALSE)

number_of_processed_contigs <- length(df)
contig_names <- fasta_data$Header

all_contig_states <- list()

# Interprete the output file
df <- load_prediction(h5_path = temp_file_binary, get_sample_position = TRUE, verbose = FALSE)
number_of_processed_contigs <- length(df)
for (i in 1:number_of_processed_contigs) {
  contig_states <- as.data.frame(df[[i]]$states)
  contig_states$sample_end_position <- df[[i]]$sample_end_position
  colnames(contig_states) <- c("is_non_virus", "is_virus", "position")
  contig_states$contig <- rep(i, nrow(contig_states))
  all_contig_states[[i]] <- contig_states
}

# Combine all data frames into one
final_df <- do.call(rbind, all_contig_states)
# Map contig names to the final data frame
final_df$contig_name <- contig_names[final_df$contig]

# Summarize data for each contig_name
contig_summary <- final_df %>%
  group_by(contig_name) %>%
  summarise(
    mean_is_virus = mean(is_virus),
    median_is_virus = median(is_virus),
    iqr_is_virus = IQR(is_virus),
    sd_is_virus = sd(is_virus),
    number_of_entries = n(),
    is_virus_binary = mean(is_virus) >= 0.5  # This will now return TRUE/FALSE
  ) %>%
  ungroup()  # Ensure the result is no longer grouped

contig_summary_df <- as.data.frame(contig_summary)

print(contig_summary_df)


# Define file paths for the summarized and non-summarized data
summarized_file_path <- file.path(output_dir, "binary_results_summarized.csv")
non_summarized_file_path <- file.path(output_dir, "binary_results.csv")


# Write the summarized data frame to CSV
write.csv(contig_summary_df, summarized_file_path, row.names = FALSE, quote = FALSE)

# Assuming 'final_df' is your non-summarized dataset
write.csv(final_df, non_summarized_file_path, row.names = FALSE, quote = FALSE)

# Optional: Message to indicate files have been written
message(paste("Summarized results written to:", summarized_file_path))
message(paste("Non-summarized results written to:", non_summarized_file_path))

# Open a PDF device to save the plots
pdf(paste0(output_dir, "/binary_results.pdf"), width =20, height = 10)

# Read the data from the non-summarized CSV file
plot_data <- final_df

# Generate a plot for each unique contig
unique_contigs <- unique(plot_data$contig_name)
for (contig_name in unique_contigs) {
  contig_data <- plot_data[plot_data$contig_name == contig_name,]
  p <- ggplot(contig_data, aes(x = position, y = is_virus)) +
    geom_line() +
    ggtitle(paste("Contig:", contig_name)) +
    xlab("Position") +
    ylab("Is Non-Virus Probability")
  p <- p + theme_classic()
  p <- p + ylim(0, 1)
  p <- p + geom_hline(yintercept =  0.5, linetype = 2)
  print(p)
}
# Close the PDF device
dev.off()


# Filter contig IDs where is_virus_binary is TRUE
viral_contigs <- contig_summary_df %>%
  filter(is_virus_binary == TRUE) %>%
  pull(contig_name)

# Filter contig IDs where is_virus_binary is FALSE
non_viral_contigs <- contig_summary_df %>%
  filter(is_virus_binary == FALSE) %>%
  pull(contig_name)

# Output summary to the screen
message(paste("Number of contigs classified as viral:", length(viral_contigs)))
message(paste("Number of contigs classified as non-viral:", length(non_viral_contigs)))

# Read the original FASTA file
fasta_data <- readFasta(arg_list$input)

# Function to safely write FASTA if not empty
safe_write_fasta <- function(fasta_subset, file_path) {
  if (!is.null(fasta_subset) && nrow(fasta_subset) > 0) {
    writeFasta(fasta_subset, file_path)
    message(paste("FASTA data written to:", file_path))
  } else {
    message(paste("No contigs found for writing to:", file_path))
  }
}

# Subset and write FASTA data for viral contigs
viral_fasta <- if (length(viral_contigs) > 0) fasta_data[fasta_data$Header %in% viral_contigs,] else NULL
safe_write_fasta(viral_fasta, file.path(output_dir, "viral_contigs.fasta"))

# Subset and write FASTA data for non-viral contigs
non_viral_fasta <- if (length(non_viral_contigs) > 0) fasta_data[fasta_data$Header %in% non_viral_contigs,] else NULL
safe_write_fasta(non_viral_fasta, file.path(output_dir, "non_viral_contigs.fasta"))

head(contig_summary_df)

custom_log("INFO", "Loading genus model", "4/8")
# arg_list$model_genus <- "~/.virusnet/virus_genus_2023-01-23.hdf5"
model_genus <- keras::load_model_hdf5(arg_list$model_genus, compile = FALSE)

custom_log("INFO", "Performing genus predictions\n", "5/8")

# predict using genus model
if (length(viral_contigs) > 0) {
  pred <- predict_model(
    output_format = "by_entry_one_file",
    model = model_genus,
    layer_name = "dense_3",
    path_input = file.path(output_dir, "viral_contigs.fasta"),
    round_digits = 3,
    step = 1000,
    batch_size = batch_size,
    verbose = FALSE,
    return_states = FALSE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = temp_file_genus,
    return_int = length(model_genus$input_shape) == 2)
}
genus_labels <- readRDS(arg_list$genus_labels)
df <- load_prediction(h5_path = temp_file_genus, get_sample_position = TRUE, verbose = FALSE)
df <- as.data.frame(df)
names(df) <- genus_labels

write.csv(df, paste0(output_dir, "/genus_output.csv"), row.names = FALSE, quote = FALSE)
