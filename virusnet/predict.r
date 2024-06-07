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
  library(ggplot2)
  library(zoo)
  library(keras)
  library(reticulate)
}))

source("virusnet/utils.r")

#use_python(paste0(Sys.getenv("PREFIX"), "/bin/python"), required = TRUE)
# use_python("/home/pmuench/micromamba/envs/build-env2/conda-bld/virusnet-0_1717753851740/_test_e/bin/python")

args <- commandArgs(trailingOnly = TRUE)

# Parse command line arguments
arg_list <- list()
for (i in seq(1, length(args), 2)) {
  arg_list[[gsub("--", "", args[i])]] <- args[i + 1]
}

print(arg_list)

message("Checking input file")

# Read the FASTA file using microseq
fasta_data <- readFasta(arg_list$input)

# Count the number of FASTA entries
num_fasta_entries <- nrow(fasta_data)

# Print the number of FASTA entries
message(paste("Number of FASTA entries in the file:", num_fasta_entries))

message("Loading genus model...")
model_genus <- keras::load_model_hdf5(arg_list$model_genus, compile = FALSE)

message("Loading binary model...")
model_binary <- keras::load_model_hdf5(arg_list$model_binary, custom_objects = custom_objects)


# predict using genus model
pred <- predict_model(
    output_format = "by_entry_one_file",
    model = model_genus,
    layer_name = "dense_3",
    path_input = arg_list$input,
    round_digits = 3,
    step = 1000,
    batch_size = 5120,
    verbose = FALSE,
    return_states = FALSE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = "output_genus.h5",
    return_int = FALSE)

head(pred)

length(model_binary$input_shape)

pred <- predict_model(
    output_format = "by_entry_one_file",
    model = model_binary,
    layer_name = "dense_26",
    path_input = arg_list$input,
    round_digits = 3,
    step = 1000,
    batch_size = 512,
    verbose = FALSE,
    return_states = FALSE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = "output_bin.h5",
    return_int = length(model_binary$input_shape) == 2)
