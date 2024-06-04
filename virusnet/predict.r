#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# Parse command line arguments
arg_list <- list()
for (i in seq(1, length(args), 2)) {
  arg_list[[gsub("--", "", args[i])]] <- args[i + 1]
}

# Load necessary library
library(keras)

# Function to load model and print summary
load_and_summarize <- function(model_path) {
  model <- keras::load_model_hdf5(model_path, compile = FALSE)
  print(summary(model))
}

# Load and summarize each model
cat("Binary Model Summary:\n")
load_and_summarize(arg_list$model_binary)

cat("\nGenus Model Summary:\n")
load_and_summarize(arg_list$model_genus)
