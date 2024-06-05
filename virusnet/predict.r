#!/usr/bin/env Rscript

# Suppress TF messages
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)
# Set library paths to the Conda environment's R library path
.libPaths(c(paste0(Sys.getenv("PREFIX"), "/lib/R/library")))

# Load libraries
suppressWarnings(suppressPackageStartupMessages({
  library(deepG)
}))

suppressPackageStartupMessages({
  library(magrittr)
  library(optparse)
  library(ggplot2)
  library(zoo)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse command line arguments
arg_list <- list()
for (i in seq(1, length(args), 2)) {
  arg_list[[gsub("--", "", args[i])]] <- args[i + 1]
}

cat("\nGenus Model Summary:\n")
model_genus <- keras::load_model_hdf5(arg_list$model_genus, compile = FALSE)
summary(model_genus)

# Load and summarize each model
cat("Binary Model Summary:\n")
model_binary <- deepG::load_cp(arg_list$model_binary)
summary(model_binary)

