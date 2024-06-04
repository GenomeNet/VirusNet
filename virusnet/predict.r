#!/usr/bin/env Rscript

library(keras)

model <- keras::load_model_hdf5("~/.virusnet/virus_genus_2023-01-23.hdf5", compile = FALSE)
summary(model)