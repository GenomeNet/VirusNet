library(reticulate)
message("Hello World")
# Get the path to the active conda environment
conda_prefix <- Sys.getenv("CONDA_PREFIX")

# Construct the path to the Python executable
python_path <- file.path(conda_prefix, "bin", "python")

# Use the Python executable from the conda environment
reticulate::use_python(python_path)

library(keras)
library(tensorflow)
model_path <- file.path("virusnet/virus_genus_2023-01-23.hdf5")
model <- keras::load_model_hdf5(model_path, compile = FALSE)
summary(model)