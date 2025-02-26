#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -x  # Print commands and their arguments as they are executed

mkdir -p $PREFIX/bin

# Install the Python script
echo "Installing the scripts..."
cp $SRC_DIR/virusnet.py $PREFIX/bin/virusnet
cp $SRC_DIR/utils.r $PREFIX/bin/utils.r
cp $SRC_DIR/setup_logger.r $PREFIX/bin/setup_logger.r
cp $SRC_DIR/models.json $PREFIX/bin/models.json
cp $SRC_DIR/predict_genus.r $PREFIX/bin/predict_genus.r
cp $SRC_DIR/predict_binary.r $PREFIX/bin/predict_binary.r
cp $SRC_DIR/predict_binary_metagenome.r $PREFIX/bin/predict_binary_metagenome.r

# Set the R library directory
mkdir -p "$PREFIX/lib/R/library"
chmod -R u+w "$PREFIX/lib/R/library"
export R_LIBS_SITE="$PREFIX/lib/R/library"

# Install deepG with verbose output and error handling
echo "Installing deepG package..."
"$BUILD_PREFIX/bin/R" --vanilla -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes'); \
  Sys.setenv(R_LIBS_USER='$PREFIX/lib/R/library'); \
  .libPaths('$PREFIX/lib/R/library'); \
  print(.libPaths()); \
  tryCatch({ \
    remotes::install_github('GenomeNet/deepg', lib='$PREFIX/lib/R/library', upgrade='never', force=TRUE, verbose=TRUE); \
    library(deepG); \
    print('deepG installed successfully!') \
  }, error = function(e) { \
    message('Error installing deepG: ', e$message); \
    print(.libPaths()); \
    system('find $PREFIX/lib/R -type d'); \
    quit(status = 1) \
  })"

# Verify installation
"$BUILD_PREFIX/bin/R" --vanilla -e "print(.libPaths()); installed.packages()"

# After installing deepG
"$BUILD_PREFIX/bin/R" --vanilla -e "packagePath <- find.package('deepG'); print(packagePath); if(dir.exists(packagePath)) { system(paste('cp -r', packagePath, '$PREFIX/lib/R/library/')) }" 