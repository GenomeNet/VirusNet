#!/bin/bash

mkdir -p $PREFIX/bin

# Install the Python script
echo "Installing the scripts..."
cp $SRC_DIR/virusnet.py $PREFIX/bin/virusnet  # copy it to the bin directory
cp $SRC_DIR/predict.r $PREFIX/bin/predict.r  # copy it to the bin directory

# Install the R package
echo "Installing the R package from GitHub..."
R -e "devtools::install_github('genomenet/deepg')"