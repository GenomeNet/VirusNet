#!/bin/bash

mkdir -p $PREFIX/bin

# Install the Python script
echo "Installing the scripts..."
cp $SRC_DIR/virusnet.py $PREFIX/bin/virusnet
cp $SRC_DIR/predict.r $PREFIX/bin/predict.r
cp $SRC_DIR/utils.r $PREFIX/bin/utils.r


# Set the R library directory
mkdir -p "$PREFIX/lib/R/library"
chmod -R u+w "$PREFIX/lib/R/library"
R_LIBS_SITE="$PREFIX/lib/R/library"
"$BUILD_PREFIX/bin/R" -e "remotes::install_github('GenomeNet/deepg', lib='$PREFIX/lib/R/library')" 