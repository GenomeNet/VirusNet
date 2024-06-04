#!/bin/bash

mkdir -p $PREFIX/bin

# Install the Python script
echo "Installing the scripts..."
cp $SRC_DIR/virusnet.py $PREFIX/bin/virusnet  # copy it to the bin directory
cp $SRC_DIR/predict.r $PREFIX/bin/predict.r  # copy it to the bin directory

echo "R is located at: $(which R)"


# Output R library paths using .libPaths()
echo "R library paths:"
$(which R) -e ".libPaths()"

# Output the R_LIBS and R_LIBS_USER environment variables
echo "R_LIBS: $R_LIBS"
echo "R_LIBS_USER: $R_LIBS_USER"

# Set the R library directory
mkdir -p "$PREFIX/lib/R/library"
chmod -R u+w "$PREFIX/lib/R/library"
R_LIBS_SITE="$PREFIX/lib/R/library"

"$BUILD_PREFIX/bin/R" -e "remotes::install_github('GenomeNet/deepg', lib='$PREFIX/lib/R/library')" # .libPaths() should be $BUILD_PREFIX/lib/R/library 


