#!/bin/bash
# run.sh
echo "hello World"
python setup.py install --single-version-externally-managed --record=record.txt

# Install the randomForest R package
R -e "install.packages('randomForest', repos='https://cran.rstudio.com/')"
R -e "install.packages('keras', repos='https://cran.rstudio.com/')"