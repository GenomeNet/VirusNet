# VirusNet

[![Build and Test Conda Package](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/version.svg)](https://anaconda.org/genomenet/virusnet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/latest_release_relative_date.svg)](https://anaconda.org/genomenet/virusnet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/downloads.svg)](https://anaconda.org/genomenet/virusnet)

## Install

Make sure you have conda installed. If you don't have it installed, you can download and install Miniconda from the official website: https://docs.conda.io/en/latest/miniconda.html

```bash
conda create -n virusnet python=3.11 -y
conda activate virusnet
conda install -c anaconda -c conda-forge -c genomenet virusnet -y
```

## Installation using Mamba (recommended)

```bash
conda install mamba -c conda-forge -y
mamba create -n virusnet python=3.11 -y
mamba activate virusnet
mamba install -c genomenet -c anaconda -c conda-forge genomenet::virusnet -y
```

## Usage

Download the models

```
virusnet download
```

Make predictions


```
virusnet predict --mode binary --input test/covid.fasta --output results_binary
```

This output will be written to the screen

```
[INFO] Checking input file 
[INFO] Number of FASTA entries in the file: 1 
[INFO] Loading binary model 
[INFO] Performing predictions 
[INFO] Using non-metagenomic mode 
[INFO] Summarized results written to: results_binary/binary_results_summarized.csv 
[INFO] Non-summarized results written to: results_binary/binary_results.csv 
[INFO] Number of contigs classified as viral: 1 
[INFO] Number of contigs classified as non-viral: 0 
[INFO] FASTA data written to: results_binary/viral_contigs.fasta 
[WARN] No contigs found 
```

```
virusnet predict --mode genus --input test/covid.fasta --output results_genus
```

Which will produce this output

```
[INFO] Checking input file 
[INFO] Number of FASTA entries in the file: 1 
[INFO] Loading genus model 
[INFO] Performing genus predictions 
Contig Summary:
Contig 1: Virus - Alphacoronavirus, Probability - 37%
```

## Development

### Building this conda package


```
# Create a fresh build environment
micromamba create -n build-env python=3.11 boa anaconda-client r-base=4.3.1 r-remotes -c conda-forge -y
micromamba activate build-env

# Clean any previous builds
rm -rf ~/micromamba/envs/build-env/conda-bld/

# Set environment variables for R
export R_MAX_VSIZE=8Gb
mkdir -p ~/.R
echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' > ~/.Rprofile

# Build the package
boa build . #  writes something like this to the screen /home/pmuench/micromamba/envs/build-env/conda-bld/noarch/virusnet-0.9.5-hadc6f11_0.tar.bz2

# Upload 
anaconda upload --user genomenet --channel genomenet /home/pmuench/micromamba/envs/build-env/conda-bld/noarch/virusnet-0.9.5-hadc6f11_0.tar.bz2

conda index ~/micromamba/envs/build-env/conda-bld/

# Deactivate the build environment
micromamba deactivate

# Clean micromamba cache to avoid any issues
micromamba clean --all -y

# Create a fresh test environment
micromamba create -n virusnet-test python=3.11 -y
micromamba activate virusnet-test

# Method 1: Install from local channel (preferred)
micromamba install /home/pmuench/micromamba/envs/build-env/conda-bld/noarch/virusnet-0.9.5-hadc6f11_0.tar.bz2 -y

virusnet -h
```