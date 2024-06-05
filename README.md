# VirusNet

[![Build and Test Conda Package](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/version.svg)](https://anaconda.org/genomenet/virusnet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/latest_release_relative_date.svg)](https://anaconda.org/genomenet/virusnet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/downloads.svg)](https://anaconda.org/genomenet/virusnet)


## Usage

Download the models

```
virusnet download
```

Make predictions

```
virusnet predict --input file.fasta
```

## Development

### Building this conda package

```
mamba create -n build-env python=3.11 -y
micromamba activate build-env
mamba install boa -y
boa build .
```
