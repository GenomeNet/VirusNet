# VirusNet

[![Build and Test Conda Package](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml)

## Build

```
conda create -n virusnet_build_env python=3.11 -y
conda mambabuild .   
conda install --use-local virusnet -y
virusnet --help # test the package
```
