# VirusNet

[![Build and Test Conda Package](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml)

## Build

```
mamba create -n build-env python=3.11 -y
micromamba activate build-env
mamba install boa -y
boa build .
```
