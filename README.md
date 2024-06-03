# VirusNet

## Build

```
conda create -n virusnet_build_env python=3.11 -y
conda mambabuild .   
conda install --use-local virusnet -y
virusnet --help # test the package
```