# Installation
Installation of development version requires `conda` or `mamba` for package management and uses a `setuptools` development install. Tested on `osx-arm64` and `linux-64`. For full fuctionality on `osx-arm64`, `sra-tools` package must be installed separately via brew (bioconda does not currently have arm-complied `sra-tools`). Following tutorials requires addition of the conda environment to jupyter.
1. Download package from github and navigate to the top level folder.
2. Create `mtbvartools` environment. Ideally, run this command with `mamba` not `conda` as the number of dependencies is large and benefits from mamba's speed - it can still take an hour or two regardless.
    ```
    mamba env create -f linux_conda_spec.yml
    ```
3. Add development version of `mtbvartools` package to the environment:
    ```
    conda activate mtbvartools && pip install -e . --config-settings editable_mode=strict
    ```
4. To follow tutorials or generally work with environment in jupyter, add environment as a kernel in an existing `jupyter` install:
    ```
    python -m ipykernel install --user --name mtbvartools --display-name "mtbvartools"
    ```
