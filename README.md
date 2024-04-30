<!--
SPDX-FileCopyrightText: 2024 Koen van Greevenbroek
SPDX-License-Identifier: CC-BY-4.0
-->

This repository contains the code and data required to reproduce the results of the paper *Can Norway save the European Union's hydrogen ambition for 2030?* by C.S.W. Cheng, K. van Greevenbroek & I. Viole.

## Instructions

In order to reproduce the main results, clone this repository and follow the instructions given at https://pypsa-eur.readthedocs.io/en/latest/installation.html (starting at installing the Python dependencies).

Instead of `environment.yaml`, use `hydrogen-exports.yaml`, i.e.
```bash 
mamba env create -f envs/hydrogen-exports.yaml
```
Activate the environment by running `conda activate hydrogen-exports`.

Execute the main workflow for Scenario A by running
```bash
snakemake --use-conda --configfile config/config.scenario_A.yaml -call -- solve_sector_networks_sweep
```
Replace `scenario_A` by `scenario_B` and `scenario_B` to run the respective scenarios.

Use the `config.regional-variation.yaml` configuration to generate the results for Figure 8.

Our analysis is contained in the Jupyter notebooks in the `notebooks` directory. Use the `envs/notebook-analysis.yaml` conda environment to run them.

## License

The present repository is based on [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur) which is licensed as free software, primarily under the MIT license; see https://pypsa-eur.readthedocs.io/en/latest/licenses.html for more information and licenses on the data included with PyPSA-Eur.

All changes in code made to PyPSA-Eur in this repository are licensed under the MIT license. This README as well as the figures produced in the Jupyter notebooks are licensed under CC-BY-4.0.
