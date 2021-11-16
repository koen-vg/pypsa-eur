# Simple analysis of representative days using pypsa-eur

We use `pypsa-eur` to get a coherent set of input time series (demand and capacity factors for wind and PV) for a model of the European energy system with relatively high spatial resolution. Then we divide the time series into days and cluster these days based on the Euclidean distance, to a varying number of clusters. The sum of the distances from days to their cluster centrers gives an indication of how well the limited number of cluster centrers ("representative days") represent the full input.

## Running the analysis

Set up the `pypsa-eur` conda environment (see the [pypsa-eur installation instructions](https://pypsa-eur.readthedocs.io/en/latest/installation.html#install-python-dependencies)) and activate it. Build the required network (without optimising it - we are only interested in the input time series) by running
```sh
snakemake -j<cores> prepare_all_networks
```
Replace `<cores>` by the number of CPU cores you want to use for the build.

Now, run the notebook `notebooks/analysis.ipynb` (again using the `pypsa-eur` conda environment) to generate the figures in `figures`.
