# Introduction

This is a modification to PyPSA-Eur which allows running the model with weather- and demand years from 1980 to 2015.

# Modifications

## Load data

Simple country-level load (demand) data is generated for the period based on temperature time series, and fed directly into the PyPSA model.

## Renewable capacity factor time series

We use the ERA5 dataset over the period to generate renewable capacity factor time series.

## Implementation details

We add a `year` wildcard almost throughout the Snakemake workflow, which selects the load and weather year. Nothing else (e.g. technology costs) are changed. The names of cutouts are changed from `europe-2013-era5` to `europe-era5_{year}` (where `{year}` is any year between 1980 and 2015) to facilitate the new wildcard. We exclusively use the ERA5 dataset, not the SARAH dataset.
