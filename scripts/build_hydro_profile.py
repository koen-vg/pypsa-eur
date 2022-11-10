#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: : 2017-2022 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Build hydroelectric inflow time-series for each country.

Relevant Settings
-----------------

.. code:: yaml

    countries:

    renewable:
        hydro:
            cutout:
            clip_min_inflow:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`, :ref:`renewable_cf`

Inputs
------

- ``data/EIA_scaled_hydro_1980_2020.csv``: Hydroelectricity net generation per country and year (`EIA <https://www.eia.gov/international/data/world/electricity/electricity-generation?pd=2&p=000000000000000000000080000000g&u=0&f=A&v=mapbubble&a=-&i=none&vo=value&t=G&g=000000000000008&l=71-0028g000017kg6940080a44000e02000g00004g8001o&l=275-00000008010000008000000000008g000g000000g0g000g0004&l=164-10002000000000020g00000000000000000000000000000000000008&l=170--110&s=315532800000&e=1577836800000>`_, calculations are described in ``../../../data/calculations_scaling_hydro_data.ods``)

    .. image:: ../img/hydrogeneration.png
        :scale: 33 %

- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``"cutouts/" + config["renewable"]['hydro']['cutout']``: confer :ref:`cutout`

Outputs
-------

- ``resources/profile_hydro.nc``:

    ===================  ================  =========================================================
    Field                Dimensions        Description
    ===================  ================  =========================================================
    inflow               countries, time   Inflow to the state of charge (in MW),
                                           e.g. due to river inflow in hydro reservoir.
    ===================  ================  =========================================================

    .. image:: ../img/inflow-ts.png
        :scale: 33 %

    .. image:: ../img/inflow-box.png
        :scale: 33 %

Description
-----------

.. seealso::
    :mod:`build_renewable_profiles`
"""

import logging
from tempfile import NamedTemporaryFile

import atlite
import country_converter as coco
import geopandas as gpd
import pandas as pd
import xarray as xr
from _helpers import configure_logging
from vresutils import hydro as vhydro

cc = coco.CountryConverter()


def get_eia_annual_hydro_generation(fn, countries):
    # in billion kWh/a = TWh/a
    df = pd.read_csv(fn, skiprows=2, index_col=1, na_values=[" ", "--"]).iloc[1:, 1:]
    df.index = df.index.str.strip()

    former_countries = {
        "Former Czechoslovakia": dict(
            countries=["Czech Republic", "Slovakia"], start=1980, end=1992
        ),
        "Former Serbia and Montenegro": dict(
            countries=["Serbia", "Montenegro"], start=1992, end=2005
        ),
        "Former Yugoslavia": dict(
            countries=[
                "Slovenia",
                "Croatia",
                "Bosnia and Herzegovina",
                "Serbia",
                "Montenegro",
                "North Macedonia",
            ],
            start=1980,
            end=1991,
        ),
    }

    for k, v in former_countries.items():
        period = [str(i) for i in range(v["start"], v["end"] + 1)]
        ratio = df.loc[v["countries"]].T.dropna().sum()
        ratio /= ratio.sum()
        for country in v["countries"]:
            df.loc[country, period] = df.loc[k, period] * ratio[country]

    baltic_states = ["Latvia", "Estonia", "Lithuania"]
    df.loc[baltic_states] = (
        df.loc[baltic_states].T.fillna(df.loc[baltic_states].mean(axis=1)).T
    )

    df.loc["Germany"] = df.filter(like="Germany", axis=0).sum()
    df.loc["Serbia"] += df.loc["Kosovo"].fillna(0.0)
    df = df.loc[~df.index.str.contains("Former")]
    df.drop(["Europe", "Germany, West", "Germany, East", "Kosovo"], inplace=True)

    df.index = cc.convert(df.index, to="iso2")
    df.index.name = "countries"

    df = df.T[countries] * 1e6  # in MWh/a

    return df


logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_hydro_profile")
    configure_logging(snakemake)

    config_hydro = snakemake.config["renewable"]["hydro"]

    # Load the provided cutout files and merge them.
    cutouts = [atlite.Cutout(c) for c in snakemake.input.cutouts]
    xdatas = [c.data for c in cutouts]
    combined_data = xr.concat(xdatas, dim="time")
    cutout = atlite.Cutout(NamedTemporaryFile().name, data=combined_data)

    countries = snakemake.config["countries"]
    country_shapes = (
        gpd.read_file(snakemake.input.country_shapes)
        .set_index("name")["geometry"]
        .reindex(countries)
    )
    country_shapes.index.name = "countries"

    fn = snakemake.input.eia_hydro_generation
    eia_stats = vhydro.get_eia_annual_hydro_generation(fn).reindex(columns=countries)

    inflow = cutout.runoff(
        shapes=country_shapes,
        smooth=True,
        lower_threshold_quantile=True,
        normalize_using_yearly=eia_stats,
    )

    if "clip_min_inflow" in config_hydro:
        inflow = inflow.where(inflow > config_hydro["clip_min_inflow"], 0)

    inflow.to_netcdf(snakemake.output[0])
