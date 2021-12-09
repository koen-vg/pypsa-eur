#!/usr/bin/env python

# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
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

- ``data/bundle/EIA_scaled_hydro_1980_2020.csv``: Hydroelectricity net generation per country and year (`EIA <https://www.eia.gov/international/data/world/electricity/electricity-generation?pd=2&p=000000000000000000000080000000g&u=0&f=A&v=mapbubble&a=-&i=none&vo=value&t=G&g=000000000000008&l=71-0028g000017kg6940080a44000e02000g00004g8001o&l=275-00000008010000008000000000008g000g000000g0g000g0004&l=164-10002000000000020g00000000000000000000000000000000000008&l=170--110&s=315532800000&e=1577836800000>`_, calculations are described in ``../../../data/calculations_scaling_hydro_data.ods``)

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
from _helpers import configure_logging

import atlite
import geopandas as gpd
from vresutils import hydro as vhydro

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_hydro_profile')
    configure_logging(snakemake)

    config = snakemake.config['renewable']['hydro']
    cutout = atlite.Cutout(snakemake.input.cutout)

    countries = snakemake.config['countries']
    country_shapes = (gpd.read_file(snakemake.input.country_shapes)
                      .set_index('name')['geometry'].reindex(countries))
    country_shapes.index.name = 'countries'

    eia_stats = vhydro.get_eia_annual_hydro_generation(
        snakemake.input.eia_hydro_generation).reindex(columns=countries)
    inflow = cutout.runoff(shapes=country_shapes,
                           smooth=True,
                           lower_threshold_quantile=True,
                           normalize_using_yearly=eia_stats)

    if 'clip_min_inflow' in config:
        inflow = inflow.where(inflow > config['clip_min_inflow'], 0)

    inflow.to_netcdf(snakemake.output[0])
