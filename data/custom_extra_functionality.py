# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


def custom_extra_functionality(n, snapshots, snakemake):
    """
    Add custom extra functionality constraints.
    """
    min_prod_const = 60000000
    n.add(
        "GlobalConstraint",
        name="min_solar_prod",
        sense=">=",
        constant=min_prod_const,
        carrier_attribute="solar",
    )

    p = n.model["Generator-p"]
    solar_i = n.generators.index[n.generators.carrier == "solar"]
    p_solar = p.sel(Generator=solar_i)
    p_solar = p_solar.sum("Generator") * n.snapshot_weightings.generators
    p_solar_total = p_solar.sum()
    n.model.add_constraints(
        p_solar_total >= min_prod_const, name="GlobalConstraint-min_solar_prod"
    )
