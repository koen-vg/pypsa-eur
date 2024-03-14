# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023-2024 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

# This script helps to generate a scenarios.yaml file for PyPSA-Eur.
# You can modify the template to your needs and define all possible combinations of config values that should be considered.

if "snakemake" in globals():
    filename = snakemake.output[0]
else:
    filename = "../config/scenarios.yaml"

import itertools

template = """
{scenario_name}:
    snapshots:
        start: "{year}-01-01 00:00"
        end: "{year}-12-31 23:00"
        inclusive: "both"
    renewable:
        onwind:
            cutout: europe-era5_{year}
        solar:
            cutout: europe-era5_{year}
        offwind-ac:
            cutout: europe-era5_{year}
        offwind-dc:
            cutout: europe-era5_{year}
        offwind-float:
            cutout: europe-era5_{year}
        hydro:
            cutout: europe-era5_{year}
    solar_thermal:
        cutout: europe-era5_{year}
    lines:
        dynamic_line_rating:
            cutout: europe-era5_{year}
    sector:
        green_imports: {green_imports}
        co2_sequestration_potential: {co2_sequestration_potential}
        co2_sequestration_cost: {co2_sequestration_cost}
        heat_demand_cutout: europe-era5_{year}
    adjustments:
        sector:
            capital_cost:
                "H2 Electrolysis": {electrolysis_cost}
                "SMR CC": {smr_cc_cost}
                "nuclear": {nuclear_cost}
            p_nom_max:
                "onwind": {land_use}
                "solar": {land_use}
"""

scenarios = {
    # Year
    "Y": {
        "1987": {
            "year": "1987",
        },
        # "2013": {
        #     "year": "2013",
        # },
        "2020": {
            "year": "2020",
        },
    },
    # Electrolysis
    "E": {
        "L": {"electrolysis_cost": 0.5},
        "H": {"electrolysis_cost": 1.5},
    },
    # CCS and adjacent technologies
    "C": {
        "L": {
            # Bring the cost of SMR CC close to the cost of SMR
            # without CC (which is 86.2% the cost of SMR CC by
            # default)
            "smr_cc_cost": 0.9,
            # An optimistic cost of co2 sequestration (assumption) [€/tCO2]
            "co2_sequestration_cost": 15,
            # Source?
            "co2_sequestration_potential": 2000,
        },
        "H": {
            "smr_cc_cost": 1.5,
            "co2_sequestration_cost": 30,  # [€/tCO2]
            "co2_sequestration_potential": 400,
        },
    },
    # Land use
    "L": {
        "L": {"land_use": 0.5},
        "H": {"land_use": 1},
    },
    # Imports
    "I": {
        "L": {
            "green_imports": "false",
        },
        "H": {
            "green_imports": "true",
        },
    },
    # Nuclear
    "N": {
        "D": {
            "nuclear_cost": 1,
        },
        # "L": {
        #     "nuclear_cost": 0.5,
        # },
        # "H": {
        #     "nuclear_cost": 1.5,
        # },
    },
}

# Create set product of all scenario outcomes
scenario_levels = {k: list(v.keys()) for k, v in scenarios.items()}
scenario_levels = [[(k, v) for v in levels] for k, levels in scenario_levels.items()]
combinations = list(itertools.product(*scenario_levels))


def merge_dicts(list_of_dicts):
    return {k: v for d in list_of_dicts for k, v in d.items()}


config_combinations = [
    merge_dicts([scenarios[s][l] for (s, l) in c]) for c in combinations
]

with open(filename, "w") as f:
    for c, config in zip(combinations, config_combinations):
        # Construct name string from c
        name = "_".join([f"{s}-{l}" for s, l in c])
        f.write(template.format(scenario_name=name, **config))
