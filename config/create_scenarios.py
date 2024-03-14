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

# Insert your config values that should be altered in the template.
# Change `config_section` and `config_section2` to the actual config sections.
template = """
run{scenario_number}:
    snapshots:
        start: {snapshot_start}
        end: {snapshot_end}
    sector:
        co2_sequestration_potential: {co2_sequestration_potential}
    adjustments:
        sector:
            capital_cost:
                "H2 Electrolysis": {electrolysis_cost}
                "SMR CC": {smr_cc_cost}
                "SMR": {smr_cost}
                "nuclear": {nuclear_cost}
            marginal_cost:
                "gas": {gas_marginal_cost}
"""

scenarios = {
    "weather_year": {
        "easy": {
            "snapshot_start": "2020-01-01",
            "snapshot_end": "2021-01-01",
        },
        "difficult": {
            "snapshot_start": "1987-01-01",
            "snapshot_end": "1988-01-01",
        },
    },
    "electrolysis": {
        "low": {"electrolysis_cost": 0.5},
        "high": {"electrolysis_cost": 2},
    },
    "CCS": {
        "low": {
            "smr_cc_cost": 0.5,
            "smr_cost": 0.5,
            "co2_sequestration_potential": 2000,
        },
        "high": {
            "smr_cc_cost": 2,
            "smr_cost": 2,
            "co2_sequestration_potential": 400,
        },
    },
    "nuclear": {
        "low": {
            "nuclear_cost": 0.5,
        },
        "high": {
            "nuclear_cost": 2,
        },
    },
    "gas": {
        "low": {
            "gas_marginal_cost": 0.5,
        },
        "high": {
            "gas_marginal_cost": 2,
        },
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
    for i, config in enumerate(config_combinations):
        f.write(template.format(scenario_number=i, **config))
