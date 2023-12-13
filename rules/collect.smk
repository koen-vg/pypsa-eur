# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


import itertools
import random
import json
import hashlib


localrules:
    all,
    cluster_networks,
    extra_components_networks,
    prepare_elec_networks,
    prepare_sector_networks,
    solve_elec_networks,
    solve_sector_networks,
    plot_networks,


rule cluster_networks:
    input:
        expand(RESOURCES + "networks/elec_s{simpl}_{clusters}.nc", **config["scenario"]),


rule extra_components_networks:
    input:
        expand(
            RESOURCES + "networks/elec_s{simpl}_{clusters}_ec.nc", **config["scenario"]
        ),


rule prepare_elec_networks:
    input:
        expand(
            RESOURCES + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule prepare_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule solve_elec_networks:
    input:
        expand(
            RESULTS + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule copy_sweep_network:
    output:
        RESULTS + "param_sweeps/{hash_digest}/{net}",
    input:
        RESULTS + "postnetworks/{net}",
    shell:
        """
        mkdir -p "$(dirname '{output}')"
        ln -sr "{input}" "{output}"
        """


def param_sweep(param_sweep, scenario, num_nets):
    def build_opt_strs(conf):
        return [conf["carrier"] + "+" + conf["attr"] + str(f) for f in conf["factors"]]

    networks = list(
        expand(
            "elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}",
            **scenario,
        )
    )

    opts = [build_opt_strs(param_sweep[n]) for n in param_sweep]
    param_space = list(itertools.product(*opts))
    random.seed(0)
    random.shuffle(param_space)
    opts = ["-".join(s) for s in param_space[:num_nets]]

    hash_digest = hashlib.md5(("".join(opts) + "".join(networks)).encode()).hexdigest()[:8]
    
    return [
        os.path.join(RESULTS, "param_sweeps", hash_digest, f"{n}-{o}_{p}.nc")
        for n in networks
        for o in opts
        for p in scenario["planning_horizons"]
    ]


rule solve_sector_networks_sweep:
    input:
        param_sweep(
            config["param_sweep"], config["scenario"], config["param_sweep_num_nets"]
        ),


rule solve_sector_networks_perfect:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc",
            **config["scenario"]
        ),


rule plot_networks:
    input:
        expand(
            RESULTS
            + "maps/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),


rule validate_elec_networks:
    input:
        expand(
            RESULTS
            + "figures/.statistics_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config["scenario"]
        ),
        expand(
            RESULTS
            + "figures/.validation_{kind}_plots_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}",
            **config["scenario"],
            kind=["production", "prices", "cross_border"]
        ),
