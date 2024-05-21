# SPDX-FileCopyrightText: : 2023-4 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


localrules:
    add_existing_baseyear,
    add_brownfield,


rule add_existing_baseyear:
    params:
        baseyear=config_provider("scenario", "planning_horizons", 0),
        sector=config_provider("sector"),
        existing_capacities=config_provider("existing_capacities"),
        costs=config_provider("costs"),
    input:
        car_ages=resources("car_ages.csv"),
        truck_ages=resources("truck_ages.csv"),
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        powerplants=resources("powerplants.csv"),
        busmap_s=resources("busmap_elec_s{simpl}.csv"),
        busmap=resources("busmap_elec_s{simpl}_{clusters}.csv"),
        clustered_pop_layout=resources("pop_layout_elec_s{simpl}_{clusters}.csv"),
        costs=lambda w: resources(
            "costs_{}.csv".format(
                config_provider("scenario", "planning_horizons", 0)(w)
            )
        ),
        cop_soil_total=resources("cop_soil_total_elec_s{simpl}_{clusters}.nc"),
        cop_air_total=resources("cop_air_total_elec_s{simpl}_{clusters}.nc"),
        existing_heating_distribution=resources(
            "existing_heating_distribution_elec_s{simpl}_{clusters}_{planning_horizons}.csv"
        ),
    output:
        RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}{near_opt}.nc",
    wildcard_constraints:
        # TODO: The first planning_horizon needs to be aligned across scenarios
        # snakemake does not support passing functions to wildcard_constraints
        # reference: https://github.com/snakemake/snakemake/issues/2703
        planning_horizons=config["scenario"]["planning_horizons"][0],  #only applies to baseyear
        near_opt=r"(_(min|max)[0-9\.]+)?",
    threads: 1
    resources:
        mem_mb=2000,
    log:
        RESULTS
        + "logs/add_existing_baseyear_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}{near_opt}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/add_existing_baseyear/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}{near_opt}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_existing_baseyear.py"


def input_profile_tech_brownfield(w):
    return {
        f"profile_{tech}": resources(f"profile_{tech}.nc")
        for tech in config_provider("electricity", "renewable_carriers")(w)
        if tech != "hydro"
    }


rule add_brownfield:
    params:
        H2_retrofit=config_provider("sector", "H2_retrofit"),
        H2_retrofit_capacity_per_CH4=config_provider(
            "sector", "H2_retrofit_capacity_per_CH4"
        ),
        threshold_capacity=config_provider("existing_capacities", " threshold_capacity"),
        snapshots=config_provider("snapshots"),
        drop_leap_day=config_provider("enable", "drop_leap_day"),
        carriers=config_provider("electricity", "renewable_carriers"),
    input:
        unpack(input_profile_tech_brownfield),
        simplify_busmap=resources("busmap_elec_s{simpl}.csv"),
        cluster_busmap=resources("busmap_elec_s{simpl}_{clusters}.csv"),
        network=RESULTS
        + "prenetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        network_p=solved_previous_horizon,  #solved network at previous time step
        costs=resources("costs_{planning_horizons}.csv"),
        cop_soil_total=resources("cop_soil_total_elec_s{simpl}_{clusters}.nc"),
        cop_air_total=resources("cop_air_total_elec_s{simpl}_{clusters}.nc"),
    output:
        RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}{near_opt}.nc",
    wildcard_constraints:
        near_opt=r"(_(min|max)[0-9\.]+)?",
    threads: 4
    resources:
        mem_mb=10000,
    log:
        RESULTS
        + "logs/add_brownfield_elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}{near_opt}.log",
    benchmark:
        (
            RESULTS
            + "benchmarks/add_brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}{near_opt}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/add_brownfield.py"


ruleorder: add_existing_baseyear > add_brownfield


rule solve_sector_network_myopic:
    params:
        solving=config_provider("solving"),
        foresight=config_provider("foresight"),
        build_year_agg=config_provider("clustering", "build_year_aggregation"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        sector=config_provider("sector"),
        custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        costs=resources("costs_{planning_horizons}.csv"),
    output:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        config=RESULTS
        + "configs/config.elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.yaml",
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_solver.log",
        memory=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_memory.log",
        python=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_python.log",
    threads: solver_threads
    resources:
        mem_mb=lambda w, attempt: (
            config_provider("solving", "mem_mb")(w) * 1.2 ** (attempt - 1)
        ),
        runtime=lambda w, attempt: (
            config_provider("solving", "runtime", default="6h")(w) * attempt
        ),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    priority: 10
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/solve_network.py"


rule near_opt_myopic:
    params:
        solving=config_provider("solving"),
        near_opt=config_provider("near_opt"),
        planning_horizons=config_provider("scenario", "planning_horizons"),
        sector=config_provider("sector"),
        build_year_agg=config_provider("clustering", "build_year_aggregation"),
        # Not supported yet:
        # custom_extra_functionality=input_custom_extra_functionality,
    input:
        network=RESULTS
        + "prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_{sense}{slack}.nc",
        network_opt=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
    output:
        RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_{sense}{slack}.nc",
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/near_opt_myopic/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_{sense}{slack}_solver.log",
        python=RESULTS
        + "logs/near_opt_myopic/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_{sense}{slack}_python.log",
    threads: solver_threads
    resources:
        mem_mb=lambda w, attempt: (
            config_provider("solving", "mem_mb")(w) * 1.2 ** (attempt - 1)
        ),
        runtime=lambda w, attempt: (
            1.5 * config_provider("solving", "runtime", default="6h")(w) * attempt
        ),
    benchmark:
        (
            RESULTS
            + "benchmarks/solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}_{sense}{slack}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/near_opt_myopic.py"
