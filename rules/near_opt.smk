def near_opt_threads(w):
    return solver_threads(w) * config["near_opt"].get("num_parallel", 1)


near_opt_mem = config["solving"]["mem"] * config["near_opt"].get("num_parallel", 1)


rule near_opt_perfect:
    params:
        near_opt_config=config["near_opt"],
        solving=config["solving"],
    input:
        network=RESULTS
        + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_brownfield_all_years.nc",
        config=RESULTS + "config.yaml",
    output:
        directory(
            RESULTS + "near_opt/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}"
        ),
    threads: near_opt_threads
    resources:
        mem_mb=near_opt_mem,
    shadow:
        "shallow"
    log:
        solver=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_near_opt_solver.log",
        python=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_near_opt_python.log",
        memory=RESULTS
        + "logs/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_near_opt_memory.log",
    benchmark:
        (
            BENCHMARKS
            + "solve_sector_network/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_near_opt"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/near_opt.py"


rule all_near_opt_perfect:
    input:
        expand(
            RESULTS + "near_opt/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}",
            **config["scenario"]
        ),


rule all_near_opt_myopic:
    input:
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_"
            + str(config["scenario"]["planning_horizons"][-1])
            + "_max{slack}.nc",
            **config["scenario"]
        ),
        expand(
            RESULTS
            + "postnetworks/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_"
            + str(config["scenario"]["planning_horizons"][-1])
            + "_min{slack}.nc",
            **config["scenario"]
        ),
