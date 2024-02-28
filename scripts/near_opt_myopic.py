import itertools
import logging
import os
from multiprocessing import get_context

import numpy as np
import pandas as pd
import pypsa
from _benchmark import memory_logger
from _helpers import configure_logging, get_opt, update_config_with_sector_opts
from solve_network import extra_functionality, prepare_network

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def near_opt(
    n,
    config,
    solving,
    opts,
    near_opt_config,
    sense,
    slack,
    **kwargs,
):
    # The following solver setup follows that of `solve_network` in `solve_network.py`:
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    kwargs["multi_investment_periods"] = config["foresight"] == "perfect"
    kwargs["solver_options"] = (
        solving["solver_options"][set_of_options] if set_of_options else {}
    )
    kwargs["solver_name"] = solving["solver"]["name"]
    kwargs["extra_functionality"] = extra_functionality
    kwargs["assign_all_duals"] = cf_solving.get("assign_all_duals", False)
    kwargs["io_api"] = cf_solving.get("io_api", None)

    model_kwargs = {}
    model_kwargs["transmission_losses"] = cf_solving.get("transmission_losses", False)
    model_kwargs["linearized_unit_commitment"] = cf_solving.get(
        "linearized_unit_commitment", False
    )

    if kwargs["solver_name"] == "gurobi":
        logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

    n.config = config
    n.opts = opts

    weights = {}
    static = near_opt_config["weights"].get("static", {})
    for c in static:
        vars = {}
        for v in static[c]:
            w = pd.Series(0, index=n.df(c).index)
            for carrier, const in static[c][v].items():
                w.loc[(n.df(c).carrier == carrier) & n.df(c).p_nom_extendable] = const
            vars[v] = w
        weights[c] = vars
    varying = near_opt_config["weights"].get("varying", {})
    for c in varying:
        vars = {}
        for v in varying[c]:
            w = pd.DataFrame(0, columns=n.df(c).index, index=n.snapshots)
            for carrier, const in varying[c][v].items():
                w.loc[:, n.df(c).carrier == carrier] = const
            vars[v] = w
        weights[c] = vars

    status, condition = pypsa.optimization.abstract.optimize_mga(
        n,
        weights=weights,
        slack=slack,
        sense=sense,
        model_kwargs=model_kwargs,
        **kwargs,
    )

    if status == "ok":
        print("Solved successfully")
        return n
    else:
        breakpoint()
        raise ValueError(
            f"Solve with condition status {status} and condition {condition}"
        )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "solve_sector_network",
            configfiles="../config/test/config.perfect.yaml",
            simpl="",
            opts="",
            clusters="37",
            ll="v1.0",
            sector_opts="CO2L0-1H-T-H-B-I-A-dist1",
            planning_horizons="2030",
        )

    # Follow `solve_network.py` for how to set up logging,
    # configuration, etc.
    configure_logging(snakemake)
    if "sector_opts" in snakemake.wildcards.keys():
        update_config_with_sector_opts(
            snakemake.config, snakemake.wildcards.sector_opts
        )

    opts = snakemake.wildcards.opts
    if "sector_opts" in snakemake.wildcards.keys():
        opts += "-" + snakemake.wildcards.sector_opts
    opts = [o for o in opts.split("-") if o != ""]
    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)

    # Note that the prepare_network function is not needed here; it's
    # good enough that this was run before the system cost
    # optimisation. Re-running could lead to land use constraint
    # issues.

    # NB: We _don't_ support custom extra functoinality in this script!
    # n.custom_extra_functionality = snakemake.params.custom_extra_functionality

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        n = near_opt(
            n,
            snakemake.config,
            snakemake.params.solving,
            opts,
            snakemake.params.near_opt,
            snakemake.wildcards.sense,
            float(snakemake.wildcards.slack),
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.export_to_netcdf(snakemake.output[0])
