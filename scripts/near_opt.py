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
    **kwargs,
):
    # The following solver setup follows that of `solve_network` in `solve_network.py`:
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

    kwargs["multi_investment_periods"] = True
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

    # TODO: deal graceful with failed solve.

    # We are now ready to really start with near-opt computations.
    num_parallel = near_opt_config.get("num_parallel", 1)
    results = []
    networks = []
    with get_context("spawn").Pool(num_parallel) as pool:
        for slack, sense in itertools.product(
            near_opt_config.get("slacks", [0.05]), ["min", "max"]
        ):
            weights = {}
            static = near_opt_config["weights"].get("static", {})
            for c in static:
                vars = {}
                for v in static[c]:
                    w = pd.Series(0, index=n.df(c).index)
                    for carrier, const in static[c][v].items():
                        w.loc[n.df(c).carrier == carrier] = const
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

            m = n.copy()
            # Add to network for extra_functionality
            m.config = config
            m.opts = opts

            networks.append(m)

            args = dict(
                weights=weights,
                slack=slack,
                sense=sense,
                model_kwargs=model_kwargs,
                **kwargs,
            )
            results.append(
                pool.apply_async(
                    pypsa.optimization.abstract.optimize_mga,
                    (m,),
                    kwds=args,
                    error_callback=print,
                )
            )
            # pypsa.optimization.abstract.optimize_mga(m, **args)

        for r in results:
            r.get()

    return networks


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

    # Add objective back to n from meta (currently not
    # imported/exported)
    n.objective = n.meta["objective"]
    n.objective_constant = n.meta["objective_constant"]

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        networks = near_opt(
            n,
            config=snakemake.config,
            solving=snakemake.params.solving,
            opts=opts,
            near_opt_config=snakemake.params.near_opt_config,
            log_fn=snakemake.log.solver,
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    # Make sure the output directory exists.
    if not os.path.exists(snakemake.output[0]):
        os.makedirs(snakemake.output[0])

    for i, m in enumerate(networks):
        m.export_to_netcdf(os.path.join(snakemake.output[0], f"{i}.nc"))
