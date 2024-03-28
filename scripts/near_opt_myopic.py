import itertools
import logging
import os

import numpy as np
import pandas as pd
import pypsa
from _benchmark import memory_logger
from _helpers import (
    configure_logging,
    set_scenario_config,
    update_config_from_wildcards,
)
from linopy import LinearExpression, QuadraticExpression, merge
from pypsa.descriptors import nominal_attrs
from solve_network import extra_functionality, prepare_network

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def optimize_mga_fixed_bound(
    n,
    obj_bound,
    weights,
    sense="min",
    model_kwargs={},
    **kwargs,
):
    """
    Run modelling-to-generate-alternatives (MGA) on network to find
    near-optimal solutions. Modification of
    pypsa.optimize.abstract.optimize_mga with direct cost bound
    argument.

    Parameters
    ----------
    n : pypsa.Network
    obj_bound : float
        Right hand side on total system cost constraint.
    weights : dict-like
        Weights for alternate objective function. The default is None, which
        minimizes generation capacity. The weights dictionary should be keyed
        with the component and variable (see ``pypsa/variables.csv``), followed
        by a float, dict, pd.Series or pd.DataFrame for the coefficients of the
        objective function. Examples:

        >>> {"Generator": {"p_nom": 1}}
        >>> {"Generator": {"p_nom": pd.Series(1, index=n.generators.index)}}
        >>> {"Generator": {"p_nom": {"gas": 1, "coal": 2}}}
        >>> {"Generator": {"p": pd.Series(1, index=n.generators.index)}
        >>> {"Generator": {"p": pd.DataFrame(1, columns=n.generators.index, index=n.snapshots)}

        Weights for non-extendable components are ignored. The dictionary does
        not need to provide weights for all extendable components.
    sense : str|int
        Optimization sense of alternate objective function. Defaults to 'min'.
        Can also be 'max'.
    model_kwargs: dict
        Keyword arguments used by `linopy.Model`, such as `solver_dir` or
        `chunk`.
    **kwargs:
        Keyword argument used by `linopy.Model.solve`, such as `solver_name`,

    Returns
    status : str
        The status of the optimization, either "ok" or one of the codes listed
        in https://linopy.readthedocs.io/en/latest/generated/linopy.constants.SolverStatus.html
    condition : str
        The termination condition of the optimization, either
        "optimal" or one of the codes listed in
        https://linopy.readthedocs.io/en/latest/generated/linopy.constants.TerminationCondition.html
    """
    if weights is None:
        weights = dict(Generator=dict(p_nom=pd.Series(1, index=n.generators.index)))

    # create basic model
    m = n.optimize.create_model(
        snapshots=n.snapshots,
        **model_kwargs,
    )

    # build budget constraint
    fixed_cost = n.statistics.installed_capex().sum()
    objective = m.objective
    if not isinstance(objective, (LinearExpression, QuadraticExpression)):
        objective = objective.expression

    m.add_constraints(objective + fixed_cost <= obj_bound, name="budget")

    # parse optimization sense
    if (
        isinstance(sense, str)
        and sense.startswith("min")
        or isinstance(sense, int)
        and sense > 0
    ):
        sense = 1
    elif (
        isinstance(sense, str)
        and sense.startswith("max")
        or isinstance(sense, int)
        and sense < 0
    ):
        sense = -1
    else:
        raise ValueError(f"Could not parse optimization sense {sense}")

    # build alternate objective
    objective = []
    for c, attrs in weights.items():
        for attr, coeffs in attrs.items():
            if isinstance(coeffs, dict):
                coeffs = pd.Series(coeffs)
            if attr == nominal_attrs[c] and isinstance(coeffs, pd.Series):
                coeffs = coeffs.reindex(n.get_extendable_i(c))
                coeffs.index.name = ""
            elif isinstance(coeffs, pd.Series):
                coeffs = coeffs.reindex(columns=n.df(c).index)
            elif isinstance(coeffs, pd.DataFrame):
                coeffs = coeffs.reindex(columns=n.df(c).index, index=n.snapshots)
            objective.append(m[f"{c}-{attr}"] * coeffs * sense)

    m.objective = merge(objective)

    status, condition = n.optimize.solve_model(**kwargs)

    # write MGA coefficients into metadata
    n.meta["obj_bound"] = obj_bound
    n.meta["sense"] = sense

    def convert_to_dict(obj):
        if isinstance(obj, pd.DataFrame):
            return obj.to_dict(orient="list")
        elif isinstance(obj, pd.Series):
            return obj.to_dict()
        elif isinstance(obj, dict):
            return {k: convert_to_dict(v) for k, v in obj.items()}
        else:
            return obj

    n.meta["weights"] = convert_to_dict(weights)

    return status, condition


def near_opt(
    n,
    config,
    solving,
    near_opt_config,
    sense,
    obj_base,
    slack,
    **kwargs,
):
    # The following solver setup follows that of `solve_network` in `solve_network.py`:
    set_of_options = solving["solver"]["options"]
    cf_solving = solving["options"]

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

    if "model_options" in solving:
        model_kwargs = model_kwargs | solving["model_options"]

    n.config = config

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

    status, condition = optimize_mga_fixed_bound(
        n,
        obj_base * (1 + slack),
        weights=weights,
        sense=sense,
        model_kwargs=model_kwargs,
        **kwargs,
    )

    if status == "ok":
        print("Solved successfully")
        return n
    if "infeasible" in condition:
        if cf_solving.get("print_infeasibilities", True):
            labels = n.model.compute_infeasibilities()
            logger.info(f"Labels:\n{labels}")
            n.model.print_infeasibilities()

    raise RuntimeError(
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
    set_scenario_config(snakemake)
    update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    n = pypsa.Network(snakemake.input.network)

    n = prepare_network(
        n,
        solve_opts,
        clusters=snakemake.wildcards.clusters,
        config=snakemake.config,
        sector=snakemake.params.sector,
        foresight="myopic",
        planning_horizons=snakemake.params.planning_horizons,
        current_horizon=snakemake.wildcards.planning_horizons,
    )

    # Base objective slack on optimal network.
    n_opt = pypsa.Network(snakemake.input.network_opt)
    obj_base = n_opt.statistics.capex().sum() + n_opt.statistics.opex().sum()
    del n_opt

    # NB: We _don't_ support custom extra functionality in this script!
    # n.custom_extra_functionality = snakemake.params.custom_extra_functionality

    with memory_logger(
        filename=getattr(snakemake.log, "memory", None), interval=30.0
    ) as mem:
        n = near_opt(
            n,
            snakemake.config,
            snakemake.params.solving,
            snakemake.params.near_opt,
            snakemake.wildcards.sense,
            obj_base,
            float(snakemake.wildcards.slack),
        )

    logger.info(f"Maximum memory usage: {mem.mem_usage}")

    n.export_to_netcdf(snakemake.output[0])
