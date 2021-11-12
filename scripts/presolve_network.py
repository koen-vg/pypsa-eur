
import os
import pandas as pd
import gurobipy
from tempfile import mkstemp

import pypsa
from pypsa.pf import (_as_snapshots, get_switchable_as_dense as get_as_dense)
from pypsa.descriptors import (get_bounds_pu, get_extendable_i, get_non_extendable_i,
                          expand_series, nominal_attrs, additional_linkports,
                          Dict, get_active_assets, get_activity_mask)
from pypsa.linopt import (linexpr, write_bound, write_constraint, write_objective,
                     set_conref, set_varref, get_con, get_var, join_exprs,
                     run_and_read_cbc, run_and_read_gurobi, run_and_read_glpk,
                     run_and_read_cplex, run_and_read_xpress,
                     define_constraints, define_variables, define_binaries,
                     align_with_static_component)
from pypsa.linopf import prepare_lopf

from solve_network import prepare_network, extra_functionality

num_presolves = 5

if __name__ == "__main__":
    """Perform a number of presolves of the given network.

    Follow the steps of the `solve_network` script and
    `pypsa.network_lopf`, up to the point of presolving a network. Do
    this several times, and write some statistics about the number of
    non-zero coefficients in the resulting LP."""

    stats = pd.DataFrame()
    for i in range(num_presolves):
        tmpdir = snakemake.config['solving'].get('tmpdir')
        if tmpdir is not None:
            Path(tmpdir).mkdir(parents=True, exist_ok=True)
        opts = snakemake.wildcards.opts.split('-')
        solve_opts = snakemake.config['solving']['options']

        n = pypsa.Network(snakemake.input.network)
        n = prepare_network(n, solve_opts)

        solver_options = snakemake.config["solving"]["solver"].copy()
        solver_options.pop('name')
        # add to network for extra_functionality
        n.config = snakemake.config
        n.opts = opts

        snapshots = n.snapshots
        n.calculate_dependent_values()
        n.determine_network_topology()
        n._multi_invest = 0

        fdp, problem_fn = prepare_lopf(n, snapshots, True, False,
                                       extra_functionality, None)

        m = gurobipy.read(problem_fn)

        for key, value in solver_options.items():
            m.setParam(key, value)

        p = m.presolve()
        stats.loc[i, "presolve_numNZs"] = p.numNZs

    stats.to_csv(snakemake.output.stats)
