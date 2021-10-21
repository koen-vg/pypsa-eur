# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Calculates the misallocation metric between two scenarios.

Relevant Settings
-----------------

.. code:: yaml

    misallocation:
        scenario_a:
        scenario_b:

.. seealso::
    TODO: reference to rule `solve_network`.

Inputs
------

- Two solved networks of the form ``results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``, possibly with differing options.

Outputs
-------

- ``results/misalloc``

Description
-----------

Compute the misallocation metric between two models, following
http://arxiv.org/abs/2007.09956. The misallocation metric compares two
different models of the same underlaying system, and is a measure for
how sensitive the model solutions are to the differences in assumption
or formulation between the two models. See the linked paper for
details.

The settings of two models to be compared are given in the
`misallocation` section of `config.yaml`. Two scenarios must be given
in the same format as the `scenario` section, and must define a single
model. (That is, don't specify multiple options or numbers of
clusters.) See the default configuration for an example.

The misallocation number is a single number which is written to
`results/misalloc`.
"""

import logging
from _helpers import configure_logging

import pypsa
from pypsa.linopf import (get_var, define_constraints, linexpr,
                          network_lopf, ilopf)
import solve_network

import pandas as pd

from pathlib import Path
import os
from vresutils.benchmark import memory_logger

logger = logging.getLogger(__name__)


def add_misalloc_bounds(n):
    """Add lower bounds on capacities to calculate misallocation metric.

    The given network `n` is assumed to have the second (solved)
    network `n.nom_min_network` attached. This network must have
    busses, generators, stores, lines and links with the same names as
    `n`. For all of these components, the lower bound on their
    installed capacity in `n` is set to their optimum capacity in
    `n.nom_min_network`.

    For example, for a generator `g` we add a constraint of the form
        n.g.p_nom => n.nom_min_network.g.p_nom_opt.
    """
    m = n.nom_min_network

    # A helper function to factor out the common code required for the
    # different types of components.
    def add_misalloc_bound_comp(comp: str,
                                m_comp_df: pd.DataFrame,
                                var: str) -> None:
        """Add misallocation bounds for one kind of component.

        The `comp` can be one of 'Generator', 'Store', 'StorageUnit',
        'Line' or 'Link', and the `var` must be a corresponding
        variable such as 'p_nom' or 'e_nom'. (See `n.variables` for
        possible choices.)

        `m_comp_df` should be the DataFrame of m with the data
        corresponding to `comp`. For example, when `comp` is
        'Generator', then `m_comp_df` should be `m.generators`.

        """
        comp_vars = get_var(n, comp, var)
        comp_min = m_comp_df.loc[comp_vars.index, var + '_opt']
        define_constraints(n, linexpr((1, comp_vars)), '>=', comp_min,
                           comp, 'misallocation_bound')

    # Define misallocation bounds for generators, stores, lines and
    # links.
    # TODO: find a way of getting just the extendable components.
    logger.info("Adding lower bounds for misallocation metric.")
    add_misalloc_bound_comp('Generator', m.generators, 'p_nom')
    add_misalloc_bound_comp('Store', m.stores, 'e_nom')
    # add_misalloc_bound_comp('StorageUnit', m.storage_units, 'p_nom')
    add_misalloc_bound_comp('Line', m.lines, 's_nom')
    add_misalloc_bound_comp('Link', m.links, 'p_nom')


def misalloc_extra_functionality(n, snapshots):
    """Apply misallocation bounds and other functionality."""
    solve_network.extra_functionality(n, snapshots)
    add_misalloc_bounds(n)


def misalloc_solve_network(n, config, opts='', **kwargs):
    """Taken from the `solve_network` module with slight modifications.

    The only difference is that we want to solve with a different
    `extra_functionality` function.
    """
    solver_options = config['solving']['solver'].copy()
    solver_name = solver_options.pop('name')
    cf_solving = config['solving']['options']
    track_iterations = cf_solving.get('track_iterations', False)
    min_iterations = cf_solving.get('min_iterations', 4)
    max_iterations = cf_solving.get('max_iterations', 6)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if cf_solving.get('skip_iterations', False):
        network_lopf(n, solver_name=solver_name, solver_options=solver_options,
                     extra_functionality=misalloc_extra_functionality, **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options,
              track_iterations=track_iterations,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              extra_functionality=misalloc_extra_functionality, **kwargs)
    return n


if __name__ == "__main__":
    # Boilerplate code similar to the code in the `solve_network`
    # script:
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_network', network='elec', simpl='',
                                  clusters='5', ll='copt', opts='Co2L-BAU-CCL-24H')
    configure_logging(snakemake)

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = [snakemake.wildcards.optsA.split('-'),
            snakemake.wildcards.optsB.split('-')]
    solve_opts = snakemake.config['solving']['options']

    fn = getattr(snakemake.log, 'memory', None)
    with memory_logger(filename=fn, interval=30.) as mem:
        # Load the two networks.
        ns = [pypsa.Network(snakemake.input.networkA),
              pypsa.Network(snakemake.input.networkB)]
        solved_ns = [pypsa.Network(snakemake.input.solved_networkA),
                     pypsa.Network(snakemake.input.solved_networkB)]

        ns = [solve_network.prepare_network(n, solve_opts) for n in ns]

        # First extract the objection function value from the input networks.
        obj_vals = [n.objective for n in solved_ns]

        # To each network, attach the solved version of the other
        # network so we can use it's optimal capacities as new lower
        # bounds.
        ns[0].nom_min_network = solved_ns[1]
        ns[1].nom_min_network = solved_ns[0]

        # Now solve each network again, but this time add lower bounds
        # on capacities. In `misalloc_solve_network`, we use the
        # `extra_functionality` function to add the lower bounds.
        for i in range(2):
            ns[i] = misalloc_solve_network(
                ns[i], config=snakemake.config, opts=opts[i],
                solver_dir=tmpdir,
                solver_logfile=snakemake.log.solver)
        ns[0].export_to_netcdf(snakemake.output.misallocated_networkA)
        ns[1].export_to_netcdf(snakemake.output.misallocated_networkB)

        # Finally calculate the misallocation metric using the new
        # objection values.
        new_obj_vals = [n.objective for n in ns]
        M = new_obj_vals[0] - obj_vals[0] + new_obj_vals[1] - obj_vals[1]
        logger.info(f"Obj of original networks: {obj_vals[0]} and {obj_vals[1]}")
        logger.info(f"Obj of bounded networks: {new_obj_vals[0]} and {new_obj_vals[1]}")
        logger.info(f"Differences: {new_obj_vals[0] - obj_vals[0]} and {new_obj_vals[1] - obj_vals[1]}")
        with open(snakemake.output.metric, 'w') as f:
            f.write(f"{M}")

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
