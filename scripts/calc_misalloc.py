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

from pathlib import Path
from vresutils.benchmark import memory_logger

import solve_network

logger = logging.getLogger(__name__)


def add_misalloc_bounds(n, snapshots):
    """Add lower bounds on capacities to calculate misallocation metric.

    The given network `n` is assumed to have the second (solved)
    network `n.nom_min_network` attached. This network must have
    busses, generators, stores, lines and links with the same names as
    `n`. For all of these components, the lower bound on their
    installed capacity in `n` is set to their optimum capacity in
    `n.nom_min_network`.

    For example, for a generator `g` we add a constraint of the form
        n.g.p_nom_min = n.nom_min_network.g.p_nom_opt.
    """
    m = n.nom_min_network

    # A helper function to factor out the common code required for the
    # different types of components.
    def add_misalloc_bound_comp(type: str, var: str) -> None:
        """Add misallocation bounds for one kind of component.

        The `type` can be one of 'Generator', 'Store', 'StorageUnit',
        'Line' or 'Link', and the `var` must be a corresponding
        variable such as 'p_nom' or 'e_nom'. The strings `'_min'` and
        `'_opt'` are appended to `var`; the results must be variables
        of the given components.

        """
        n_comps = get_var(n, type, var + '_min')
        m_comps = get_var(m, type, var + '_opt')
        lhs = linexpr((1, n_comps), (-1, m_comps))
        define_constraints(n, lhs, '>=', 0, type, 'misallocation_bound')

    # Define misallocation bounds for generators, stores, lines and
    # links.
    logger.info("Adding lower bounds for misallocation metric.")
    add_misalloc_bound_comp('Generator', 'p_nom')
    add_misalloc_bound_comp('Store', 'e_nom')
    add_misalloc_bound_comp('StorageUnit', 'p_nom')
    add_misalloc_bound_comp('Line', 's_nom')
    add_misalloc_bound_comp('Link', 'p_nom')


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
                     extra_functionality=add_misalloc_bounds, **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options,
              track_iterations=track_iterations,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              extra_functionality=add_misalloc_bounds, **kwargs)
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
    # opts = snakemake.wildcards.opts.split('-')
    # solve_opts = snakemake.config['solving']['options']

    fn = getattr(snakemake.log, 'memory', None)
    with memory_logger(filename=fn, interval=30.) as mem:

        # Check that we received exactly two networks as input to
        # compare.
        if len(snakemake.input) != 2:
            raise ValueError("calc_misalloc didn't get exactly two networks to "
                             "compare.")

        # Load the two networks.
        ns = [pypsa.Network(fn) for fn in snakemake.input]

        # First extract the objection function value from the input networks.
        obj_vals = [n.objective for n in ns]

        # To each network, attach the other network so we can use it's
        # optimal capacities as new lower bounds.
        ns[0].nom_min_network = ns[1]
        ns[1].nom_min_network = ns[0]

        # Now solve each network again, but this time add lower bounds
        # on capacities. Note that since these networks have already
        # be solved previously, we don't need to prepare them as in
        # the `solve_network` rule. We also don't need to pass the
        # `solve_network.extra_functionality` function, since it has
        # likewise already been applied.
        #
        # We use the `extra_functionality` function to add the lower
        # bounds.
        for i in range(2):
            ns[i] = solve_network.solve_network(
                ns[i], config=snakemake.config, # opts=opts,
                solver_dir=tmpdir,
                solver_logfile=snakemake.log.solver)

        # Finally calculate the misallocation metric using the new
        # objection values.
        new_obj_vals = [n.objective for n in ns]
        M = new_obj_vals[0] - obj_vals[0] + new_obj_vals[1] - obj_vals[1]
        with open(snakemake.output[0], 'w') as f:
            f.write(f"{M}")

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
