# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
Solves linear optimal power flow for a network iteratively while updating reactances.

Relevant Settings
-----------------

.. code:: yaml

    solving:
        tmpdir:
        options:
            formulation:
            clip_p_max_pu:
            load_shedding:
            noisy_costs:
            nhours:
            min_iterations:
            max_iterations:
            skip_iterations:
            track_iterations:
        solver:
            name:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`electricity_cf`, :ref:`solving_cf`, :ref:`plotting_cf`

Inputs
------

- ``networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: confer :ref:`prepare`

Outputs
-------

- ``results/networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: Solved PyPSA network including optimisation results

    .. image:: ../img/results.png
        :scale: 40 %

Description
-----------

Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.
The optimization is based on the ``pyomo=False`` setting in the :func:`network.lopf` and  :func:`pypsa.linopf.ilopf` function.
Additionally, some extra constraints specified in :mod:`prepare_network` are added.

Solving the network in multiple iterations is motivated through the dependence of transmission line capacities and impedances.
As lines are expanded their electrical parameters change, which renders the optimisation bilinear even if the power flow
equations are linearized.
To retain the computational advantage of continuous linear programming, a sequential linear programming technique
is used, where in between iterations the line impedances are updated.
Details (and errors made through this heuristic) are discussed in the paper

- Fabian Neumann and Tom Brown. `Heuristics for Transmission Expansion Planning in Low-Carbon Energy System Models <https://arxiv.org/abs/1907.10548>`_), *16th International Conference on the European Energy Market*, 2019. `arXiv:1907.10548 <https://arxiv.org/abs/1907.10548>`_.

.. warning::
    Capital costs of existing network components are not included in the objective function,
    since for the optimisation problem they are just a constant term (no influence on optimal result).

    Therefore, these capital costs are not included in ``network.objective``!

    If you want to calculate the full total annual system costs add these to the objective value.

.. tip::
    The rule :mod:`solve_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`solve_network`.

"""

import logging
from _helpers import configure_logging

import numpy as np
import pandas as pd
import re

import pypsa
from pypsa.linopf import (get_var, define_constraints, linexpr, join_exprs,
                          network_lopf, ilopf)

from pathlib import Path
from vresutils.benchmark import memory_logger

logger = logging.getLogger(__name__)


def prepare_network(n, solve_opts):

    if 'clip_p_max_pu' in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df>solve_opts['clip_p_max_pu'], other=0., inplace=True)

    if solve_opts.get('load_shedding'):
        n.add("Carrier", "Load")
        buses_i = n.buses.query("carrier == 'AC'").index
        n.madd("Generator", buses_i, " load",
               bus=buses_i,
               carrier='load',
               sign=1e-3, # Adjust sign to measure p and p_nom in kW instead of MW
               marginal_cost=1e2, # Eur/kWh
               # intersect between macroeconomic and surveybased
               # willingness to pay
               # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
               p_nom=1e9 # kW
               )

    if solve_opts.get('noisy_costs'):
        for t in n.iterate_components(n.one_port_components):
            #if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if 'marginal_cost' in t.df:
                t.df['marginal_cost'] += (1e-2 + 2e-3 *
                                          (np.random.random(len(t.df)) - 0.5))

        for t in n.iterate_components(['Line', 'Link']):
            t.df['capital_cost'] += (1e-1 +
                2e-2*(np.random.random(len(t.df)) - 0.5)) * t.df['length']

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760. / nhours

    return n


def add_CCL_constraints(n, config):
    agg_p_nom_limits = config['electricity'].get('agg_p_nom_limits')

    try:
        agg_p_nom_minmax = pd.read_csv(agg_p_nom_limits,
                                       index_col=list(range(2)))
    except IOError:
        logger.exception("Need to specify the path to a .csv file containing "
                          "aggregate capacity limits per country in "
                          "config['electricity']['agg_p_nom_limit'].")
    logger.info("Adding per carrier generation capacity constraints for "
                "individual countries")

    gen_country = n.generators.bus.map(n.buses.country)
    # cc means country and carrier
    p_nom_per_cc = (pd.DataFrame(
                    {'p_nom': linexpr((1, get_var(n, 'Generator', 'p_nom'))),
                    'country': gen_country, 'carrier': n.generators.carrier})
                    .dropna(subset=['p_nom'])
                    .groupby(['country', 'carrier']).p_nom
                    .apply(join_exprs))
    minimum = agg_p_nom_minmax['min'].dropna()
    if not minimum.empty:
        minconstraint = define_constraints(n, p_nom_per_cc[minimum.index],
                                           '>=', minimum, 'agg_p_nom', 'min')
    maximum = agg_p_nom_minmax['max'].dropna()
    if not maximum.empty:
        maxconstraint = define_constraints(n, p_nom_per_cc[maximum.index],
                                           '<=', maximum, 'agg_p_nom', 'max')


def add_EQ_constraints(n, o, scaling=1e-1):
    float_regex = "[0-9]*\.?[0-9]+"
    level = float(re.findall(float_regex, o)[0])
    if o[-1] == 'c':
        ggrouper = n.generators.bus.map(n.buses.country)
        lgrouper = n.loads.bus.map(n.buses.country)
        sgrouper = n.storage_units.bus.map(n.buses.country)
    else:
        ggrouper = n.generators.bus
        lgrouper = n.loads.bus
        sgrouper = n.storage_units.bus
    load = n.snapshot_weightings.generators @ \
           n.loads_t.p_set.groupby(lgrouper, axis=1).sum()
    inflow = n.snapshot_weightings.stores @ \
             n.storage_units_t.inflow.groupby(sgrouper, axis=1).sum()
    inflow = inflow.reindex(load.index).fillna(0.)
    rhs = scaling * ( level * load - inflow )
    lhs_gen = linexpr((n.snapshot_weightings.generators * scaling,
                       get_var(n, "Generator", "p").T)
              ).T.groupby(ggrouper, axis=1).apply(join_exprs)
    lhs_spill = linexpr((-n.snapshot_weightings.stores * scaling,
                         get_var(n, "StorageUnit", "spill").T)
                ).T.groupby(sgrouper, axis=1).apply(join_exprs)
    lhs_spill = lhs_spill.reindex(lhs_gen.index).fillna("")
    lhs = lhs_gen + lhs_spill
    define_constraints(n, lhs, ">=", rhs, "equity", "min")


def add_BAU_constraints(n, config):
    mincaps = pd.Series(config['electricity']['BAU_mincapacities'])
    lhs = (linexpr((1, get_var(n, 'Generator', 'p_nom')))
           .groupby(n.generators.carrier).apply(join_exprs))
    define_constraints(n, lhs, '>=', mincaps[lhs.index], 'Carrier', 'bau_mincaps')


def add_SAFE_constraints(n, config):
    peakdemand = (1. + config['electricity']['SAFE_reservemargin']) *\
                  n.loads_t.p_set.sum(axis=1).max()
    conv_techs = config['plotting']['conv_techs']
    exist_conv_caps = n.generators.query('~p_nom_extendable & carrier in @conv_techs')\
                       .p_nom.sum()
    ext_gens_i = n.generators.query('carrier in @conv_techs & p_nom_extendable').index
    lhs = linexpr((1, get_var(n, 'Generator', 'p_nom')[ext_gens_i])).sum()
    rhs = peakdemand - exist_conv_caps
    define_constraints(n, lhs, '>=', rhs, 'Safe', 'mintotalcap')


def add_time_agg_constraints(n, h, time_agg_config):
    """
    Add constraints to reduce the time resolution of some technologies to `h`
    hours (corresponding to the `nHSUB` wildcard).

    The impacted technologies are given in the `time_agg_config` dictionary.

    Parameters
    ----------
    n : pypsa.Network
    h : int
        The new time resolution for the given technologies.
    time_agg_config : dict
        A configuration dictionary giving the technologies whose operations to
        aggregate.
    """

    techs = time_agg_config['aggregate_uniform']

    logger.info(f"Adding constraints to reduce time resolution of {techs} to\
                {h} hours.")

    h_hour_bins = uniform_snapshot_bins(n, h)
    aggregate_links_p(n, h_hour_bins, techs)


def add_day_night_constraints(n, time_agg_config):
    """
    Add constraints to aggregate the operations of certain
    technologies during the day and night, as specified in the
    `time_aggregation` section of the config file.

    Parameters
    ----------
    n : pypsa.Network
    time_agg_config : dict
    """
    times = time_agg_config['time_boundaries']

    # Aggregate day time operations.
    day_techs = time_agg_config['aggregate_day']
    if day_techs:
        start, end = times['day_start'], times['day_end']
        day_bins = day_snapshot_bins(n, start, end)
        aggregate_links_p(n, day_bins, day_techs)

    # Aggregate night time operations.
    night_techs = time_agg_config['aggregate_night']
    if night_techs:
        start, end = times['night_start'], times['night_end']
        night_bins = night_snapshot_bins(n, start, end)
        aggregate_links_p(n, night_bins, night_techs)


def uniform_snapshot_bins(n, h):
    """Create uniformly sized intervals of snapshots of a network.
    
    Divide all snapshots of `n` into `h`-hour intervals, and return an
    iterable over those intervals.

    Parameters
    ----------
    n : pypsa.Network
    h : int

    Returns
    -------
    iterable of pd.Series
        An iterable of subsets of `n.snapshots`.
    """
    intervals = pd.interval_range(n.snapshots[0] - pd.Timedelta("1 days"),
                                  n.snapshots[-1] + pd.Timedelta("1 days"),
                                  freq=f'{h}H',
                                  closed='left')
    # TODO: in the configuration, there is an option about closed-ness
    # of intervals. Should we consider it?

    # Collect the snapshots into bins corresponding to the storage intervals.
    bins = {i: [] for i in intervals}
    interval_iter = iter(intervals)
    current_interval = next(interval_iter)
    for s in n.snapshots:
        while s not in current_interval:
            current_interval = next(interval_iter)
        bins[current_interval].append(s)

    return [i for i in bins.values() if i]


def night_snapshot_bins(n, night_start='18:00:00', night_end='06:00:00'):
    """
    Creates an interval for each night (from `night_start` to `night_end`), and
    returns a dictionary mapping each such interval to the snapshots of `n`
    contained in it.

    Parameters
    ----------
    n : pypsa.Network
    night_start : str
        The time of day in hh:mm:ss format at which to start night intervals.
    night_end : str
        The time of day in hh:mm:ss format at which to end night intervals.
    """
    days = pd.date_range(start=n.snapshots[0].round('D'),
                         end=n.snapshots[-1].round('D'),
                         freq='D')

    def night_inter(day):
        return pd.Interval(left=day + pd.Timedelta(night_start)
                                - pd.Timedelta('1D'),
                           right=day + pd.Timedelta(night_end),
                           closed='both')
    night_intervals = pd.IntervalIndex(list(map(night_inter, days)))

    # Collect the snapshots into bins corresponding to the intervals.
    bins = {i: [] for i in night_intervals}
    for s in n.snapshots:
        # Is s during the night?
        if night_intervals.contains(s).max():
            # In which night is s?
            night_i = night_intervals.get_loc(s)
            # Add s to the right bin
            bins[night_intervals.values[night_i]].append(s)

    return bins.values()


def day_snapshot_bins(n, day_start='10:00:00', day_end='15:00:00'):
    """
    Create an interval for each day (between `day_start` and
    `day_end`), and returns a dictionary mapping each such interval to
    the snapshots of `n` contained in it.

    Parameters
    ----------
    n : pypsa.Network
    day_start : str
        The time of day in hh:mm:ss format at which to start day intervals.
    day_end : str
        The time of day in hh:mm:ss format at which to end day intervals.
    """
    days = pd.date_range(start=n.snapshots[0].round('D'),
                         end=n.snapshots[-1].round('D'),
                         freq='D')

    def day_inter(day):
        return pd.Interval(left=day + pd.Timedelta(day_start),
                           right=day + pd.Timedelta(day_end),
                           closed='both')
    day_intervals = pd.IntervalIndex(list(map(day_inter, days)))

    # Collect the snapshots into bins corresponding to the intervals.
    bins = {i: [] for i in day_intervals}
    for s in n.snapshots:
        # Is s during the day?
        if day_intervals.contains(s).max():
            # In which day is s?
            day_i = day_intervals.get_loc(s)
            # Add s to the right bin
            bins[day_intervals.values[day_i]].append(s)

    return bins.values()


def aggregate_links_p(n, snapshot_collection, carriers):
    """Aggregate operations of links in a PyPSA network.

    For each link in `n` with carrier in `carriers`, aggregate the
    power p along the link over the given `snapshot_collection`. In
    particular, `snapshot_collection` is a list of subsets of
    `n.snapshots`, and p is aggregated over each subset in
    `snapshot_collection` individually.
    
    Specifically, we add
        p_{i+1} - p_{i} = 0
    as a constraint for all i in each subsets of `snapshots`.

    Parameters
    ----------
    n : pypsa.Network
    snapshot_collection : List[pd.Series]
        a list of subsets of `n.snapshots`.
    carrier : iterable
        a subset of `n.carriers`.

    """
    # Get the index of all links with the given carrier.
    store_link_i = n.links.loc[n.links['carrier'].isin(carriers)].index

    # Add a constraint for each set of snapshots in snapshot_collection
    for snapshots in snapshot_collection:
        # Get the power time series variables for the above links. 'p' is
        # a Dataframe indexed over snapshots, with the relevant links as
        # columns.
        p = get_var(n, 'Link', 'p').loc[snapshots, store_link_i]
        # Create a table of 1-term linear expressions of all except the
        # first row of 'p'.
        lhs, *axes = linexpr((1, p.loc[snapshots[1:]]), return_axes=True)
        # Subtract from each entry in 'lhs' the variable "above" (one
        # snapshot before) it. This is achieved by 'shift()'ing the values
        # in the dataframe 'p'.
        lhs += linexpr((-1, p.shift().loc[snapshots[1:]]))\
            .reindex(index=axes[0], columns=axes[1]).values
        define_constraints(n, lhs, "=", 0, 'Link', 'aggregation')


def add_battery_constraints(n):
    nodes = n.buses.index[n.buses.carrier == "battery"]
    if nodes.empty or ('Link', 'p_nom') not in n.variables.index:
        return
    link_p_nom = get_var(n, "Link", "p_nom")
    lhs = linexpr((1,link_p_nom[nodes + " charger"]),
                  (-n.links.loc[nodes + " discharger", "efficiency"].values,
                   link_p_nom[nodes + " discharger"].values))
    define_constraints(n, lhs, "=", 0, 'Link', 'charger_ratio')


def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to ``pypsa.linopf.network_lopf``.
    If you want to enforce additional custom constraints, this is a good location to add them.
    The arguments ``opts`` and ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config
    if 'BAU' in opts and n.generators.p_nom_extendable.any():
        add_BAU_constraints(n, config)
    if 'SAFE' in opts and n.generators.p_nom_extendable.any():
        add_SAFE_constraints(n, config)
    if 'CCL' in opts and n.generators.p_nom_extendable.any():
        add_CCL_constraints(n, config)
    for o in opts:
        if "EQ" in o:
            add_EQ_constraints(n, o)
    for o in opts:
        m = re.match(r'^(\d)+hsub$', o, re.IGNORECASE)
        if m is not None:
            h = int(m.group(1))
            add_time_agg_constraints(n, h, config['time_aggregation'])
            break
    add_battery_constraints(n)
    add_day_night_constraints(n, config['time_aggregation'])


def solve_network(n, config, opts='', **kwargs):
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
                     extra_functionality=extra_functionality, **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options,
              track_iterations=track_iterations,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              extra_functionality=extra_functionality, **kwargs)
    return n


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_network', network='elec', simpl='',
                                  clusters='5', ll='copt', opts='Co2L-BAU-CCL-24H')
    configure_logging(snakemake)

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.opts.split('-')
    solve_opts = snakemake.config['solving']['options']

    fn = getattr(snakemake.log, 'memory', None)
    with memory_logger(filename=fn, interval=30.) as mem:
        n = pypsa.Network(snakemake.input[0])
        n = prepare_network(n, solve_opts)
        n = solve_network(n, config=snakemake.config, opts=opts,
                          solver_dir=tmpdir,
                          solver_logfile=snakemake.log.solver)
        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
