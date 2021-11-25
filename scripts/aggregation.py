"""Utility functions for variable aggregation in a pypsa network."""


from typing import (List, Collection)
import re

import pandas as pd

import pypsa
from pypsa.linopf import (get_var, define_constraints, linexpr)


def aggregate_links_p(n: pypsa.Network,
                      snapshot_collection: Collection[Collection[pd.Timestamp]],
                      carriers: Collection[str],
                      strategy: str) -> None:
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
    n : 
    snapshot_collection : 
        a list of subsets of `n.snapshots`.
    carrier : 
        a subset of `n.carriers`.
    strategy :
        either 'con' or 'sub', indicating whether to aggregate using
        constraints or substitution, respectively.
    """
    if not snapshot_collection:
        raise ValueError("Cannot aggregate over empty collection of snapshots.")

    # Get the index of all links with the given carrier.
    link_i = n.links.loc[n.links['carrier'].isin(carriers)].index

    # Collect all decision variables for those links.
    vars = get_var(n, 'Link', 'p').loc[:, link_i]

    if strategy == 'con':
        aggregate_by_constraints(n, vars, snapshot_collection)
    elif strategy == 'sub':
        aggregate_by_substitution(n, vars, snapshot_collection)
    else:
        raise ValueError(f"Aggregation strategy {strategy} not recognised.")


def aggregate_by_substitution(
        n: pypsa.Network,
        vars: pd.DataFrame,
        snapshot_collection: Collection[Collection[pd.Timestamp]]) -> None:
    """Aggregate the given variables by substitution."""
    # First compile a dictionary of substitutions.
    subs = {}
    for snapshots in snapshot_collection:
        vars_ = vars.loc[snapshots]
        for col in vars_.columns:
            var_array = vars_[col].values
            for var in var_array[1:]:
                subs['x' + str(var)] = 'x' + str(var_array[0])

    def sub_fun(m):
        if m.group(0) in subs:
            return subs[m.group(0)]
        else:
            return m.group(0)

    # Now go through each LP file and perform the variable substitutions.
    for f in [n.objective_f, n.constraints_f, n.bounds_f, n.binaries_f]:
        f.seek(0)
        content = f.read()
        content = re.sub(r'\bx\d+\b',
                         sub_fun, content)
        f.seek(0)
        f.truncate()
        f.write(content)


def aggregate_by_constraints(
        n: pypsa.Network,
        vars: pd.DataFrame,
        snapshot_collection: Collection[Collection[pd.Timestamp]]) -> None:
    """Aggregate the given variables by adding constraints."""
    # Add all constraints in one go by looping over the snapshot
    # collection. For each set of snapshots (over which we must
    # aggregate), subtract the last n-1 variables from the first n-1
    # variable on order to set each variable equal to the next. The
    # constraints must be added this way since only constants are
    # allowed on the right hand side.
    lhs = pd.concat([linexpr((1, vars.loc[sns[1:]].values),
                             (-1, vars.loc[sns[:-1]].values))
                     for sns in snapshot_collection],
                    axis=0)
    # Indexing the constraints properly makes it possible to unpack
    # the associated dual variables after solving.
    lhs.index = pd.Index([sns[0] for sns in snapshot_collection])
    define_constraints(n, lhs, "=", 0, 'Link', 'link aggregation')


def uniform_snapshot_bins(n: pypsa.Network, h: int) -> List[List[pd.Timestamp]]:
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
    intervals = pd.interval_range(n.snapshots[0] - pd.Timedelta(f"{h}H"),
                                  n.snapshots[-1] + pd.Timedelta(f"{h}H"),
                                  freq=f'{h}H',
                                  closed='left')

    # Collect the snapshots into bins corresponding to the storage intervals.
    bins = {i: [] for i in intervals}
    interval_iter = iter(intervals)
    current_interval = next(interval_iter)
    for s in n.snapshots:
        while s not in current_interval:
            current_interval = next(interval_iter)
        bins[current_interval].append(s)

    return [i for i in bins.values() if len(i) >= 2]


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


def add_battery_constraints(n):
    nodes = n.buses.index[n.buses.carrier == "battery"]
    if nodes.empty or ('Link', 'p_nom') not in n.variables.index:
        return
    link_p_nom = get_var(n, "Link", "p_nom")
    lhs = linexpr((1,link_p_nom[nodes + " charger"]),
                  (-n.links.loc[nodes + " discharger", "efficiency"].values,
                   link_p_nom[nodes + " discharger"].values))
    define_constraints(n, lhs, "=", 0, 'Link', 'charger_ratio')
