# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from os.path import normpath, exists
from shutil import copyfile, move

from scripts._helpers import parse_year_wildcard

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "config.yaml"

COSTS="data/costs.csv"
ATLITE_NPROCESSES = config['atlite'].get('nprocesses', 4)


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    # The {year} wildcard represents a set of years and consists of a
    # number of single years or ranges (of the form 2000-2020) all
    # separated by `+`s.
    year="([0-9]+)(-[0-9]+)?(\+([0-9]+)(-[0-9]+)?)*",


rule cluster_all_networks:
    input: expand("networks/elec_{year}_s{simpl}_{clusters}.nc", **config['scenario'])


rule extra_components_all_networks:
    input: expand("networks/elec_{year}_s{simpl}_{clusters}_ec.nc", **config['scenario'])


rule prepare_all_networks:
    input: expand("networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc", **config['scenario'])


rule solve_all_networks:
    input: expand("results/networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc", **config['scenario'])


if config['enable'].get('prepare_links_p_nom', False):
    rule prepare_links_p_nom:
        output: 'data/links_p_nom.csv'
        log: 'logs/prepare_links_p_nom.log'
        threads: 1
        resources: mem_mb=500
        script: 'scripts/prepare_links_p_nom.py'


datafiles = ['ch_cantons.csv', 'je-e-21.03.02.xls', 'eez/World_EEZ_v8_2014.shp',
            'hydro_capacities.csv', 'naturalearth/ne_10m_admin_0_countries.shp', 
            'NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp', 'nama_10r_3popgdp.tsv.gz', 
            'nama_10r_3gdp.tsv.gz', 'corine/g250_clc06_V18_5.tif']


if not config.get('tutorial', False):
    datafiles.extend(["natura/Natura2000_end2015.shp", "GEBCO_2014_2D.nc"])


if config['enable'].get('retrieve_databundle', True):
    rule retrieve_databundle:
        output: expand('data/bundle/{file}', file=datafiles)
        log: "logs/retrieve_databundle.log"
        conda: "envs/environment.yaml"
        script: 'scripts/retrieve_databundle.py'
    

rule build_powerplants:
    input:
        base_network="networks/base.nc",
        custom_powerplants="data/custom_powerplants.csv"
    output: "resources/powerplants.csv"
    log: "logs/build_powerplants.log"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=500
    script: "scripts/build_powerplants.py"


rule base_network:
    input:
        eg_buses='data/entsoegridkit/buses.csv',
        eg_lines='data/entsoegridkit/lines.csv',
        eg_links='data/entsoegridkit/links.csv',
        eg_converters='data/entsoegridkit/converters.csv',
        eg_transformers='data/entsoegridkit/transformers.csv',
        parameter_corrections='data/parameter_corrections.yaml',
        links_p_nom='data/links_p_nom.csv',
        links_tyndp='data/links_tyndp.csv',
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        europe_shape='resources/europe_shape.geojson'
    output: "networks/base.nc"
    log: "logs/base_network.log"
    conda: "envs/environment.yaml"
    benchmark: "benchmarks/base_network"
    threads: 1
    resources: mem_mb=500
    script: "scripts/base_network.py"


rule build_shapes:
    input:
        naturalearth='data/bundle/naturalearth/ne_10m_admin_0_countries.shp',
        eez='data/bundle/eez/World_EEZ_v8_2014.shp',
        nuts3='data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp',
        nuts3pop='data/bundle/nama_10r_3popgdp.tsv.gz',
        nuts3gdp='data/bundle/nama_10r_3gdp.tsv.gz',
        ch_cantons='data/bundle/ch_cantons.csv',
        ch_popgdp='data/bundle/je-e-21.03.02.xls'
    output:
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        europe_shape='resources/europe_shape.geojson',
        nuts3_shapes='resources/nuts3_shapes.geojson'
    log: "logs/build_shapes.log"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=500
    script: "scripts/build_shapes.py"


rule build_bus_regions:
    input:
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        base_network="networks/base.nc"
    output:
        regions_onshore="resources/regions_onshore.geojson",
        regions_offshore="resources/regions_offshore.geojson"
    log: "logs/build_bus_regions.log"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=1000
    script: "scripts/build_bus_regions.py"

if config['enable'].get('build_cutout', False):
    rule build_cutout:
        input: 
            regions_onshore="resources/regions_onshore.geojson",
            regions_offshore="resources/regions_offshore.geojson"
        output: "cutouts/{cutout}.nc"
        log: "logs/build_cutout/{cutout}.log"
        benchmark: "benchmarks/build_cutout_{cutout}"
        conda: "envs/environment.yaml"
        threads: ATLITE_NPROCESSES
        resources: mem_mb=ATLITE_NPROCESSES * 1000
        script: "scripts/build_cutout.py"


# For now, the user has to download the required cutouts manually and
# place them in the "cutouts" directory.

# if config['enable'].get('retrieve_cutout', True):
#     rule retrieve_cutout:
#         input: HTTP.remote("zenodo.org/record/6382570/files/{cutout}.nc", keep_local=True, static=True)
#         output: "cutouts/{cutout}.nc"
#         shell: "mv {input} {output}"

if config['enable'].get('build_natura_raster', False):
    rule build_natura_raster:
        input:
            natura="data/bundle/natura/Natura2000_end2015.shp",
            cutouts=expand("cutouts/{cutouts}.nc", **config['atlite'])
        output: "resources/natura.tiff"
        log: "logs/build_natura_raster.log"
        conda: "envs/environment.yaml"
        script: "scripts/build_natura_raster.py"


if config['enable'].get('retrieve_natura_raster', True):
    rule retrieve_natura_raster:
        input: HTTP.remote("zenodo.org/record/4706686/files/natura.tiff", keep_local=True, static=True)
        output: "resources/natura.tiff"
        run: move(input[0], output[0])


def renewable_profiles_cutouts(wildcards):
    # Parse the set of years over which the required network is defined.
    years = parse_year_wildcard(wildcards.year)
    # In case the year boundary is not the 1st of January, we also
    # need the cutout for year y+1 for each y in `years`.
    if config['snapshots'].get('year_boundary', '01-01') != '01-01':
        years = set(years + [y + 1 for y in years])
    return [
        f"cutouts/{config['renewable'][wildcards.technology]['cutout']}_{y}.nc"
        for y in years
    ]


def build_renewable_profiles_memory(wildcards):
    """Estimate memory requirements of `build_renewable_profiles`"""
    num_years = len(parse_year_wildcard(wildcards.year))
    return ATLITE_NPROCESSES * (4750 + 250 * num_years)


rule build_renewable_profiles:
    input:
        corine="data/bundle/corine/g250_clc06_V18_5.tif",
        natura="resources/natura.tiff",
        gebco=lambda w: ("data/bundle/GEBCO_2014_2D.nc"
                         if "max_depth" in config["renewable"][w.technology].keys()
                         else []),
        country_shapes='resources/country_shapes.geojson',
        offshore_shapes='resources/offshore_shapes.geojson',
        regions=lambda w: ("resources/regions_onshore.geojson"
                                   if w.technology in ('onwind', 'solar')
                                   else "resources/regions_offshore.geojson"),
        cutouts=renewable_profiles_cutouts
    output: profile="resources/profile_{technology}_{year}.nc",
    log: "logs/build_renewable_profile_{technology}_{year}.log"
    benchmark: "benchmarks/build_renewable_profiles_{technology}_{year}"
    conda: "envs/environment.yaml"
    threads: ATLITE_NPROCESSES
    resources: mem_mb=build_renewable_profiles_memory
    wildcard_constraints: technology="(?!hydro).*" # Any technology other than hydro
    script: "scripts/build_renewable_profiles.py"



ruleorder: build_hydro_profile > build_renewable_profiles

def hydro_profiles_cutouts(wildcards):
    if "hydro" in config["renewable"]:
        years = parse_year_wildcard(wildcards.year)
        # In case the year boundary is not the 1st of January, we also
        # need the cutout for year y+1 for each y in `years`.
        if config['snapshots'].get('year_boundary', '01-01') != '01-01':
            years = set(years + [y + 1 for y in years])
        return [
            f"cutouts/{config['renewable']['hydro']['cutout']}_{y}.nc"
            for y in years
        ]
    else:
        return "config['renewable']['hydro']['cutout'] not configured"

rule build_hydro_profile:
    input:
        country_shapes='resources/country_shapes.geojson',
        eia_hydro_generation='data/EIA_scaled_hydro_1980_2020.csv',
        cutouts=hydro_profiles_cutouts,
    output: 'resources/profile_hydro_{year}.nc'
    log: "logs/build_hydro_profile_{year}.log"
    benchmark: "benchmarks/build_renewable_profiles_hydro_{year}"
    conda: "envs/environment.yaml"
    resources: mem_mb=build_renewable_profiles_memory
    script: 'scripts/build_hydro_profile.py'


def add_electricity_memory(wildcards):
    """Estimate memory requirements of `add_electricity`"""
    num_years = len(parse_year_wildcard(wildcards.year))
    return 3000 + num_years * 2000


rule add_electricity:
    input:
        base_network='networks/base.nc',
        tech_costs=COSTS,
        regions="resources/regions_onshore.geojson",
        powerplants='resources/powerplants.csv',
        hydro_capacities='data/bundle/hydro_capacities.csv',
        geth_hydro_capacities='data/geth2015_hydro_capacities.csv',
        load='data/europe_demand_artificial_1980-2020.csv',
        nuts3_shapes='resources/nuts3_shapes.geojson',
        **{f"profile_{tech}": "resources/profile_" + str(tech) + "_{year}.nc"
           for tech in config['renewable']}
    output: "networks/elec_{year}.nc"
    log: "logs/add_electricity_{year}.log"
    benchmark: "benchmarks/add_electricity_{year}"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=add_electricity_memory
    script: "scripts/add_electricity.py"


def simplify_memory(wildcards):
    """Estimate memory requirements of `simplify_network`"""
    num_years = len(parse_year_wildcard(wildcards.year))
    return 2000 + num_years * 2000


rule simplify_network:
    input:
        network='networks/elec_{year}.nc',
        network_constant=f"networks/elec_{config['net_clustering_year']}.nc",
        tech_costs=COSTS,
        regions_onshore="resources/regions_onshore.geojson",
        regions_offshore="resources/regions_offshore.geojson"
    output:
        network='networks/elec_{year}_s{simpl}.nc',
        # Note that the following output files should actually for
        # equal for all years, but we can't remove the {year} wildcard
        # from the names since that would lead to output filename
        # conflicts between different years.
        regions_onshore="resources/regions_onshore_elec_{year}_s{simpl}.geojson",
        regions_offshore="resources/regions_offshore_elec_{year}_s{simpl}.geojson",
        busmap='resources/busmap_elec_{year}_s{simpl}.csv',
        connection_costs='resources/connection_costs_{year}_s{simpl}.csv'
    log: "logs/simplify_network/elec_{year}_s{simpl}.log"
    benchmark: "benchmarks/simplify_network/elec_{year}_s{simpl}"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=simplify_memory
    script: "scripts/simplify_network.py"


rule cluster_network:
    input:
        network='networks/elec_{year}_s{simpl}.nc',
        network_constant='networks/elec_' + str(config['net_clustering_year']) + '_s{simpl}.nc',
        regions_onshore="resources/regions_onshore_elec_{year}_s{simpl}.geojson",
        regions_offshore="resources/regions_offshore_elec_{year}_s{simpl}.geojson",
        busmap=ancient('resources/busmap_elec_{year}_s{simpl}.csv'),
        custom_busmap=("data/custom_busmap_elec_s{simpl}_{clusters}.csv"
                       if config["enable"].get("custom_busmap", False) else []),
        tech_costs=COSTS
    output:
        network='networks/elec_{year}_s{simpl}_{clusters}.nc',
        regions_onshore="resources/regions_onshore_elec_{year}_s{simpl}_{clusters}.geojson",
        regions_offshore="resources/regions_offshore_elec_{year}_s{simpl}_{clusters}.geojson",
        busmap="resources/busmap_elec_{year}_s{simpl}_{clusters}.csv",
        linemap="resources/linemap_elec_{year}_s{simpl}_{clusters}.csv"
    log: "logs/cluster_network/elec_{year}_s{simpl}_{clusters}.log"
    benchmark: "benchmarks/cluster_network/elec_{year}_s{simpl}_{clusters}"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=simplify_memory
    script: "scripts/cluster_network.py"


rule add_extra_components:
    input:
        network='networks/elec_{year}_s{simpl}_{clusters}.nc',
        tech_costs=COSTS,
    output: 'networks/elec_{year}_s{simpl}_{clusters}_ec.nc'
    log: "logs/add_extra_components/elec_{year}_s{simpl}_{clusters}.log"
    benchmark: "benchmarks/add_extra_components/elec_{year}_s{simpl}_{clusters}_ec"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=3000
    script: "scripts/add_extra_components.py"


rule prepare_network:
    input: 'networks/elec_{year}_s{simpl}_{clusters}_ec.nc', tech_costs=COSTS
    output: 'networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc'
    log: "logs/prepare_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.log"
    benchmark: "benchmarks/prepare_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    conda: "envs/environment.yaml"
    threads: 1
    resources: mem_mb=4000
    script: "scripts/prepare_network.py"


def memory(w):
    factor = 3.
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)h$', o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split('-'):
        m = re.match(r'^(\d+)seg$', o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith('m'):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


rule solve_network:
    input: "networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
    output: "results/networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
    log:
        solver=normpath("logs/solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"),
        python="logs/solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
        memory="logs/solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_memory.log"
    benchmark: "benchmarks/solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    conda: "envs/environment.yaml"
    threads: 4
    resources: mem_mb=memory
    shadow: "minimal"
    script: "scripts/solve_network.py"


rule solve_operations_network:
    input:
        unprepared="networks/elec_{year}_s{simpl}_{clusters}_ec.nc",
        optimized="results/networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc"
    output: "results/networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op.nc"
    log:
        solver=normpath("logs/solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_solver.log"),
        python="logs/solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_python.log",
        memory="logs/solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_memory.log"
    benchmark: "benchmarks/solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}"
    conda: "envs/environment.yaml"
    threads: 4
    resources: mem_mb=(lambda w: 5000 + 372 * int(w.clusters))
    shadow: "minimal"
    script: "scripts/solve_operations_network.py"


rule plot_network:
    input:
        network="results/networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        tech_costs=COSTS
    output:
        only_map="results/plots/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}.{ext}",
        ext="results/plots/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_ext.{ext}"
    log: "logs/plot_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_{ext}.log"
    conda: "envs/environment.yaml"
    script: "scripts/plot_network.py"


def input_make_summary(w):
    # It's mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return ([COSTS] +
            expand("results/networks/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
                   ll=ll,
                   **{k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
                      for k in ["simpl", "clusters", "opts"]}))


rule make_summary:
    input: input_make_summary
    output: directory("results/summaries/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}")
    log: "logs/make_summary/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.log",
    conda: "envs/environment.yaml"
    script: "scripts/make_summary.py"


rule plot_summary:
    input: "results/summaries/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}"
    output: "results/plots/summary_{summary}_elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}"
    log: "logs/plot_summary/{summary}_elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{ext}.log"
    conda: "envs/environment.yaml"
    script: "scripts/plot_summary.py"


def input_plot_p_nom_max(w):
    return [("networks/elec_{year}_s{simpl}{maybe_cluster}.nc"
             .format(maybe_cluster=('' if c == 'full' else ('_' + c)), **w))
            for c in w.clusts.split(",")]


rule plot_p_nom_max:
    input: input_plot_p_nom_max
    output: "results/plots/elec_{year}_s{simpl}_cum_p_nom_max_{clusts}_{techs}_{country}.{ext}"
    log: "logs/plot_p_nom_max/elec_{year}_s{simpl}_{clusts}_{techs}_{country}_{ext}.log"
    conda: "envs/environment.yaml"
    script: "scripts/plot_p_nom_max.py"

