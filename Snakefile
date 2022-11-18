# SPDX-FileCopyrightText: : 2017-2022 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from os.path import normpath, exists
from shutil import copyfile, move

from scripts._helpers import parse_year_wildcard

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
CDIR = RDIR if not run.get("shared_cutouts") else ""

COSTS = "resources/" + RDIR + "costs.csv"
ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 4)


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
    input:
        expand(
            "networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}.nc",
            **config["scenario"]
        ),


rule extra_components_all_networks:
    input:
        expand(
            "networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}_ec.nc",
            **config["scenario"]
        ),


rule prepare_all_networks:
    input:
        expand(
            "networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule solve_all_networks:
    input:
        expand(
            "results/networks/"
            + RDIR
            + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


if config["enable"].get("prepare_links_p_nom", False):

    rule prepare_links_p_nom:
        output:
            "data/links_p_nom.csv",
        log:
            "logs/" + RDIR + "prepare_links_p_nom.log",
        threads: 1
        resources:
            mem_mb=500,
        script:
            "scripts/prepare_links_p_nom.py"


datafiles = [
    "ch_cantons.csv",
    "je-e-21.03.02.xls",
    "eez/World_EEZ_v8_2014.shp",
    "hydro_capacities.csv",
    "naturalearth/ne_10m_admin_0_countries.shp",
    "NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp",
    "nama_10r_3popgdp.tsv.gz",
    "nama_10r_3gdp.tsv.gz",
    "corine/g250_clc06_V18_5.tif",
]


if not config.get("tutorial", False):
    datafiles.extend(["natura/Natura2000_end2015.shp", "GEBCO_2014_2D.nc"])


if config["enable"].get("retrieve_databundle", True):

    rule retrieve_databundle:
        output:
            expand("data/bundle/{file}", file=datafiles),
        log:
            "logs/" + RDIR + "retrieve_databundle.log",
        resources:
            mem_mb=1000,
        conda:
            "envs/environment.fixed.yaml"
        script:
            "scripts/retrieve_databundle.py"


rule retrieve_load_data:
    input:
        HTTP.remote(
            "data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_singleindex.csv",
            keep_local=True,
            static=True,
        ),
    output:
        "data/load_raw.csv",
    resources:
        mem_mb=5000,
    run:
        move(input[0], output[0])


rule build_load_data:
    input:
        "data/load_raw.csv",
    output:
        "resources/" + RDIR + "load.csv",
    log:
        "logs/" + RDIR + "build_load_data.log",
    resources:
        mem_mb=5000,
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/build_load_data.py"


rule build_powerplants:
    input:
        base_network="networks/" + RDIR + "base.nc",
        custom_powerplants="data/custom_powerplants.csv",
    output:
        "resources/" + RDIR + "powerplants.csv",
    log:
        "logs/" + RDIR + "build_powerplants.log",
    conda:
        "envs/environment.fixed.yaml"
    threads: 1
    resources:
        mem_mb=5000,
    script:
        "scripts/build_powerplants.py"


rule base_network:
    input:
        eg_buses="data/entsoegridkit/buses.csv",
        eg_lines="data/entsoegridkit/lines.csv",
        eg_links="data/entsoegridkit/links.csv",
        eg_converters="data/entsoegridkit/converters.csv",
        eg_transformers="data/entsoegridkit/transformers.csv",
        parameter_corrections="data/parameter_corrections.yaml",
        links_p_nom="data/links_p_nom.csv",
        links_tyndp="data/links_tyndp.csv",
        country_shapes="resources/" + RDIR + "country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "offshore_shapes.geojson",
        europe_shape="resources/" + RDIR + "europe_shape.geojson",
    output:
        "networks/" + RDIR + "base.nc",
    log:
        "logs/" + RDIR + "base_network.log",
    benchmark:
        "benchmarks/" + RDIR + "base_network"
    threads: 1
    resources:
        mem_mb=500,
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/base_network.py"


rule build_shapes:
    input:
        naturalearth="data/bundle/naturalearth/ne_10m_admin_0_countries.shp",
        eez="data/bundle/eez/World_EEZ_v8_2014.shp",
        nuts3="data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp",
        nuts3pop="data/bundle/nama_10r_3popgdp.tsv.gz",
        nuts3gdp="data/bundle/nama_10r_3gdp.tsv.gz",
        ch_cantons="data/bundle/ch_cantons.csv",
        ch_popgdp="data/bundle/je-e-21.03.02.xls",
    output:
        country_shapes="resources/" + RDIR + "country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "offshore_shapes.geojson",
        europe_shape="resources/" + RDIR + "europe_shape.geojson",
        nuts3_shapes="resources/" + RDIR + "nuts3_shapes.geojson",
    log:
        "logs/" + RDIR + "build_shapes.log",
    threads: 1
    resources:
        mem_mb=500,
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/build_shapes.py"


rule build_bus_regions:
    input:
        country_shapes="resources/" + RDIR + "country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "offshore_shapes.geojson",
        base_network="networks/" + RDIR + "base.nc",
    output:
        regions_onshore="resources/" + RDIR + "regions_onshore.geojson",
        regions_offshore="resources/" + RDIR + "regions_offshore.geojson",
    log:
        "logs/" + RDIR + "build_bus_regions.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/build_bus_regions.py"


if config["enable"].get("build_cutout", False):

    rule build_cutout:
        input:
            regions_onshore="resources/" + RDIR + "regions_onshore.geojson",
            regions_offshore="resources/" + RDIR + "regions_offshore.geojson",
        output:
            "cutouts/" + CDIR + "{cutout}_{year}.nc",
        log:
            "logs/" + CDIR + "build_cutout/{cutout}_{year}.log",
        benchmark:
            "benchmarks/" + CDIR + "build_cutout_{cutout}_{year}"
        threads: ATLITE_NPROCESSES
        resources:
            mem_mb=ATLITE_NPROCESSES * 1000,
        conda:
            "envs/environment.fixed.yaml"
        script:
            "scripts/build_cutout.py"


if config["enable"].get("retrieve_cutout", True):

    rule retrieve_cutout:
        input:
            HTTP.remote(
                "https://ns9999k.webs.sigma2.no/10.11582_2022.00034/{cutout}.nc",
                keep_local=True,
                static=True,
            ),
        output:
            "cutouts/" + CDIR + "{cutout}.nc",
        log:
            "logs/" + CDIR + "retrieve_cutout_{cutout}.log",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


if config["enable"].get("retrieve_cost_data", True):

    rule retrieve_cost_data:
        input:
            HTTP.remote(
                f"raw.githubusercontent.com/PyPSA/technology-data/{config['costs']['version']}/outputs/costs_{config['costs']['year']}.csv",
                keep_local=True,
            ),
        output:
            COSTS,
        log:
            "logs/" + RDIR + "retrieve_cost_data.log",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


def raster_cutouts(_):
    """Return a list of cutout names used to compute rasters."""
    # Select the clustering year as the "canonical" year from which to
    # compute rasters.
    y = config["net_clustering_year"]
    # Collect the (year-independent) cutout names used to computed
    # renewable potential.
    names = set(
        [
            config["renewable"][tech]["cutout"]
            for tech in config["electricity"]["renewable_carriers"]
        ]
    )
    return [f"cutouts/{CDIR}{name}_{y}.nc" for name in names]


if config["enable"].get("build_natura_raster", False):

    rule build_natura_raster:
        input:
            natura="data/bundle/natura/Natura2000_end2015.shp",
            cutouts=raster_cutouts,
        output:
            "resources/" + RDIR + "natura.tiff",
        resources:
            mem_mb=5000,
        log:
            "logs/" + RDIR + "build_natura_raster.log",
        conda:
            "envs/environment.fixed.yaml"
        script:
            "scripts/build_natura_raster.py"


if config["enable"].get("retrieve_natura_raster", True):

    rule retrieve_natura_raster:
        input:
            HTTP.remote(
                "zenodo.org/record/4706686/files/natura.tiff",
                keep_local=True,
                static=True,
            ),
        output:
            "resources/" + RDIR + "natura.tiff",
        resources:
            mem_mb=5000,
        run:
            move(input[0], output[0])


rule retrieve_ship_raster:
    input:
        HTTP.remote(
            "https://zenodo.org/record/6953563/files/shipdensity_global.zip",
            keep_local=True,
            static=True,
        ),
    output:
        "data/shipdensity_global.zip",
    resources:
        mem_mb=5000,
    run:
        move(input[0], output[0])


rule build_ship_raster:
    input:
        ship_density="data/shipdensity_global.zip",
        cutouts=raster_cutouts,
    output:
        "resources/" + RDIR + "shipdensity_raster.nc",
    log:
        "logs/" + RDIR + "build_ship_raster.log",
    resources:
        mem_mb=5000,
    benchmark:
        "benchmarks/" + RDIR + "build_ship_raster"
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/build_ship_raster.py"


def renewable_profiles_cutouts(wildcards):
    # Parse the set of years over which the required network is defined.
    years = parse_year_wildcard(wildcards.year)
    # In case the year boundary is not the 1st of January, we also
    # need the cutout for year y+1 for each y in `years`.
    if config["snapshots"].get("year_boundary", "01-01") != "01-01":
        years = set(years + [y + 1 for y in years])
    return [
        "cutouts/"
        + CDIR
        + f"{config['renewable'][wildcards.technology]['cutout']}_{y}.nc"
        for y in years
    ]


def build_renewable_profiles_memory(wildcards):
    """Estimate memory requirements of `build_renewable_profiles`"""
    num_years = len(parse_year_wildcard(wildcards.year))
    return 3000 + 150 * num_years


rule build_renewable_profiles:
    input:
        base_network="networks/" + RDIR + "base.nc",
        corine="data/bundle/corine/g250_clc06_V18_5.tif",
        natura=lambda w: (
            "resources/" + RDIR + "natura.tiff"
            if config["renewable"][w.technology]["natura"]
            else []
        ),
        gebco=lambda w: (
            "data/bundle/GEBCO_2014_2D.nc"
            if "max_depth" in config["renewable"][w.technology].keys()
            else []
        ),
        ship_density=lambda w: (
            "resources/" + RDIR + "shipdensity_raster.nc"
            if "ship_threshold" in config["renewable"][w.technology].keys()
            else []
        ),
        country_shapes="resources/" + RDIR + "country_shapes.geojson",
        offshore_shapes="resources/" + RDIR + "offshore_shapes.geojson",
        regions=lambda w: (
            "resources/" + RDIR + "regions_onshore.geojson"
            if w.technology in ("onwind", "solar")
            else "resources/" + RDIR + "regions_offshore.geojson"
        ),
        cutouts=renewable_profiles_cutouts,
    output:
        profile="resources/" + RDIR + "profile_{technology}_{year}.nc",
    log:
        "logs/" + RDIR + "build_renewable_profile_{technology}_{year}.log",
    benchmark:
        "benchmarks/" + RDIR + "build_renewable_profiles_{technology}_{year}"
    threads: ATLITE_NPROCESSES
    resources:
        mem_mb=build_renewable_profiles_memory,
    wildcard_constraints:
        technology="(?!hydro).*",  # Any technology other than hydro
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/build_renewable_profiles.py"


ruleorder: build_hydro_profile > build_renewable_profiles


def hydro_profiles_cutouts(wildcards):
    if "hydro" in config["renewable"]:
        years = parse_year_wildcard(wildcards.year)
        # In case the year boundary is not the 1st of January, we also
        # need the cutout for year y+1 for each y in `years`.
        if config["snapshots"].get("year_boundary", "01-01") != "01-01":
            years = set(years + [y + 1 for y in years])
        return [
            "cutouts/" + CDIR + f"{config['renewable']['hydro']['cutout']}_{y}.nc"
            for y in years
        ]
    else:
        return "config['renewable']['hydro']['cutout'] not configured"


rule build_hydro_profile:
    input:
        country_shapes="resources/" + RDIR + "country_shapes.geojson",
        eia_hydro_generation="data/EIA_scaled_hydro_1980_2020.csv",
        cutouts=hydro_profiles_cutouts,
    output:
        "resources/" + RDIR + "profile_hydro_{year}.nc",
    log:
        "logs/" + RDIR + "build_hydro_profile_{year}.log",
    conda:
        "envs/environment.fixed.yaml"
    resources:
        mem_mb=build_renewable_profiles_memory,
    script:
        "scripts/build_hydro_profile.py"


def add_electricity_memory(wildcards):
    """Estimate memory requirements of `add_electricity`"""
    num_years = len(parse_year_wildcard(wildcards.year))
    return 3000 + num_years * 2000


rule add_electricity:
    input:
        **{
            f"profile_{tech}": "resources/" + RDIR + f"profile_{tech}" + "_{year}.nc"
            for tech in config["renewable"]
        },
        **{
            f"conventional_{carrier}_{attr}": fn
            for carrier, d in config.get("conventional", {None: {}}).items()
            for attr, fn in d.items()
            if str(fn).startswith("data/")
        },
        base_network="networks/" + RDIR + "base.nc",
        tech_costs=COSTS,
        regions="resources/" + RDIR + "regions_onshore.geojson",
        powerplants="resources/" + RDIR + "powerplants.csv",
        hydro_capacities="data/bundle/hydro_capacities.csv",
        geth_hydro_capacities="data/geth2015_hydro_capacities.csv",
        load="data/europe_demand_artificial_1980-2020.csv",
        nuts3_shapes="resources/" + RDIR + "nuts3_shapes.geojson",
    output:
        "networks/" + RDIR + "elec_{year}.nc",
    log:
        "logs/" + RDIR + "add_electricity_{year}.log",
    benchmark:
        "benchmarks/" + RDIR + "add_electricity_{year}"
    threads: 1
    resources:
        mem_mb=add_electricity_memory,
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/add_electricity.py"


def simplify_memory(wildcards):
    """Estimate memory requirements of `simplify_network`"""
    num_years = len(parse_year_wildcard(wildcards.year))
    return 2000 + num_years * 2000


rule simplify_network:
    input:
        network="networks/" + RDIR + "elec_{year}.nc",
        network_constant="networks/" + RDIR + f"elec_{config['net_clustering_year']}.nc",
        tech_costs=COSTS,
        regions_onshore="resources/" + RDIR + "regions_onshore.geojson",
        regions_offshore="resources/" + RDIR + "regions_offshore.geojson",
    output:
        network="networks/" + RDIR + "elec_{year}_s{simpl}.nc",
        # Note that the following output files should actually for
        # equal for all years, but we cannot remove the {year} wildcard
        # from the names since that would lead to output filename
        # conflicts between different years.
        regions_onshore="resources/"
        + RDIR
        + "regions_onshore_elec_{year}_s{simpl}.geojson",
        regions_offshore="resources/"
        + RDIR
        + "regions_offshore_elec_{year}_s{simpl}.geojson",
        busmap="resources/" + RDIR + "busmap_elec_{year}_s{simpl}.csv",
        connection_costs="resources/" + RDIR + "connection_costs_{year}_s{simpl}.csv",
    log:
        "logs/" + RDIR + "simplify_network/elec_{year}_s{simpl}.log",
    benchmark:
        "benchmarks/" + RDIR + "simplify_network/elec_{year}_s{simpl}"
    conda:
        "envs/environment.fixed.yaml"
    threads: 1
    resources:
        mem_mb=simplify_memory,
    script:
        "scripts/simplify_network.py"


rule cluster_network:
    input:
        network="networks/" + RDIR + "elec_{year}_s{simpl}.nc",
        network_constant="networks/"
        + RDIR
        + "elec_"
        + str(config["net_clustering_year"])
        + "_s{simpl}.nc",
        regions_onshore="resources/"
        + RDIR
        + "regions_onshore_elec_{year}_s{simpl}.geojson",
        regions_offshore="resources/"
        + RDIR
        + "regions_offshore_elec_{year}_s{simpl}.geojson",
        busmap=ancient("resources/" + RDIR + "busmap_elec_{year}_s{simpl}.csv"),
        custom_busmap=(
            "data/custom_busmap_elec_s{simpl}_{clusters}.csv"
            if config["enable"].get("custom_busmap", False)
            else []
        ),
        tech_costs=COSTS,
    output:
        network="networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}.nc",
        regions_onshore="resources/"
        + RDIR
        + "regions_onshore_elec_{year}_s{simpl}_{clusters}.geojson",
        regions_offshore="resources/"
        + RDIR
        + "regions_offshore_elec_{year}_s{simpl}_{clusters}.geojson",
        busmap="resources/" + RDIR + "busmap_elec_{year}_s{simpl}_{clusters}.csv",
        linemap="resources/" + RDIR + "linemap_elec_{year}_s{simpl}_{clusters}.csv",
    log:
        "logs/" + RDIR + "cluster_network/elec_{year}_s{simpl}_{clusters}.log",
    benchmark:
        "benchmarks/" + RDIR + "cluster_network/elec_{year}_s{simpl}_{clusters}"
    conda:
        "envs/environment.fixed.yaml"
    threads: 1
    resources:
        mem_mb=simplify_memory,
    script:
        "scripts/cluster_network.py"


rule add_extra_components:
    input:
        network="networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}.nc",
        tech_costs=COSTS,
    output:
        "networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}_ec.nc",
    log:
        "logs/" + RDIR + "add_extra_components/elec_{year}_s{simpl}_{clusters}.log",
    benchmark:
        "benchmarks/" + RDIR + "add_extra_components/elec_{year}_s{simpl}_{clusters}_ec"
    conda:
        "envs/environment.fixed.yaml"
    threads: 1
    resources:
        mem_mb=3000,
    script:
        "scripts/add_extra_components.py"


rule prepare_network:
    input:
        "networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}_ec.nc",
        tech_costs=COSTS,
    output:
        "networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        "logs/"
        + RDIR
        + "prepare_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.log",
    benchmark:
        (
            "benchmarks/"
            + RDIR
            + "prepare_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    conda:
        "envs/environment.fixed.yaml"
    threads: 1
    resources:
        mem_mb=4000,
    script:
        "scripts/prepare_network.py"


def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


rule solve_network:
    input:
        "networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        "results/networks/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    log:
        solver=normpath(
            "logs/"
            + RDIR
            + "solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_solver.log"
        ),
        python="logs/"
        + RDIR
        + "solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_python.log",
        memory="logs/"
        + RDIR
        + "solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_memory.log",
    benchmark:
        (
            "benchmarks/"
            + RDIR
            + "solve_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    conda:
        "envs/environment.fixed.yaml"
    threads: 4
    resources:
        mem_mb=memory,
    shadow:
        "minimal"
    script:
        "scripts/solve_network.py"


rule solve_operations_network:
    input:
        unprepared="networks/" + RDIR + "elec_{year}_s{simpl}_{clusters}_ec.nc",
        optimized="results/networks/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
    output:
        "results/networks/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op.nc",
    log:
        solver=normpath(
            "logs/"
            + RDIR
            + "solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_solver.log"
        ),
        python="logs/"
        + RDIR
        + "solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_python.log",
        memory="logs/"
        + RDIR
        + "solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_op_memory.log",
    benchmark:
        (
            "benchmarks/"
            + RDIR
            + "solve_operations_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}"
        )
    conda:
        "envs/environment.fixed.yaml"
    threads: 4
    resources:
        mem_mb=(lambda w: 5000 + 372 * int(w.clusters)),
    shadow:
        "minimal"
    script:
        "scripts/solve_operations_network.py"


rule plot_network:
    input:
        network="results/networks/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        tech_costs=COSTS,
    output:
        only_map="results/plots/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}.{ext}",
        ext="results/plots/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_ext.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_network/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_{ext}.log",
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/plot_network.py"


def input_make_summary(w):
    # It is mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return [COSTS] + expand(
        "results/networks/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        ll=ll,
        **{
            k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
            for k in ["simpl", "clusters", "opts"]
        }
    )


rule make_summary:
    input:
        input_make_summary,
    output:
        directory(
            "results/summaries/"
            + RDIR
            + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}"
        ),
    log:
        "logs/"
        + RDIR
        + "make_summary/elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.log",
    resources:
        mem_mb=500,
    conda:
        "envs/environment.fixed.yaml"
    script:
        "scripts/make_summary.py"


rule plot_summary:
    input:
        "results/summaries/"
        + RDIR
        + "elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}",
    output:
        "results/plots/"
        + RDIR
        + "summary_{summary}_elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_summary/{summary}_elec_{year}_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{ext}.log",
    conda:
        "envs/environment.fixed.yaml"
    resources:
        mem_mb=500,
    script:
        "scripts/plot_summary.py"


def input_plot_p_nom_max(w):
    return [
        (
            "results/networks/"
            + RDIR
            + "elec_{year}_s{simpl}{maybe_cluster}.nc".format(
                maybe_cluster=("" if c == "full" else ("_" + c)), **w
            )
        )
        for c in w.clusts.split(",")
    ]


rule plot_p_nom_max:
    input:
        input_plot_p_nom_max,
    output:
        "results/plots/"
        + RDIR
        + "elec_{year}_s{simpl}_cum_p_nom_max_{clusts}_{techs}_{country}.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_p_nom_max/elec_{year}_s{simpl}_{clusts}_{techs}_{country}_{ext}.log",
    conda:
        "envs/environment.fixed.yaml"
    resources:
        mem_mb=500,
    script:
        "scripts/plot_p_nom_max.py"
