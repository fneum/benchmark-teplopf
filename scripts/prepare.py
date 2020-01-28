"""
Prepares a network for solving, including filtering regions, time series aggregation and load shedding.
"""

import pypsa
import pandas as pd
import tsam.timeseriesaggregation as tsam
import logging

logger = logging.getLogger(__name__)
logger.setLevel("INFO")


def filter_countries(n, countries=["DE"]):

    if type(countries) is str:
        countries = [countries]

    def filter_branch(x):
        return x.bus1[:2] not in countries or x.bus0[:2] not in countries

    def filter_oneport(x):
        return x.bus[:2] not in countries

    def filter_bus(x):
        return x.name[:2] not in countries

    to_delete = n.storage_units.loc[n.storage_units.apply(filter_oneport, axis=1)].index
    n.mremove(
        "StorageUnit",
        n.storage_units.loc[n.storage_units.apply(filter_oneport, axis=1)].index,
    )

    to_delete = n.generators.loc[n.generators.apply(filter_oneport, axis=1)].index
    n.mremove("Generator", to_delete)

    to_delete = n.loads.loc[n.loads.apply(filter_oneport, axis=1)].index
    n.mremove("Load", to_delete)

    to_delete = n.links.loc[n.links.apply(filter_branch, axis=1)].index
    n.mremove("Link", to_delete)

    to_delete = n.lines.loc[n.lines.apply(filter_branch, axis=1)].index
    n.mremove("Line", to_delete)

    to_delete = n.buses.loc[n.buses.apply(filter_bus, axis=1)].index
    n.mremove("Bus", to_delete)

    return n


def cluster_snapshots(n, noTypicalPeriods, hoursPerPeriod):

    # prepare time series for aggregation
    timeseries = pd.concat(
        [n.generators_t.p_max_pu, n.loads_t.p_set], axis=1, sort=False
    )

    # perform aggregation
    aggregation = tsam.TimeSeriesAggregation(
        timeseries, noTypicalPeriods=noTypicalPeriods, hoursPerPeriod=hoursPerPeriod
    )
    typPeriods = aggregation.createTypicalPeriods()

    # process output of aggregation to fit pypsa format
    if hoursPerPeriod == 1:
        new_index = ["H{}".format(str(t).zfill(4)) for t in range(0, noTypicalPeriods)]
    else:
        new_index = [
            "D{} {}:00".format(str(d).zfill(2), str(h).zfill(2))
            for d in range(0, noTypicalPeriods)
            for h in range(0, hoursPerPeriod)
        ]

    typPeriods = typPeriods.reset_index(drop=True)
    typPeriods = typPeriods.rename(dict(zip(typPeriods.index, new_index)), axis="index")
    cluster_sizes = list(aggregation.clusterPeriodNoOccur.values())
    new_weightings = [
        float(cluster_sizes[i])
        for i in range(0, len(cluster_sizes))
        for dupl in range(0, hoursPerPeriod)
    ]

    # overwrite time series data of network
    n.generators_t.p_max_pu = typPeriods.loc[:, n.generators_t.p_max_pu.columns]
    n.loads_t.p_set = typPeriods.loc[:, n.loads_t.p_set.columns]
    n.snapshots = new_index
    n.snapshot_weightings = pd.Series(new_weightings, new_index, name="weightings")
    n.snapshot_weightings.index.name = "name"

    assert (
        n.snapshot_weightings.sum() == 8760.0
    ), "Snapshots do not add up to 8760 hours."

    if len(n.storage_units.index) != 0:
        logger.warning(
            "Removing storage units since they are not supported with time series clustering."
        )
        n.mremove("StorageUnit", n.storage_units.index)

    if len(n.stores.index) != 0:
        logger.warning(
            "Removing stores since they are not supported with time series clustering."
        )
        n.mremove("Store", n.store.index)

    return n


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])

    countries = snakemake.wildcards.countries.split("-")

    filter_countries(n, countries=countries)

    # >>> workaround fixed in https://github.com/PyPSA/pypsa-eur/commit/97985e3c5219a75392291125914225d116ba81e9
    n.lines.loc[n.lines.num_parallel == 0.0, "num_parallel"] = 1.0
    n.lines.s_nom_extendable = True
    n.calculate_dependent_values()
    # <<<

    target_snapshots = snakemake.wildcards.snapshots

    if not (target_snapshots == "8760" or target_snapshots == ""):
        logger.info("Clustering is required.")
        n = cluster_snapshots(n, int(target_snapshots), 1)

    if snakemake.config["load_shedding"]:
        n.add("Carrier", "load")
        n.madd(
            "Generator",
            n.buses.index,
            " load",
            bus=n.buses.index,
            carrier="load",
            marginal_cost=1e5,
            p_nom=1e6,
        )

    n.export_to_netcdf(snakemake.output[0])
