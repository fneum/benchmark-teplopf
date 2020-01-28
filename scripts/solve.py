"""
Runs experiment and tracks time and memory. Substitutes interconnectors with HVAC lines.
"""

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level="INFO")

from vresutils.benchmark import memory_logger
from vresutils.benchmark import timer

import pypsa
import pandas as pd
import numpy as np


def substitute_interconnectors(n):
    # currently specific to lithuania and poland
    # for synchronizing other countries adaptations are necessary
    # TODO: fix

    links_to_drop = None
    synchronizers = None
    for c in range(1, snakemake.config["synchronizer_circuits"] + 1):
        fields = ["bus0", "bus1", "length", "capital_cost", "terrain_factor"]
        level_synchronizers = n.links.loc[
            [
                ((b0[:2] == "PL") & (b1[:2] == "LT"))
                | ((b1[:2] == "PL") & (b0[:2] == "LT"))
                for b0, b1 in zip(n.links.bus0, n.links.bus1)
            ]
        ]

        if links_to_drop is None:
            links_to_drop = level_synchronizers.index

        level_synchronizers["type"] = "Al/St 240/40 4-bundle 380.0"
        level_synchronizers["v_nom"] = 380.0
        level_synchronizers["s_nom"] = (
            c
            * np.sqrt(3)
            * level_synchronizers.type.map(n.line_types.i_nom)
            * level_synchronizers.v_nom
        )
        level_synchronizers["num_parallel"] = float(c)
        level_synchronizers.index = level_synchronizers.apply(
            lambda l: "sync-c{}-{}-{}".format(c, l.bus0, l.bus1), axis=1
        )
        synchronizers = (
            synchronizers.append(level_synchronizers)
            if synchronizers is not None
            else level_synchronizers
        )
    synchronizers["operative"] = False
    synchronizers["s_nom_extendable"] = True
    synchronizers["s_max_pu"] = 0.7
    synchronizers["s_nom_max"] = synchronizers.s_nom
    synchronizers["s_nom_min"] = 0.0

    n.lines = n.lines.append(synchronizers)
    n.links.drop(links_to_drop, axis=0, inplace=True)

    n.calculate_dependent_values()

    return n.lines, n.links


if __name__ == "__main__":

    lp_solver_options = snakemake.config["lp_solver"]
    lp_solver_name = lp_solver_options.pop("name")

    milp_solver_options = snakemake.config["milp_solver"]
    milp_solver_name = milp_solver_options.pop("name")
    milp_solver_options["MIPGap"] = float(snakemake.wildcards.gap)

    times = pd.Series()

    with memory_logger(filename=snakemake.log.memory, interval=5.0) as mem:

        for i in range(snakemake.config["iterations"]):

            logger.info("Start benchmark iteration {}".format(i))

            with timer("loading network") as t:
                n = pypsa.Network(snakemake.input[0])
                unit_s_nom = (
                    np.sqrt(3) * n.lines.type.map(n.line_types.i_nom) * n.lines.v_nom
                ).fillna(1700.0)
                n.lines.s_nom_max = n.lines.s_nom + unit_s_nom * float(
                    snakemake.wildcards.candidates
                )
                n.links.p_nom_max = n.links.p_nom + snakemake.config[
                    "link_unit_p_nom"
                ] * float(snakemake.wildcards.candidates)
            times["loading {}".format(i)] = t.usec

            if snakemake.wildcards.problem == "TEP":
                with timer("prepare tep") as t:
                    n.lopf(
                        solver_name=lp_solver_name,
                        solver_options=lp_solver_options,
                        formulation="kirchhoff",
                    )
                    # fix generation and storage fleet
                    n.generators.p_nom = n.generators.p_nom_opt
                    n.generators.p_nom_extendable = False
                    n.storage_units.p_nom = n.storage_units.p_nom_opt
                    n.storage_units.p_nom_extendable = False
                    n.stores.e_nom = n.stores.e_nom_opt
                    n.stores.e_nom_extendable = False
                    n.links.loc[n.links.carrier != "DC", "p_nom"] = n.links.loc[
                        n.links.carrier != "DC", "p_nom_opt"
                    ]
                    n.links.loc[n.links.carrier != "DC", "p_nom_extendable"] = False
                    line_expansion_volume = (
                        (n.lines.s_nom_opt - n.lines.s_nom) * n.lines.length
                    ).sum() / 1e6  # TWkm
                    link_expansion_volume = (
                        (n.links.p_nom_opt - n.links.p_nom) * n.links.length
                    ).sum() / 1e6  # TWkm
                    logger.info(
                        "Volume of continuous transmission line expansion is {0:.2f} TWkm for links and {1:.2f} TWkm for lines.".format(
                            link_expansion_volume, line_expansion_volume
                        )
                    )
                times["prepare tep {}".format(i)] = t.usec

            with timer("infer candidates") as t:
                n.lines = pypsa.tepopf.infer_candidates_from_existing(n)
                n.lines, n.links = substitute_interconnectors(n)
            times["candidates {}".format(i)] = t.usec

            with timer("building pyomo model") as t:
                pypsa.tepopf.network_teplopf_build_model(
                    n, formulation=snakemake.wildcards.formulation
                )
            times["building {}".format(i)] = t.usec

            with timer("solving model") as t:
                pypsa.opf.network_lopf_prepare_solver(n, solver_name=milp_solver_name)
                status, termination_condition = pypsa.opf.network_lopf_solve(
                    n,
                    formulation=snakemake.wildcards.formulation,
                    solver_options=milp_solver_options,
                    solver_logfile=snakemake.log.solver,
                    candidate_lines=True,
                )
            times["solving {}".format(i)] = t.usec

            with timer("exporting network") as t:
                n.export_to_netcdf(snakemake.output[0])
            times["exporting {}".format(i)] = t.usec

        df = pd.DataFrame(dict(step=times.index, secs=times.values / 1e6))
        df.to_csv(snakemake.log.times, index=False)

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
