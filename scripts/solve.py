import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level="INFO")

from vresutils.benchmark import memory_logger
from vresutils.benchmark import timer

import pypsa
import pandas as pd

solver_options = snakemake.config['solver']
solver_name = solver_options.pop('name')
solver_options['MIPGap'] = float(snakemake.wildcards.gap)

times = pd.Series()

with memory_logger(filename=snakemake.log.memory, interval=5.) as mem:

    for i in range(snakemake.config['iterations']):

        logger.info("Start benchmark iteration {}".format(i))

        with timer("loading network") as t:
            n = pypsa.Network(snakemake.input[0])
            n.lines.loc[n.lines.num_parallel==0., 'num_parallel'] = 1.
            n.lines.s_nom_extendable=True
            n.lines.s_nom_max = n.lines.s_nom + 1700 * float(snakemake.wildcards.lv)#n.lines.s_nom * float(snakemake.wildcards.lv)
            n.links.p_nom_max = n.links.p_snom + 4000
            n.global_constraints.constant = n.global_constraints.constant / 2
            n.calculate_dependent_values()
        times["loading {}".format(i)] = t.usec

        with timer("infer candidates") as t:
            n.lines = pypsa.tepopf.infer_candidates_from_existing(n)
        times["candidates {}".format(i)] = t.usec

        with timer("building pyomo model") as t:
            pypsa.tepopf.network_teplopf_build_model(n,
                formulation=snakemake.wildcards.formulation)
        times["building {}".format(i)] = t.usec

        with timer("solving model") as t:
            pypsa.opf.network_lopf_prepare_solver(n, solver_name)
            pypsa.opf.network_lopf_solve(n, formulation=snakemake.wildcards.formulation,
                solver_options=solver_options,
                solver_logfile=snakemake.log.solver)
        times["solving {}".format(i)] = t.usec

        with timer("exporting network") as t:
            n.export_to_netcdf(snakemake.output[0])
        times["exporting {}".format(i)] = t.usec

    df = pd.DataFrame(dict(step=times.index, secs=times.values/1e6))
    df.to_csv(snakemake.log.times, index=False)

logger.info("Maximum memory usage: {}".format(mem.mem_usage))