""""
Joins the summaries of individual runs into one pd.DataFrame and .csv file.
"""

import pandas as pd
import numpy as np

if __name__ == "__main__":

    df = None

    for scenario in snakemake.input:
        series = pd.read_csv(scenario, header=None, index_col=0).squeeze()
        series.name = scenario[18:-4]
        if df is None:
            df = pd.DataFrame(series)
        else:
            df = df.join(series)

    df = df.T

    dtypelist = [
        np.int,  # Constrs
        np.int,  # Vars
        np.int,  # NZs
        np.float,  # presol_time
        np.int,  # pConstrs
        np.int,  # pVars
        np.int,  # pNZs
        np.int,  # CVars
        np.int,  # IVars
        np.int,  # BVars
        np.int,  # pCVars
        np.int,  # PIVars
        np.int,  # PBVars
        np.str,  # solved_with
        np.float,  # obj_upper_bound
        np.float,  # obj_lower_bound
        np.float,  # gap
        np.float,  # matrix_range_lower
        np.float,  # matrix_range_upper
        np.float,  # obj_range_lower
        np.float,  # obj_range_upper
        np.float,  # bounds_range_lower
        np.float,  # bounds_range_upper
        np.int,  # N
        np.int,  # L
        np.int,  # T
        np.float,  # pot_circuits
        np.str,  # formulation
        np.float,  # target_gap
        np.int,  # threads
        np.int,  # walltime
        np.str,  # region
        np.float,  # co2limit
        np.float,  # new_circuits
        np.float,  # rel_exp_volume
        np.float,  # abs_exp_volume
        np.float,  # peak_mem
        np.float,
        np.float,
        np.float,
    ]
    dtypes = dict(zip(df.columns, dtypelist))
    df = df.astype(dtypes)

    df.to_csv(snakemake.output[0])
