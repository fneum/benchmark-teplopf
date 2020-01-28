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
        np.int,
        np.int,
        np.int,
        np.float,
        np.int,
        np.int,
        np.int,
        np.int,
        np.int,
        np.int,
        np.int,
        np.int,
        np.int,
        np.str,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.int,
        np.int,
        np.int,
        np.float,
        np.str,
        np.float,
        np.int,
        np.int,
        np.int,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
        np.float,
    ]
    dtypes = dict(zip(df.columns, dtypelist))
    df = df.astype(dtypes)

    df.to_csv(snakemake.output[0])
