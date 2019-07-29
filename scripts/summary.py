import pandas as pd
import numpy as np

concat_list = []
for scenario in snakemake.input:
    series = pd.read_csv(scenario, header=None, index_col=0).squeeze()
    concat_list.append(series)

df = pd.condat(concat_list, axis=1)

df.to_csv(snakemake.output[0])
