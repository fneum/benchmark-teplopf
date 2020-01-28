"""
Extracts logged data from multiple sources such as from pypsa.Network and the solver logfile.
""""

import logging
logging.basicConfig(level=logging.WARNING)

import numpy as np
import re, os
from textwrap import dedent

import pandas as pd
import pypsa

def extract_from_log(logfile, solver_name):

    if solver_name == 'gurobi':
        return extract_from_log_gurobi(logfile)
    else:
        raise NotImplementedError("Log extraction not defined for solver {}".format(solver_name))

# adapted from https://github.com/FRESNA/benchmark-lopf
def extract_from_log_gurobi(logfile):

    stats = {}

    # Extract info from log file
    with open(logfile) as log_fp:
        logs = log_fp.read()
    
    m = re.search(r"Optimize a model with (\d+) rows, (\d+) columns and (\d+) nonzeros", logs)
    if m is not None:
        stats['Constrs'] = int(m.group(1))
        stats['Vars'] = int(m.group(2))
        stats['NZs'] = int(m.group(3))
    else:
        stats['Constrs'] = np.nan
        stats['Vars'] = np.nan
        stats['NZs'] = np.nan

    m = re.search(r"Presolve time: ([\d\.]+)s", logs)
    if m is not None:
        stats['presol_time'] = float(m.group(1))
    else:
        stats['presol_time'] = np.nan

    m = re.search(r"Presolved: (\d+) rows, (\d+) columns, (\d+) nonzeros", logs)
    if m is not None:
        stats['pConstrs'] = int(m.group(1))
        stats['pVars'] = int(m.group(2))
        stats['pNZs'] = int(m.group(3))
    else:
        stats['pConstrs'] = np.nan
        stats['pVars'] = np.nan
        stats['pNZs'] = np.nan

    ms = re.findall(r"Variable types: (\d+) continuous, (\d+) integer \((\d+) binary\)", logs)
    if ms is not None:
        pre = ''
        for m in ms:
            stats['{}CVars'.format(pre)] = int(m[0])
            stats['{}IVars'.format(pre)] = int(m[1])
            stats['{}BVars'.format(pre)] = int(m[2])
            pre = 'p'
    else:
        stats['CVars'] = np.nan
        stats['IVars'] = np.nan
        stats['BVars'] = np.nan
        stats['pCVars'] = np.nan
        stats['pIVars'] = np.nan
        stats['pBVars'] = np.nan

    m = re.search(r"Solved with ([a-z ]+)", logs)
    if m is not None:
        stats['solved_with'] = m.group(1)
    else:
        stats['solved_with'] = np.nan

    rfloate = "([+-]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"
    m = re.search(r"Best objective {r}, best bound {r}, gap {r}%".format(r=rfloate), logs)
    if m is not None:
        stats['obj_upper_bound'] = float(m.group(1))
        stats['obj_lower_bound'] = float(m.group(3))
        stats['gap'] = float(m.group(5)) / 100
    else:
        stats['obj_upper_bound'] = np.nan
        stats['obj_lower_bound'] = np.nan
        stats['gap'] = np.nan

    rinte = "([0-9]+[e][+-][0-9]+)"
    m = re.search(r"  Matrix range     \[{r}, {r}\]".format(r=rinte), logs)
    if m is not None:
        stats['matrix_range_lower'] = float(m.group(1))
        stats['matrix_range_upper'] = float(m.group(2))
    m = re.search(r"  Objective range  \[{r}, {r}\]".format(r=rinte), logs)
    if m is not None:
        stats['obj_range_lower'] = float(m.group(1))
        stats['obj_range_upper'] = float(m.group(2))
    m = re.search(r"  Bounds range     \[{r}, {r}\]".format(r=rinte), logs)
    if m is not None:
        stats['bounds_range_lower'] = float(m.group(1))
        stats['bounds_range_upper'] = float(m.group(2))
    m = re.search(r"  RHS range        \[{r}, {r}\]".format(r=rinte), logs)
    if m is not None:
        stats['obj_range_lower'] = float(m.group(1))
        stats['obj_range_upper'] = float(m.group(2))

    return stats

if __name__ == '__main__':

    network = pypsa.Network(snakemake.input.network)

    solve_conf = snakemake.config['milp_solver']
    solver_name = solve_conf['name']

    stats = extract_from_log(snakemake.input.solver, solver_name)

    stats['N'] = len(network.buses)
    stats['L'] = len(network.lines)
    stats['T'] = len(network.snapshots)
    stats['pot_circuits'] = snakemake.wildcards.candidates
    stats['formulation'] = snakemake.wildcards.formulation
    stats['target_gap'] = snakemake.wildcards.gap
    stats['threads'] = solve_conf['threads']
    stats['walltime'] = solve_conf['TimeLimit']
    stats['region'] = snakemake.wildcards.region

    for o in snakemake.wildcards.opts.split('-'):
        if 'Co2L' in o:
            m = re.findall("[0-9]*\.?[0-9]+$", o)
            stats['co2limit'] = float(m[0])

    clines = network.lines.loc[~network.lines.operative]
    elines = network.lines.loc[network.lines.operative]
    ctwkm = (clines.length * clines.s_nom_opt).sum() / 1e6
    etwkm = (elines.length * elines.s_nom_opt).sum() / 1e6

    stats['new_circuits'] = len(clines.loc[clines.s_nom_opt>0])
    stats['rel_exp_volume'] = ctwkm / etwkm
    stats['abs_exp_volume'] = ctwkm

    peak_mem = pd.read_csv(snakemake.input.memory,
                           sep=' ', index_col=1, usecols=[1,2],
                           header=None).max().values[0]
    stats['peak_mem'] = peak_mem
    
    stats_df = pd.Series(stats)

    times_df = pd.read_csv(snakemake.input.times, index_col=0).squeeze()

    stats_df = stats_df.append(times_df)

    stats_df.to_csv(snakemake.output[0])