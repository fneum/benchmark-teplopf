import logging
logging.basicConfig(level=logging.WARNING)

import numpy as np
import re, os
from textwrap import dedent

def extract_from_log(logfile, solver_name):

    if solver_name == 'gurobi':
        return extract_from_log_gurobi(logfile)
    elif solver_name == 'cplex':
        return extract_from_log_cplex(logfile)
    else:
        raise NotImplementedError("Log extraction not defined for solver {}".format(solver_name))

# https://github.com/FRESNA/benchmark-lopf
def extract_from_log_cplex(logfile):

    #stdout = logs

    m = re.search("^Read time = ([\d\.]+) sec\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['read_time'] = float(m.group(1))

    m = re.search(r"^Variables\s+:\s+(\d+).*\n"
                r".*\n.*\n"
                r"Linear constraints\s+:\s+(\d+).*\n"
                r"\s+Nonzeros\s+:\s+(\d+)\n", stdout, flags=re.MULTILINE)
    if m is not None:
        nums["Vars"] = int(m.group(1))
        nums["Constrs"] = int(m.group(2))
        nums["NZs"] = int(m.group(3))

    m = re.search(r"^Presolve time = ([\d\.]+) sec\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums["presol_time"] = float(m.group(1))

    m = re.search(r"^Reduced LP has (\d+) rows, (\d+) columns, and (\d+) nonzeros\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['pConstrs'] = int(m.group(1))
        nums['pVars'] = int(m.group(2))
        nums['pNZs'] = int(m.group(3))
    else:
        nums['pConstrs'] = np.nan
        nums['pVars'] = np.nan
        nums['pNZs'] = np.nan

    m = re.search(r"^([ a-zA-Z]+) solved model\.", stdout, flags=re.MULTILINE)
    if m is not None:
        nums['solved_with'] = m.group(1).lower()

    return nums

# https://github.com/FRESNA/benchmark-lopf
def extract_from_log_gurobi(logfile):

    # Extract info from log file
    with open(log_fn) as log_fp:
        logs = log_fp.read()
    m = re.search(r"Presolve time: ([\d\.]+)s", logs)
    if m is not None:
        nums['presol_time'] = float(m.group(1))
        nums['sol_time'] -= nums['presol_time']
    m = re.search(r"Presolved: (\d+) rows, (\d+) columns, (\d+) nonzeros", logs)
    if m is not None:
        nums['pConstrs'] = int(m.group(1))
        nums['pVars'] = int(m.group(2))
        nums['pNZs'] = int(m.group(3))
    else:
        nums['pConstrs'] = np.nan
        nums['pVars'] = np.nan
        nums['pNZs'] = np.nan

    m = re.search(r"Solved with ([a-z ]+)", logs)
    if m is not None:
        nums['solved_with'] = m.group(1)

    nums['Constrs'] = model.NumConstrs
    nums['Vars'] = model.NumVars
    nums['NZs'] = model.NumNZs

    nums['objective'] = getattr(model, 'ObjVal', np.nan)

    return nums