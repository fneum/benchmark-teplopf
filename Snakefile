configfile: "config.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    snapshots="[0-9]*",
    lv="[0-9\.]*",
    gap="[+-e0-9\.]*",
    formulation="(kirchhoff|angles)"

rule solve_network:
    input: "networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    output: 
        network="results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.nc",
    log: 
        memory="benchmarks/memorylogger_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.log",
        times="benchmarks/times_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.csv",
        solver="benchmarks/solver_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.log"
    threads: config["solver"]["Threads"]+1
    benchmark: "benchmarks/snakemake_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.txt"
    script: "scripts/solve.py"

rule summarize_individual:
    input:
        network="results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.nc",
        memory="benchmarks/memorylogger_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.log",
        times="benchmarks/times_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.csv",
        solver="benchmarks/solver_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.log"
    output: "results/summaries/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.csv"
    script: "scripts/logs.py"

rule summarize_all:
    input:
        expand("results/summaries/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.csv",
                **config['scenario'])
    output: "results/summaries/joint_summary.csv"
    script: "scripts/summary.py"

rule solve_all:
    input:
        expand("results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.nc",
                **config['scenario'])