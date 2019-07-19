configfile: "config.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    snapshots="[0-9]*",
    lv="[0-9\.]*",
    gap="[0-9\.]*",
    formulation="(kirchhoff|angles)"

rule solve_network:
    input: "networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    output: 
        network="results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.nc",
    log: 
        memory="benchmarks/memorylogger_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.log",
        times="benchmarks/times_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.csv"
    benchmark: "benchmarks/snakemake_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.txt"
    script: "scripts/solve.py"

rule solve_all:
    input:
        expand("results/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_{formulation}.nc",
                **config['scenario'])