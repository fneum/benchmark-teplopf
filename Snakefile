configfile: "config.yaml"

wildcard_constraints:
    clusters="[0-9]+m?",
    opts="[-+a-zA-Z0-9\.]*",
    snapshots="[0-9]*",
    candidates="[0-9\.]*",
    gap="[+-e0-9\.]*",
    formulation="(kirchhoff|angles)",
    problem="(TEP|GTEP)",
    region="[-a-zA-Z]*"
 

subworkflow pypsaeur:
    workdir: "pypsa-eur"
    snakefile: "pypsa-eur/Snakefile"
    configfile: "config.pypsaeur.yaml"


rule prepare_network:
    input: pypsaeur("networks/elec_s_{clusters}_ec_lcopt_{opts}.nc")
    output: "networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}.nc"
    script: "scripts/prepare.py"


rule solve_network:
    input: "networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}.nc"
    output: 
        network="results/networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.nc",
    log: 
        memory="benchmarks/memorylogger_elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.log",
        times="benchmarks/times_elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.csv",
        solver="benchmarks/solver_elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.log"
    threads: max(config["milp_solver"]["threads"], config["lp_solver"]["threads"])
    resources: mem=20000
    benchmark: "benchmarks/snakemake_elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.txt"
    script: "scripts/solve.py"


rule summarize_individual:
    input:
        network="results/networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.nc",
        memory="benchmarks/memorylogger_elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.log",
        times="benchmarks/times_elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.csv",
        solver="benchmarks/solver_elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.log"
    output: "results/summaries/elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.csv"
    script: "scripts/logs.py"


rule prepare_all_networks:
    input:
        expand("networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}.nc", **config['scenarios-space']),
        expand("networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}.nc", **config['scenarios-time'])


rule solve_all:
    input:
        expand("results/networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.nc", **config['scenarios-space']),
        expand("results/networks/elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.nc", **config['scenarios-time'])


rule summarize_all:
    input:
        expand("results/summaries/elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.csv", **config['scenarios-space']),
        expand("results/summaries/elec_s_{clusters}_t{snapshots}_r{region}_{opts}_c{candidates}_g{gap}_p{problem}_{formulation}.csv", **config['scenarios-time'])
    output: "results/summaries/joint_summary.csv"
    script: "scripts/summary.py"