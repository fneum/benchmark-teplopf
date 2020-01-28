configfile: "config.yaml"

wildcard_constraints:
    clusters="[0-9]+m?|all",
    ll="(v|c)([0-9\.]+|opt|all)|all",
    opts="[-+a-zA-Z0-9\.]*",
    snapshots="[0-9]*",
    lv="[0-9\.]*",
    gap="[+-e0-9\.]*",
    formulation="(kirchhoff|angles)",
    problem="(TEP|GTEP)",
    region="[a-zA-Z]*"
 
# note: two separate subworkflows are necessary since 
# each region is clustered to a specific number of nodes

subworkflow pypsaeur_germany:
    workdir: "pypsa-eur-germany"
    snakefile: "pypsa-eur-germany/Snakefile"
    configfile: "config_pypsaeur_germany.yaml"


subworkflow pypsaeur_baltics:
    workdir: "pypsa-eur-baltics"
    snakefile: "pypsa-eur-baltics/Snakefile"
    configfile: "config_pypsaeur_baltics.yaml"


rule prepare_network_germany:
    input: pypsaeur_germany("networks/elec_s_{clusters}_l{ll}_{opts}_ec.nc")
    output: "networks/germany/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    script: "scripts/prepare.py"


rule prepare_network_baltics:
    input: pypsaeur_baltics("networks/elec_s_{clusters}_l{ll}_{opts}_ec.nc")
    output: "networks/baltics/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    script: "scripts/prepare.py"


rule solve_network:
    input: "networks/{region}/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc"
    output: 
        network="results/{region}/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.nc",
    log: 
        memory="benchmarks/{region}/memorylogger_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.log",
        times="benchmarks/{region}/times_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.csv",
        solver="benchmarks/{region}/solver_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.log"
    threads: max(config["milp_solver"]["threads"], config["lp_solver"]["threads"])
    resources: mem=20000
    benchmark: "benchmarks/{region}/snakemake_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.txt"
    script: "scripts/solve.py"


rule summarize_individual:
    input:
        network="results/{region}/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.nc",
        memory="benchmarks/{region}/memorylogger_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.log",
        times="benchmarks/{region}/times_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.csv",
        solver="benchmarks/{region}/solver_elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.log"
    output: "results/{region}/summaries/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.csv"
    script: "scripts/logs.py"


rule prepare_all_networks:
    input:
        expand("networks/{region}/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc", **config['scenarios-space']),
        expand("networks/{region}/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}.nc", **config['scenarios-time'])


rule solve_all:
    input:
        expand("results/{region}/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.nc", **config['scenarios-space']),
        expand("results/{region}/networks/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.nc", **config['scenarios-time'])


rule summarize_all:
    input:
        expand("results/{region}/summaries/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.csv", **config['scenarios-space']),
        expand("results/{region}/summaries/elec_s_{clusters}_l{ll}_t{snapshots}_{opts}_lv{lv}_gap{gap}_p{problem}_{formulation}.csv", **config['scenarios-time'])
    output: "results/summaries/joint_summary.csv"
    script: "scripts/summary.py"