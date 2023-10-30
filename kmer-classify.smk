import pandas as pd
import csv

# load from configfile
configfile: "inputs/roux2017.config.yml"

basename = config['basename']
out_dir = config['output_dir']
db_basename = config['database']['basename'] 
db_fasta = config['database']['fromfile_csv'] 
TAX_FILE = config['database']['taxonomy']
samples = config['samples_csv']
params = config['sourmash_params']

logs_dir = f'{out_dir}/logs'

onstart:
    print("------------------------------")
    print("sourmash taxonomic classification workflow")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")

#build parameter combinations for easy expansion below
dna_params = expand(f"{basename}-x-{db_basename}.dna.k{{k}}-sc{{scaled}}.t{{thresh}}", k=params['dna']['ksize'], 
                                                                                       scaled=params['dna']['scaled'],
                                                                                       thresh=params['dna']['threshold_bp'])
#prot_params = expand(f"{basename}-x-{db_basename}.protein.k{{k}}-sc{{scaled}}.t{{thresh}}", k=params['protein']['ksize'], 
                                                                      #                      scaled=params['protein']['scaled'],
                                                                      #                      thresh=params['protein']['threshold_bp'])
all_params = dna_params #+ prot_params

rule all:
    input:
        expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.csv", search_params = all_params), 
        # expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.with-lineages.csv", search_params = all_params), 


def build_param_str(moltype):
    if moltype not in ['dna', 'protein']:
        raise ValueError(f"moltype {moltype} not recognized")
    ksizes = params[moltype]['ksize']
    scaled = min(params[moltype]['scaled'])
    k_params = ",".join([f"k={k}" for k in ksizes])
    param_str = f"{moltype},{k_params},scaled={scaled},abund"
    return param_str


rule sketch_database:
    input:
        input_fasta = db_fasta,
    output:
        db_zip = "sourmash-db/{db_basename}.{moltype}.zip",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/sketch/{{db_basename}}.{{moltype}}.sketch.log"
    threads: 1 # if single file, manysketch doesn't actually parallelize
    params:
        param_str = lambda w: build_param_str(w.moltype),
    shell:
        """
        sourmash scripts manysketch -p {params.param_str} -o {output.db_zip} \
                                   {input.input_fasta} 2> {log}
        """

rule index_database:
    input:
        db_zip = "sourmash-db/{db_basename}.{moltype}.zip"
    output:
        current_file = "sourmash-db/{db_basename}.{moltype}-k{ksize}.rocksdb/CURRENT",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/index/{{db_basename}}.{{moltype}}-k{{ksize}}.index.log"
    params:
        db_rocksdb = "sourmash-db/{db_basename}.{moltype}-k{ksize}.rocksdb",
    threads: 1
    shell:
        """
        sourmash scripts index {input.db_zip} -m {wildcards.moltype} \
                               -k {wildcards.ksize} --scaled 1 -o {params.db_rocksdb} 2> {log}
        """

rule sketch_samples:
    input: 
        sample_csv = samples,
    output:
        f"{out_dir}/{{basename}}.{{moltype}}.zip",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    threads: 1,
    log: f"{logs_dir}/sketch/{{basename}}.{{moltype}}.sketch.log"
    shell:
        """
        sourmash scripts manysketch -p {threads} -o {output} \
                                    {input.sample_csv} 2> {log}
        """

rule sourmash_fastmultigather:
    input:
        query_zip = f"{out_dir}/{{basename}}.{{moltype}}.zip",
        database  = "sourmash-db/{db_basename}.{moltype}-k{ksize}.rocksdb/CURRENT",
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        time=60,
        partition="low2",
    threads: 100
    log: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    params:
        db_dir = lambda w: f'output.vmr/{w.db_basename}.{w.moltype}-k{w.ksize}.rocksdb',
    shell:
        """
        echo "DB: {params.db_dir}"
        echo "DB: {params.db_dir}" > {log}

        sourmash scripts fastmultigather --threshold {wildcards.threshold} \
                          --moltype {wildcards.moltype} --ksize {wildcards.ksize} \
                          --scaled {wildcards.scaled} {input.query_zip} \
                          {params.db_dir} -o {output} 2>> {log}
        """


## to do: make lineages file for refseq db (or will full genbank one work?)

rule tax_annotate:
    input:
        gather_csv = f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
        lineages = TAX_FILE,
    output:
        f"{out_dir}/fastmultigather/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.with-lineages.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/annotate/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    params:
        output_dir = f"{out_dir}/fastmultigather",
    shell:
        """
        sourmash tax annotate -g {input.gather_csv} \
                              -t {input.lineages} \
                              -o {params.output_dir} 2> {log}
        """
