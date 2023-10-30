import pandas as pd
import csv

out_dir = "output.roux-2017"
logs_dir = f'{out_dir}/logs'
basename="roux-2017"
db_basename = "RefSeq_v219"
db_fasta = '/home/amhorst/2023-sourmash-viruses/RefSeq_v219.fa'

samples = 'inputs/roux-2017.samples.csv'
data_dir = '/home/amhorst/2023-sourmash-viruses/test_dataset_Roux2017/raw_reads'

params = {"dna": {"ksize": [7,9,11,13,15,17,19,21,23,25,27,29,31], "scaled": 1, "threshold_bp": 0},
          "protein": {"ksize": [5, 7], "scaled": 1, "threshold_bp": 0}} #, 6,10

TAX_FILE = 'inputs/vmr_MSL38_v1.taxonomy.csv' # ICTV VMR taxonomy file (replace with RefSeq_v219 taxonomy file)

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
prot_params = expand(f"{basename}-x-{db_basename}.protein.k{{k}}-sc{{scaled}}.t{{thresh}}", k=params['protein']['ksize'], 
                                                                                            scaled=params['protein']['scaled'],
                                                                                            thresh=params['protein']['threshold_bp'])
# all_params = prot_params
all_params = dna_params #+ prot_params


rule all:
    input:
        expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.csv", search_params = all_params), 
        # expand(f"{out_dir}/fastmultigather/{{search_params}}.gather.with-lineages.csv", search_params = all_params), 


rule sketch_database:
    input:
        input_fasta = db_fasta,
    output:
        db_zip = "refseq-dbs/{db_basename}.{moltype}.zip",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/sketch/{{db_basename}}.{{moltype}}.sketch.log"
    params:
    shell:
        """
        sourmash scripts manysketch -p {param_str} -o {output.db_zip} \
                                   {input.input_fasta} 2> {log}
        """

rule index_database:
    input:
        db_zip = "refseq-dbs/{db_basename}.{moltype}.zip"
    output:
        current_file = "refseq-dbs/{db_basename}.{moltype}-k{ksize}.rocksdb/CURRENT",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/index/{{db_basename}}.{{moltype}}-k{{ksize}}.index.log"
    params:
        db_rocksdb = "refseq-dbs/{db_basename}.{moltype}-k{ksize}.rocksdb",
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
    log: f"{logs_dir}/sketch/{{basename}}.{{moltype}}.sketch.log"
    shell:
        """
        sourmash scripts manysketch -p {threads} -o {output} \
                                    {input.sample_csv} 2> {log}
        """

rule sourmash_fastmultigather:
    input:
        query_zip = f"{out_dir}/{{basename}}.{{moltype}}.zip",
        database  = "refseq-dbs/{db_basename}.{moltype}-k{ksize}.rocksdb/CURRENT",
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
