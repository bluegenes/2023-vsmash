import pandas as pd
import csv
from collections import defaultdict
from snakemake.shell import shell

# load from configfile
configfile: "inputs/roux2017.config.yml"

basename = config['basename']
out_dir = config['output_dir']
db_basename = config['database']['basename'] 
dna_zip = config['database'].get('dna_zip')
protein_zip = config['database'].get('protein_zip')
db_fasta = config['database'].get('fasta')

if dna_zip and not os.path.exists(f"sourmash-db/{db_basename}.dna.zip"):
    # copy file
    shell(f"cp {dna_zip} sourmash-db/{db_basename}.dna.zip")
if protein_zip and not os.path.exists(f"sourmash-db/{db_basename}.protein.zip"):
    shell(f"{protein_zip} sourmash-db/{db_basename}.protein.zip")

db_fromfile = config['database']['fromfile_csv']
TAX_FILE = config['database']['taxonomy']
samples = config['samples_csv']

sampleDF = pd.read_csv(samples)
sample_names = sampleDF['name'].tolist()
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
dna_params = expand(f"{db_basename}.dna.k{{k}}-sc{{scaled}}.t{{thresh}}", k=params['dna']['ksize'], 
                                                                                       scaled=params['dna']['scaled'],
                                                                                       thresh=params['dna']['threshold_bp'])
#prot_params = expand(f"{db_basename}.protein.k{{k}}-sc{{scaled}}.t{{thresh}}", k=params['protein']['ksize'], 
                                                                      #                      scaled=params['protein']['scaled'],
                                                                      #                      thresh=params['protein']['threshold_bp'])
all_params = dna_params #+ prot_params


#wildcard_constraints:
    # basename = '\w+',
   # db_basename = '\w+',


rule all:
    input:
        expand(f"{out_dir}/fastmultigather/{basename}-x-{{search_params}}.gather.csv", search_params = all_params),
        # add sample gather
        # expand(f"{out_dir}/databases/{{sample}}-x-{{search_params}}.zip", sample = sample_names,
        #                                                                   search_params = all_params),
        expand(f"{out_dir}/gather/{{sample}}-x-{{search_params}}.gather.csv", sample = sample_names,
                                                                              search_params = all_params),
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
        fasta = db_fasta,
#        fromfile = db_fromfile,
    output:
        db_zip = "sourmash-db/{db_basename}.{moltype}.zip",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/sketch/{{db_basename}}.{{moltype}}.sketch.log"
    benchmark: f"{logs_dir}/sketch/{{db_basename}}.{{moltype}}.sketch.benchmark"
    threads: 1 # if single file, manysketch doesn't actually parallelize
    params:
        param_str = lambda w: build_param_str(w.moltype),
    shell:
        """
        sourmash sketch dna --singleton -p {params.param_str} -o {output.db_zip} \
                             {input.fasta} 2> {log}
        """

rule index_database:
    input:
        db_zip = ancient("sourmash-db/{db_basename}.{moltype}.zip"),
    output:
        current_file = "sourmash-db/{db_basename}.{moltype}-k{ksize}.rocksdb/CURRENT",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/index/{{db_basename}}.{{moltype}}.k{{ksize}}.index.log"
    benchmark: f"{logs_dir}/index/{{db_basename}}.{{moltype}}.k{{ksize}}.index.benchmark"
    params:
        db_rocksdb = "sourmash-db/{db_basename}.{moltype}-k{ksize}.rocksdb",
        moltype = lambda w: w.moltype.upper() if w.moltype == 'dna' else w.moltype,
        scaled = lambda w: min(params[w.moltype]['scaled']),
    threads: 1
    shell:
        """
        sourmash scripts index {input.db_zip} -m {params.moltype} \
                               -k {wildcards.ksize} --scaled {params.scaled} \
                               -o {params.db_rocksdb} 2> {log}
        """

rule sketch_samples:
    input: 
        sample_csv = samples,
    output:
        f"{out_dir}/{basename}.{{moltype}}.zip",
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    threads: 14, # one per sample
    log: f"{logs_dir}/sketch/{basename}.{{moltype}}.sketch.log"
    benchmark: f"{logs_dir}/sketch/{basename}.{{moltype}}.sketch.benchmark"
    params:
        param_str = lambda w: build_param_str(w.moltype),
    shell:
        """
        sourmash scripts manysketch -p {params.param_str} -o {output} \
                                    {input.sample_csv} -c {threads} 2> {log}
        """

rule sourmash_fastmultigather:
    input:
        query_zip = ancient(f"{out_dir}/{basename}.{{moltype}}.zip"),
        database  = "sourmash-db/{db_basename}.{moltype}-k{ksize}.rocksdb/CURRENT",
    output:
        f"{out_dir}/fastmultigather/{basename}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        time=60,
        partition="low2",
    threads: 100
    log: f"{logs_dir}/fastmultigather/{basename}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/fastmultigather/{basename}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    params:
        db_dir = lambda w: f'sourmash-db/{w.db_basename}.{w.moltype}-k{w.ksize}.rocksdb',
        moltype = lambda w: w.moltype.upper() if w.moltype == 'dna' else w.moltype,
    shell:
        """
        echo "DB: {params.db_dir}"
        echo "DB: {params.db_dir}" > {log}

        sourmash scripts fastmultigather --threshold {wildcards.threshold} \
                          --moltype {params.moltype} --ksize {wildcards.ksize} \
                          --scaled {wildcards.scaled} {input.query_zip} -c {threads} \
                          {params.db_dir} -o {output} 2>> {log}
        """

rule split_results_by_sample:
    input: 
        sample_csv = samples,
        gather_csv = f"{out_dir}/fastmultigather/{basename}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
    output:
        sample_fmg = expand("{out_dir}/fastmultigather/sample/{sample}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv", out_dir = out_dir,sample = sample_names)
    run:
        #make dictionary of sample to filename
        sample_filemap = {}
        for fn in output.sample_fmg:
            sample = os.path.basename(fn).split('-x-')[0]
            sample_filemap[sample] = fn

        # split input csv into lists of rows per sample if match_name is equal to one of the sample names, then write each sample file
        with open(input.gather_csv) as f:
            reader = csv.DictReader(f)
            sample_rows = defaultdict(list)
            for row in reader:
                if row['query_name'] in sample_names:
                    sample_rows[row['query_name']].append(row)
            for sample in sample_names:
                with open(sample_filemap[sample], 'w') as out:
                    writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
                    writer.writeheader()
                    for row in sample_rows[sample]:
                        writer.writerow(row)

rule picklist_db:
    input:
        database_zip = ancient(f"sourmash-db/{db_basename}.{{moltype}}.zip"),
        sample_picklist =  f"{out_dir}/fastmultigather/sample/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
    output:
        picklist_zip = f"{out_dir}/databases/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.zip"
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    log: f"{logs_dir}/picklist/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/picklist/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    shell:
        """
        sourmash sig cat -o {output} --ksize {wildcards.ksize} \
                         --picklist {input.sample_picklist}:match_name:name \
                         {input.database_zip} 2> {log}
        """

rule sourmash_gather:
    input:
        query_zip = ancient(f"{out_dir}/{basename}.{{moltype}}.zip"),
        picklist_zip = f"{out_dir}/databases/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.zip",
        # database = "sourmash-db/{db_basename}.{moltype}.zip", 
        sample_picklist =  f"{out_dir}/fastmultigather/sample/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
    output:
        csv=f"{out_dir}/gather/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.csv",
        txt=f"{out_dir}/gather/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.gather.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        time=60,
        partition="low2",
    threads: 1
    log: f"{logs_dir}/gather/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.log"
    benchmark: f"{logs_dir}/gather/{{sample}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    # conda: "conf/env/branchwater.yml"
    conda: "pyo3-branch",
    params:
        alpha_cmd = lambda w: f"--{w.moltype}",
        name_regex = lambda w: f"{w.sample}[^\d]"
    shell:
        '''
        echo "DB: {input.picklist_zip}"
        echo "DB: {input.picklist_zip}" > {log}
        sourmash sig grep {params.name_regex:q} {input.query_zip} --ksize {wildcards.ksize} {params.alpha_cmd} | \
                          sourmash gather - {input.picklist_zip} \
                          --threshold {wildcards.threshold} \
                          {params.alpha_cmd} --ksize {wildcards.ksize} \
                          --scaled {wildcards.scaled}  \
                           -o {output.csv} > {output.txt} 2>> {log}
        '''
                          #--picklist {input.sample_picklist}:match_md5:md5 \
        # echo "DB: {input.database}"
        # echo "DB: {input.database}" > {log}
        #    echo "DB: {input.picklist_zip}"
        # echo "DB: {input.picklist_zip}" > {log}
        # sourmash sig extract --name "{wildcards.sample}" {input.query_zip} --ksize {wildcards.ksize} {params.alpha_cmd} | \
        #                   sourmash gather - {input.picklist_zip} \
        #                   --threshold {wildcards.threshold} \
        #                   {params.alpha_cmd} --ksize {wildcards.ksize} \
        #                   --scaled {wildcards.scaled}  \
        #                    -o {output} 2>> {log}


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
    benchmark: f"{logs_dir}/annotate/{{basename}}-x-{{db_basename}}.{{moltype}}.k{{ksize}}-sc{{scaled}}.t{{threshold}}.benchmark"
    params:
        output_dir = f"{out_dir}/fastmultigather",
    shell:
        """
        sourmash tax annotate -g {input.gather_csv} \
                              -t {input.lineages} \
                              -o {params.output_dir} 2> {log}
        """
