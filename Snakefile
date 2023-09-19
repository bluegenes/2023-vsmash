import os
import pandas as pd

outdir = config.get('outdir', 'output.vsmash')
logs_dir = os.path.join(outdir, 'logs')

# read in excel sheet to get sample accessions
SAMPLES = ["SRR3458562", "SRR3458563", "SRR3458564", "SRR3458565", "SRR3458566", "SRR3458567", "SRR3458568", "SRR3458569"]

rule all:
    input:
        expand(os.path.join(outdir, 'raw', '{sample}_1.fastq.gz'), sample=SAMPLES),

####################
# download reads #
####################

rule kingfisher_download:
    output:
        r1=os.path.join(outdir, 'raw', '{sample}_1.fastq.gz'),
        r2=os.path.join(outdir, 'raw', '{sample}_2.fastq.gz'),
    threads: 10,
    resources:
        mem_mb = 3000,
        time = 100,
        partition = 'high2',
    conda: "conf/env/kingfisher.yml",
    params:
        kf_dir = os.path.join(outdir, 'raw'),
    log: os.path.abspath(os.path.join(logs_dir, 'kingfisher_download', '{sample}.log')),
    benchmark: os.path.abspath(os.path.join(logs_dir, 'kingfisher_download', '{sample}.benchmark')),
    shell:
        """
        cd {params.kf_dir}
        kingfisher get --download-threads {threads} -t {threads} \
                       -r {wildcards.sample} -m ena-ftp aws-http 2> {log}
        cd -
        """