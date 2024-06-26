import re
import os
import pandas as pd

refseq69_roux = 'List_phages_RefSeq_69.txt'
out_dir = 'output.refseq69'
logs_dir = os.path.join(out_dir, 'logs')
basename = 'refseq69_phages'

# regex pattern to match the ref| identifier and following text (name)
# acc_pattern = r'ref\|([^|]+)\|' # just ident
pattern = r'ref\|([^|]+)\| (.+)$' # ident + name
# Note: the name text has qualifiers -- just leave in name
# - ', complete genome' ends 1455 of 1532 entries
# - others have : ', complete sequence', 'genomic sequence', or nothing


ACCESSIONS = {}

with open(refseq69_roux, 'r') as inF:
    for line in inF:
        match = re.search(pattern, line)
        if match:
            # Extract the matched identifier
            ident = match.group(1)
            # extract name
            name = match.group(2).replace(',', ';')
            ACCESSIONS[ident] = ident + ' ' + name
        else:
            print(f'No match found in line: {line}')

rule all:
    input:
        os.path.join(out_dir, f"{basename}.taxonomy.csv"),
        expand(f"{out_dir}/{{db_basename}}.{{moltype}}.zip", db_basename=basename, moltype=['dna', 'protein']),
        expand(f"{out_dir}/{{db_basename}}.{{moltype}}-sc{{scaled}}.zip", db_basename=basename, moltype=['dna', 'protein'], scaled=[2, 10, 50, 100, 1000]),

class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern
    
    def get_filenames(self, basename=None, moltype=None):
        df = pd.read_csv(f"{out_dir}/{basename}.fromfile.csv")
        filename_col = 'genome_filename'
        if moltype == "protein":
            filename_col = 'protein_filename'
        # filter df to get non-empties in relevant column
        fastas = df[filename_col][df[filename_col].notnull()].tolist()
        return fastas

    def __call__(self, w):
        global checkpoints
        # wait for the results of 'check_fromfile'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_fromfile.get(**w)
        #expand the pattern
        fastas = self.get_filenames(**w)

        pattern = expand(self.pattern, fn =fastas,**w)
        return pattern


# download FASTA and collect info
rule download_by_accession:
    output: 
        nucl=protected(os.path.join(out_dir, "genomic/{acc}.fna.gz")),
        fileinfo=protected(os.path.join(out_dir, "fileinfo/{acc}.fileinfo.csv")),
    params:
        prot=protected(os.path.join(out_dir, "protein/{acc}.faa.gz")),
        prot_out_dir=os.path.join(out_dir, 'protein') 
    conda: "conf/env/biopython.yml"
    log: os.path.join(logs_dir, "downloads", "{acc}.log")
    benchmark: os.path.join(logs_dir, "downloads", "{acc}.benchmark")
    threads: 1
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    shell:
        """
        mkdir -p {params.prot_out_dir} # since protein output is in params, need to manually make sure this dir is created
        python genbank_nuccore.py {wildcards.acc} --nucleotide {output.nucl} --protein {params.prot} --fileinfo {output.fileinfo} 2> {log}
        """

rule translate_protein:
    input:
        fileinfo=os.path.join(out_dir, "fileinfo/{acc}.fileinfo.csv"),
        nucl=os.path.join(out_dir, "genomic/{acc}.fna.gz"),
    output:
        translated=protected(os.path.join(out_dir, "translate/{acc}.faa.gz")),
    conda: "conf/env/seqkit.yml"
    log: os.path.join(logs_dir, 'seqkit', '{acc}.log')
    threads: 1
    shell:
        """
        seqkit translate {input.nucl} | gzip > {output} 2> {log}
        """

rule aggregate_fileinfo_to_fromfile:
    input: 
        fileinfo=expand(os.path.join(out_dir, "fileinfo/{acc}.fileinfo.csv"), acc=ACCESSIONS)
    output:
        csv = protected(os.path.join(out_dir, "{basename}.fromfile.csv"))
    run:
        with open(str(output.csv), "w") as outF:
            header = 'name,genome_filename,protein_filename'
            outF.write(header + '\n')
            for inp in input:
                with open(str(inp)) as inF:
                    # outF.write(inF.read())
                    # if protein_filename doesn't exist, use translated
                    line = inF.read().split(',')
                    if len(line) != 4:
                        import pdb; pdb.set_trace()
                    # name,dna,prot,_tax = inF.read().split(',')
                    name,dna,prot,_tax = line
                    this_acc = name.split(' ')[0]
                    if this_acc not in prot:
                        prot = f"{out_dir}/translate/{this_acc}.faa.gz"
                    outF.write(f"{name},{dna},{prot}\n")


rule aggregate_fileinfo_to_taxonomy:
    input: 
        fileinfo=expand(os.path.join(out_dir, "fileinfo/{acc}.fileinfo.csv"), acc=ACCESSIONS)
    output:
        csv = protected(os.path.join(out_dir, "{basename}.taxonomy.csv"))
    run:
        with open(str(output.csv), "w") as outF:
            header = 'ident,superkingdom,phylum,class,order,family,genus,species,name'
            outF.write(header + '\n')
            for inp in input:
                with open(str(inp)) as inF:
                    name,_dna,_prot,tax = inF.read().split(',')
                    tax = tax.replace(' ', '').replace(';', ',').strip() # remove spaces and replace ';' with ',', strip '\n'
                    while len(tax.split(',')) < 7:
                        tax += ','
                    this_acc = name.split(' ')[0]
                    outF.write(f"{this_acc},{tax},{name}\n")


# Define the checkpoint function that allows us to read the fromfile.csv
checkpoint check_fromfile:
    input: os.path.join(out_dir, f"{basename}.fromfile.csv"),
    output: touch(os.path.join(out_dir,".check_fromfile"))

paramD = {"dna": "dna,k=9,k=11,k=13,k=15,k=17,k=19,k=21,k=31,scaled=1,abund", "protein": "protein,k=5,k=6,k=7,k=8,k=9,k=10,scaled=1,abund"}
rule sketch_fromfile:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
    output: os.path.join(out_dir, "{basename}.{moltype}.zip")
    params:
        lambda w: paramD[w.moltype]
    threads: 10
    resources:
        mem_mb=3000,
        runtime=60,
        time=90,
        partition="low2",
    # conda: "conf/env/sourmash.yml"
    conda: "conf/env/branchwater.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.benchmark")
    wildcard_constraints:
        moltype="dna|protein",
    shell:
        """
        sourmash scripts manysketch {input.fromfile} -p {params} \
                                    -o {output} -c {threads} 2> {log}
        """


rule downsample_zip:
    input: os.path.join(out_dir, "{basename}.{moltype}.zip")
    output: os.path.join(out_dir, "{basename}.{moltype}-sc{scaled}.zip")
    threads: 1
    conda: "conf/env/branchwater.yml"
    log: os.path.join(logs_dir, "downsample", "{basename}.{moltype}.sc{scaled}.log")
    benchmark: os.path.join(logs_dir, "downsample", "{basename}.{moltype}.sc{scaled}.benchmark")
    params:
        alpha_cmd=lambda w: f"--{w.moltype}" if w.moltype == "protein" else "",
    shell:
        """
        sourmash sig downsample {input} --scaled {wildcards.scaled} {params.alpha_cmd} -o {output} 2> {log}
        """
