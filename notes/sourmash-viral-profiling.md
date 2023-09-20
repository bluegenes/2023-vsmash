# Sourmash for viral taxonomic profiling/classification

[![hackmd-github-sync-badge](https://hackmd.io/-5V2nkVyRgObc_dFkWBeyA/badge)](https://hackmd.io/-5V2nkVyRgObc_dFkWBeyA)

> a preprint (below) uses sourmash through WhatThePhage and claims it performs poorly for viral classification. The workflow conducts contig-level classification using `k21,scaled100` and `sourmash search` to the phage database using `jaccard` (not `containment`)

their commands:

- database preparation:
    ```
    sourmash compute --scaled 100 -k 21 --singleton \
    --seed 42 -p 8 -o phages.sig ${references}
    
    sourmash index phages.sbt.zip phages.sig
    ```
- classification:
    ```
    sourmash compute -p ${task.cpus} \
     --scaled 100 -k 21 \${fastafile}

    sourmash search -k 21 \
    ${signature} phages.sbt.zip -o \
    ${signature}.temporary
    ```

> Notes:
> - code context: [WtP sourmash taxonomic classification](https://github.com/replikation/What_the_Phage/blob/18b39e060edf0001a2d0dfc07748005681cc0c00/workflows/process/phage_tax_classification/sourmash_for_tax.nf#L4)
> - [What the Phage manuscript](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9673492/) describes using sourmash 2.0.1 - is sourmash pinned to 2.x in WtP?

**Goal: Assess taxonomic profiling (+classification?) using `gather`--> `tax` workflow workflow on mock and real datasets**

## Benchmarking Reference
 - [Benchmarking Bioinformatic Virus Identification Tools Using Real-World Metagenomic Data across Biomes](https://www.biorxiv.org/content/10.1101/2023.04.26.538077v2)
     > As most viruses remain uncultivated, metagenomics is currently the main method for virus discovery. Detecting viruses in metagenomic data is not trivial. In the past few years, many bioinformatic virus identification tools have been developed for this task, making it challenging to choose the right tools, parameters, and cutoffs. As all these tools  measure  different  biological  signals,  and  use  different  algorithms  and  training/reference  databases,  it  is imperative  to  conduct  an  independent  benchmarking  to  give  users  objective  guidance.  We  compared  the performance of ten state-of-the-art virus identification tools in thirteen modes on eight paired viral and microbial datasets from three distinct biomes, including a new complex dataset from Antarctic coastal waters. The tools had highly variable true positive rates (0 – 68%) and false positive rates (0 – 15%). PPR-Meta best distinguished viral  from  microbial  contigs,  followed  by  DeepVirFinder,  VirSorter2,  and  VIBRANT.  Different  tools  identified different subsets of the benchmarking data and **all tools, except for Sourmash, found unique viral contigs**. Tools performance could be improved with adjusted parameter cutoffs, indicating that adjustment of parameter cutoffs before usage should be considered. Together, our independent benchmarking provides guidance on choices of bioinformatic virus identification tools and gives suggestions for parameter adjustments for viromics researchers.

## Test Datasets

### Mock/synthetic datasets:
- `SRR3458562`-`SRR3458569`
- used for benchmarking ["Genome Detective"](https://academic.oup.com/bioinformatics/article/35/5/871/5075035) software
    > We first validated Genome Detective using a synthetic virus dataset (NCBI SRA: SRR3458562-SRR3458569), originally prepared to optimize laboratory-based virus extraction procedures, in which viruses were carefully selected to cover the range of naturally occurring diversity (Conceição-Neto et al., 2015). This published dataset also includes carefully validated quantitative results, confirmed with quantitative PCR. Genome Detective identified all of the viruses in the synthetic dataset. We then validated Genome Detective with real clinical datasets. In total, we analyzed 208 datasets, which are available via Sequence Read Archive (SRA) or the European Nucleotide Archive. We then compared our results to the published results and found a >95% concordance, successfully identifying 257 viral species (Table 1 and Supplementary Table). These included single viruses with unsegmented (HIV) and segmented genomes (Influenza A, Rotavirus, MERS) from amplicon-based NGS sequenced as well as unbiased metagenomic datasets (Table 1). Overall, precision, sensitivity and specificity were high, with the exception of 20 metagenomic datasets from human fecal (ERR233412-ERR233431), which had scarce viral reads (Supplementary Table).
- Reference for the design of these mock datasets: [Modular approach to customise sample preparation procedures for viral metagenomics: a reproducible protocol for virome analysis](https://www.nature.com/articles/srep16532)

### Real metagenome/virome datasets

use biome datasets from: [Benchmarking Bioinformatic Virus Identification Tools Using Real-World Metagenomic Data across Biomes](https://www.biorxiv.org/content/10.1101/2023.04.26.538077v2)
> We collected a total of 48 metagenome datasets, including eight paired viral and microbial datasets from each biome (Supplementary Table S4).

> We used unique contigs with lengths of at least 1,500 bp from the viral and microbial size fractions as ground truth positives and negatives, respectively. Viral contigs that were identified as viral and non-viral by the tools were regarded as true positives and false negatives, respectively. Microbial contigs that were identified as viral and non-viral were regarded as false positives and true negatives, respectively (Figure 1).

- full github associated with the paper: https://github.com/MGXlab/virus_identification_tools_benchmarking


## Approach

- apply gather-> tax workflow with appropriate parameters
- use mock datasets, real datasets --> true positives, false positives, F1 score, etc. Follow benchmarking preprint & Portik et al methods for assessing datasets.

### 1. Sourmash profiling /classification workflow:
- gather--> tax approach on whole dataset, not contig-level (profiling)
    - for "taxonomic classification" (reads or contigs), we could add a contig/read level gather --> tax step. 
        - use the profiling results as a smaller database, then 
        - assign reads/contigs at high resolution (sc 10?)
    - needs `singleton` `manysearch` - should work with pyo3-branchwater workflow, i think? Though we can't currently `manysketch` with `--singleton`...

### 2. Download, sketch reference dbs (DNA, Protein?)
- What the Phage
- Human-associated viral dbs
    - Gut phage database (human)
    - gut virome database (human)
- Refseq / genbank
- PIGEON (or part of PIGEON)
- any additional used in benchmarking paper

### 3. Assess performance
- develop notebook for estimating TP, FP, F1score, etc for mock datasets. Plots!
- optional, NTP: turn this into plugin so it can be used for more tax classification assessments

### 4. Sweep sketch/search params
- dna
    - k=8-15,21,31 (annie suggests k=8)
    - scaled 50-200?
- protein
    - k=4-15? (ntp: try k=7)
    - scaled?? (10-200)

## 5. Products
- paper/tutorial/blog post?
- submit PR with updated workflow to What the Phage?
- submit comment to benchmarking preprint?

---

## To Do:

- tessa: snakefile + download, sketch, gather
    - [ ] Get the reference dbs on FARM
    - [ ] Snakefile: sketch dbs across ksizes, scaled
    - [ ] Snakefile: gather virome --> all ref dbs
- annie: find info on community composition of mock datasets from https://academic.oup.com/bioinformatics/article/35/5/871/5075035 / https://www.nature.com/articles/srep16532

---

## Done

- [x] Mock communities and real virome data? ID some representative datasets.  
    - [x] We need to create some sort of mock virome community, where we know what reads are in there (which viruses), and see if we can retrieve those kmers from a reference db. I don't know exactly how to make a mock virome. or a 'metagenome' of which we know exactly what sequences are in there. 
**Ideas?**