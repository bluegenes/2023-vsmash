#! /usr/bin/env python
import sys
import gzip
import os
import csv
import argparse

from Bio import Entrez, SeqIO
from Bio.Seq import UndefinedSequenceError

# Set your email address (required by NCBI)
Entrez.email = "ntpierce@ucdavis.edu"

def main(args):
    """
    Downloads nucleotide and/or protein sequences for the given GenBank accession number
    and saves them as gzip-compressed FASTA files with the specified filenames, or with
    default filenames based on the accession number.

    Args:
        accession (str): The GenBank accession number to download.
        nucleotide_file (str): The filename for the nucleotide FASTA file, or None to use the default filename.
        protein_file (str): The filename for the protein FASTA file, or None to use the default filename.

    Returns:
        None
    """
    accession = args.accession

    # Use Entrez to download the record
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")

    # Use SeqIO to parse the record and extract the organism name and sequences
    record = SeqIO.read(handle, "genbank")
    organism_name = record.annotations["organism"]
    # nucl
    nucleotide_sequence=None
    try:
        seq = str(record.seq)
        nucleotide_sequence = seq
    except UndefinedSequenceError:
        sys.stderr.write(f"Could not find sequence for accession {accession}: record.seq is undefined.")
    # prot
    protein_sequences = []
    for feature in record.features:
        if feature.type == "CDS":
            if "translation" in feature.qualifiers:
                #import pdb;pdb.set_trace()
                #protein_sequences.append((feature.qualifiers["translation"][0], feature.qualifiers["locus_tag"][0]))
                protein_sequences.append((feature.qualifiers["translation"][0], feature.qualifiers["protein_id"][0]))

    # Determine the output filenames
    nucleotide_file = args.nucleotide
    protein_file = args.protein
    fileinfo_file = args.fileinfo
    protein_exists=False
    if args.nucleotide is None and nucleotide_sequence is not None:
        nucleotide_file = f"{accession}.fna.gz"
    if args.protein is None and len(protein_sequences) > 0:
        protein_file = f"{accession}.faa.gz"
    if args.fileinfo is None:
        fileinfo_file = f"{accession}.fileinfo.csv"

    # Save the sequences to gzip-compressed FASTA files
    if nucleotide_sequence is not None:
        with gzip.open(nucleotide_file, "wt") as f:
            f.write(f">{accession}\n{nucleotide_sequence}\n")
    if len(protein_sequences) > 0:
        protein_exists = True
    if protein_exists:
        with gzip.open(protein_file, "wt") as f:
            #for i, (protein_sequence, locus_tag) in enumerate(protein_sequences):
            #    f.write(f">{accession}_CDS{i+1}|{locus_tag}\n{protein_sequence}\n")
            for i, (protein_sequence, protein_id) in enumerate(protein_sequences):
                f.write(f">{accession}_CDS{i+1}|{protein_id}\n{protein_sequence}\n")

    # Output file details on to a CSV file
    name = accession + ' ' + organism_name
    sys.stderr.write(f'name: {name}')
    with open(fileinfo_file, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        if protein_exists and nucleotide_sequence is not None:
            writer.writerow([name, nucleotide_file, protein_file])
        elif nucleotide_sequence is not None:
            writer.writerow([name, nucleotide_file,""])
        else:
            writer.writerow([name, "" ,""])



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    parser = argparse.ArgumentParser(description="Download nucleotide and/or protein sequences for a GenBank accession")
    parser.add_argument("accession", type=str, help="The GenBank accession number to download")
    parser.add_argument("--nucleotide", type=str, help="The filename for the nucleotide FASTA file")
    parser.add_argument("--protein", type=str, help="The filename for the protein FASTA file")
    parser.add_argument("--fileinfo", type=str, help="Filename details for downstream. Only includes protein file if it existed for download.")

    # Parse command line arguments
    args = parser.parse_args()

    # Download sequences for the specified accession
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
