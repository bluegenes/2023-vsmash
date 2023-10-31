import sys, os
import argparse
import csv
import screed
import time

from Bio import Entrez, SeqIO

def get_accessions_from_fasta(fasta):
    """
    Parse fasta file with screed to get accession numbers
    :param fasta: fasta file
    :return: list of accession numbers
    """
    # format gi|2497518365|ref|NC_074657.1| Nitrosopumilus spindle-shaped virus isolate NSV2, partial genome
    accInfo = []
    for record in screed.open(fasta):
        # extract accession number, e.g. NC_074657.1 from gi|2497518365|ref|NC_074657.1| Nitrosopumilus spindle-shaped virus isolate NSV2, partial genome
        _, _, _, accession,name = record.name.split('|')
        accInfo.append([accession, name])
    return accInfo


def main(args):

    Entrez.email = args.email
    # open input fasta and parse with screed to get accession numbers
    info = []
    for fasta in args.input:
        info.extend(get_accessions_from_fasta(fasta))
    n_acc = len(info)
    print(f"Found {n_acc} accessions")
    # number of seconds to sleep for between requests
    sleep_seconds = 15
    # open output file and use csvwriter
    with open(args.output, 'w') as out:
        taxranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'clade']
        header = ['accession', 'name', 'taxid', 'scientific_name'] + taxranks
        writer = csv.DictWriter(out, delimiter=',', fieldnames=header)
        # for each accession number, get taxonomic info from NCBI
        reporting_threshold = max(round(n_acc / 20), 1) # 5 percent
        for n, (acc,name) in enumerate(info):
            # report every 5 percent
            if n % reporting_threshold == 0:
                print(f"Getting taxonomic info for {n} of {n_acc} accessions ({round(n/n_acc*100)}%)")
            # get taxonomic info from NCBI
            genbank_record = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(genbank_record, "genbank")
            # to get taxid, search with organism name (bc scientific name not always properly populated)
            organism = record.annotations['organism'] 
            time.sleep(sleep_seconds)  # Sleep to avoid too many requests
            handle = Entrez.esearch(db="taxonomy", term=organism)
            record = Entrez.read(handle)
            taxid = record['IdList'][0]

            # Add a delay before fetching taxonomy information
            time.sleep(sleep_seconds)  # Sleep to avoid too many requests

            # Use taxid to fetch the taxonomic information with ranks
            try:
                taxonomy_summary = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                taxonomy_record = Entrez.read(taxonomy_summary)
                # Extract the scientific name and taxonomic lineage with ranks
                scientific_name = taxonomy_record[0]["ScientificName"]
            except:
                print(f"Error getting taxonomic info for accession {acc} with taxid {taxid}")
                sys.exit(1)
            
            accInfo = {"accession": acc,
                       "name": name.strip(),
                       "taxid": taxid,
                       "scientific_name": scientific_name}
            lineage_ex = []
            lineage_ex = taxonomy_record[0].get("LineageEx")
            if not lineage_ex:
                lineage = taxonomy_record[0].get("Lineage")
                print(f"LineageEx not found for accession {acc} with taxid {taxid}")
            for taxon in lineage_ex:
                rank = taxon["Rank"]
                tname = taxon["ScientificName"]
                accInfo[rank] = tname
            with open(args.output, 'a') as out:
                try:
                    writer.writerow(accInfo)
                except ValueError as e:
                    print(f"Error writing accession {acc} to file {args.output}")
                    print(e)
                    print(accInfo)
                    sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get taxonomic info from NCBI for genome fasta(s).')
    parser.add_argument('-i', '--input', required=True, nargs='+', help='Input fasta files')
    parser.add_argument('-o', '--output', required=True, help='csv file with taxonomic info')
    parser.add_argument('--email', help='Email address for NCBI query', default= "ntpierce@ucdavis.edu")
    args = parser.parse_args()

    main(args)