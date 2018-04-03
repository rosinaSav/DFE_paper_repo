'''
Author: Rosina Savisaar.
Extract exon end/core regions from feature set.
'''

from bedtools_games import setup
from housekeeping import parse_arguments

def main():
    
    description = "Extract exon end/core regions from feature set."
    args = parse_arguments(description, ["genome", "features_file", "families_file", "dataset", "start_only"], flags = [4])
    genome, features_file, families_file, dataset, start_only = args.genome, args.features_file, args.families_file, args.dataset, args.start_only

    genome = "hg38"
    features_file = "general/Homo_sapiens.GRCh38.85.gtf"
    families_file = "general/filtered_hg38_85_pc_multiexon_families.txt"
    dataset = "filtered_hg38_85_pc_multiexon"

    #prepare feature set, get necessary genomic features
    fs = setup(features_file, genome, dataset, families_file = families_file)
    exons = fs.get_exons()
    CDS = fs.get_CDS()

    #pick a random member from each family
    picked = fs.pick_random_members()
    exons = {i: exons[i] for i in picked}
    CDS = {i: CDS[i] for i in picked}

    if start_only:
        #only the 5' end
        fs.get_exon_beginnings(exons, CDS, file_prefix = "general/{0}".format(dataset), write_to_fasta = True)
    else:
        #both flanks and the core
        fs.get_exon_cores_and_flanks(exons, CDS, file_prefix = "general/{0}".format(dataset), write_to_fasta = True)

if __name__ == "__main__":
    main()
