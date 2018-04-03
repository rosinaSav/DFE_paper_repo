'''
Author: Rosina Savisaar.
Prepare input file for running MultiDFEest.
'''

from housekeeping import list_to_dict, parse_arguments, print_elements, run_process
from INSIGHT import shuffle_dictionaries
import read_and_write as rw

def get_SFS(positions, SNPs, to_remove, N):
    '''
    Given two dictionaries with transcripts IDs as keys and one with hit positions as values
    and another with subdictionaries as values (containing SNP positions as keys along with
    their minor allele counts as values), generate a site frequency spectrum.
    to_remove contains positions where SNP positions were filtered out.
    These must be ignored when making the SFS.
    '''
    N = int(N)
    #note that for multiDFEest, the SFS must have a 0 category
    SFS = {i: 0 for i in range(N + 1)}
    for trans in positions:
        if trans in to_remove:
            to_remove[trans] = [i for i in to_remove[trans] if i not in SNPs[trans]]
            for position in positions[trans]:
                if position not in to_remove[trans]:
                    if position not in SNPs[trans]:
                        #update 0 category
                        SFS[0] = SFS[0] + 1
                    else:
                        SFS[SNPs[trans][position]] = SFS[SNPs[trans][position]] + 1
        else:
            SFS[0] = SFS[0] + len(positions[trans])
    SFS = [SFS[i] for i in sorted(SFS.keys())]
    return(SFS)

def parse_pos(file):
    '''
    Parse a hits/controls positions file.
    '''
    pos_list = rw.read_many_fields(file, "\t")
    pos_dict = list_to_dict(pos_list, 0, 1)
    pos_dict = {i: [int(j) for j in pos_dict[i].split(",") if j != ""] for i in pos_dict}
    return(pos_dict)

def main():
    description = "Prepare input file for running MultiDFEest."
    args = parse_arguments(description, ["hit_file", "control_file", "SNPs_file_prefix", "N", "output_file", "per_chrom_files", "shuffle"], ints = [3], flags = [5, 6])
    hit_file, control_file, SNPs_file_prefix, N, output_file, per_chrom_files, shuffle = args.hit_file, args.control_file, args.SNPs_file_prefix, args.N, args.output_file, args.per_chrom_files, args.shuffle

    hits = parse_pos(hit_file)
    controls = parse_pos(control_file)

    if shuffle:
        hits, controls = shuffle_dictionaries(hits, controls)

    SNPs = {}
    to_remove_all = {}
    #if the data is stored chromosome by chromosome, rather than all combined
    if per_chrom_files:
        for chrom in range(1, 23):
            try:
                SNPs_file = "{0}{1}.bed".format(SNPs_file_prefix, str(chrom))
                current_SNPs = rw.read_many_fields(SNPs_file, "\t")
                to_remove = list_to_dict(current_SNPs, 0, 2)
                to_remove = {i: to_remove[i].split(",") for i in to_remove}
                current_SNPs = list_to_dict(current_SNPs, 0, 1)
                for trans in current_SNPs:
                    if trans in controls:
                        SNPs[trans] = {}
                        trans_SNPs = current_SNPs[trans]
                        if trans_SNPs:
                            trans_SNPs = [i.split(",") for i in trans_SNPs.split("|")]
                            #this is where you get the allele count
                            trans_SNPs = list_to_dict(trans_SNPs, 0, 3)
                            trans_SNPs = {int(i): int(trans_SNPs[i]) for i in trans_SNPs}
                            SNPs[trans] = trans_SNPs
                        to_remove_all[trans] = [int(i) for i in to_remove[trans] if i not in ["error", ""]]
            except FileNotFoundError:
                pass
    else:
        SNPs_file = SNPs_file_prefix
        current_SNPs = rw.read_many_fields(SNPs_file, "\t")
        to_remove = list_to_dict(current_SNPs, 0, 2)
        to_remove = {i: to_remove[i].split(",") for i in to_remove}
        current_SNPs = list_to_dict(current_SNPs, 0, 1)
        counter = 0
        for trans in current_SNPs:
            if trans in controls:
                SNPs[trans] = {}
                trans_SNPs = current_SNPs[trans]
                if trans_SNPs:
                    trans_SNPs = [i.split(",") for i in trans_SNPs.split("|")]
                    #this is where you get the allele count
                    trans_SNPs = list_to_dict(trans_SNPs, 0, 3)
                    trans_SNPs = {int(i): int(trans_SNPs[i]) for i in trans_SNPs}
                    SNPs[trans] = trans_SNPs
                to_remove_all[trans] = [int(i) for i in to_remove[trans] if i not in ["error", ""]]

    hit_SFS = get_SFS(hits, SNPs, to_remove_all, N)
    control_SFS = get_SFS(controls, SNPs, to_remove_all, N)

    with open(output_file, "w") as file:
        file.write("{0}\n".format(N))
        file.write(" ".join([str(i) for i in hit_SFS]))
        file.write("\n")
        file.write(" ".join([str(i) for i in control_SFS]))
        file.write("\n")

if __name__ == "__main__":
    main()
