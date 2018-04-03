'''
Author: Rosina Savisaar.
Run INSIGHT on a set of sequences and a set of sites.
'''

from bedtools_games import bed_to_CDS_indices, convert_coords, Feature_Set, map_relative_to_feature, write_features_to_bed
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
import conservation
import copy
import ensembl_ops as eo
from housekeeping import flatten, list_to_dict, make_dir, parse_arguments, pause_script, print_elements, remove_file, run_in_parallel, run_process, update_counter
import nucleotide_comp as nc
import numpy as np
import operator
import os
import random
import re
import read_and_write as rw
import scipy.stats

#globals
#this is to have the 4 bases in alphabetical order, which isn't normally the case
base_order = sorted(nc._canon_bases_)
lambda_regex = re.compile("(?<=\(homo:)[\d\.\-e]*(?=\,pan)")

def filter_on_SNPs(current_CDSs, current_SNPs, current_SNPs_to_remove, combined_dict, pos_dict, hit_dict):
    '''
    Filter out hit/control positions that overlap with SNP positions that were filtered out
    (because you don't want a SNP that was filtered out to count as a monomorphic site).
    '''
    for trans in current_CDSs:
        #get current SNPs
        if trans in current_SNPs:
            trans_SNPs = list(current_SNPs[trans].keys())
        else:
            trans_SNPs = []
        temp_pos = []
        temp_hits = []
        #put everything that doesn't overlap a filtered out SNP into a temp list
        for position in combined_dict[trans]:
            if (trans not in current_SNPs_to_remove) or (position not in current_SNPs_to_remove[trans]) or (position in trans_SNPs):
                if position in pos_dict[trans]:
                    temp_pos.append(position)
                else:
                    temp_hits.append(position)
        #controls
        pos_dict[trans] = sorted(temp_pos)
        #hits
        hit_dict[trans] = sorted(temp_hits)
        #the two combined
        combined_dict[trans] = sorted(temp_pos + temp_hits)
    return(pos_dict, hit_dict, combined_dict)

def get_lambda(lambda_file_outroot, phy_file, subst_model, min_inf = None):
    '''
    Calculate lambda input parameter for INSIGHT.
    '''
    lambda_file = "{0}.mod".format(lambda_file_outroot)
    #to make sure you catch it if the phyloFit process fails
    remove_file(lambda_file)
    #from UCSC
    tree_file = "DFE/UCSC_model.mod"
    #subst_model is JC69, for instance
    #scale-only, cause you don't want it to estimate a new tree, just to scale the whole thing
    arguments = ["phyloFit", "--init-model", tree_file, "--out-root", lambda_file_outroot, "--subst-mod", subst_model,
                           "--msa-format", "PHYLIP", "--scale-only", phy_file]
    #must be set to False for testing
    if min_inf:
        arguments.extend(["-I", min_inf])
    results = run_process(arguments)
    with open(lambda_file) as file:
        lambda_b = file.read()
    lambda_b = re.findall(lambda_regex, lambda_b)[0]
    return(lambda_b)

def get_MSA(chroms, chrom_dict, control_file, hit_file, CDSs, lengths, names, seqs, clean_names, freq_threshold, dataset, suffix, genome, output_folder, general_folder, n, SNP_file_name_prefix, CDS_SNP_file_name_prefix, MSA_file_name_prefix, new_SNPs, new_MSA, shuffle, remove_GT, big_tree, degen_hits = None, degen_controls = None, hit_reduce = 0, control_reduce = 0):
    '''
    Get SNP data and MSA, and prepare input files for INSIGHT based on the data retrieved.
    '''
    subst_model = "JC69"
    #n is population size, a is a parameter used for calculating theta for INSIGHT
    a = np.sum([(1/i) for i in range(1, n)])

    #so you wouldn't have to reprocess the SNPs if it's already been done
    if SNP_file_name_prefix == "None":
        SNP_file_name_prefix = "{0}/{1}_SNPs_tabix".format(general_folder, dataset)
    if CDS_SNP_file_name_prefix == "None":
        CDS_SNP_file_name_prefix = "{0}/{1}_SNPs_relative_more_info".format(general_folder, dataset)
        
    neutral_output = []
    hit_output = []
    chroms_to_keep = []
    hit_counts = {}
    control_counts = {}
    for chrom in sorted(chroms):

        #species in MSA
        phylip_data = {"homo_sapiens": [], "pongo_abelii": [], "macaca_mulatta": [], "pan_troglodytes": []}
        if big_tree:
            phylip_data = {"gorilla_gorilla": [], "callithrix_jacchus": [], "papio_anubis": [], "chlorocebus_sabaeus": [], "homo_sapiens": [], "pongo_abelii": [], "macaca_mulatta": [], "pan_troglodytes": []}

        print(chrom)

        #read in hits/controls for current chromosome
        hit_dict = parse_pos_file(chrom_dict, chrom, hit_file)
        pos_dict = parse_pos_file(chrom_dict, chrom, control_file)

        #if running a negative control
        if shuffle:
            print("Shuffling dictionaries...")
            hit_dict, pos_dict = shuffle_dictionaries(hit_dict, pos_dict)

        #reduce the sizes of the hit/control lists by a specified fraction
        if hit_reduce > 0:
            hit_dict = reduce_dict(hit_dict, hit_reduce)
            pos_dict = reduce_dict(pos_dict, control_reduce)

        #total number of neutral sites
        total_neutral = np.sum([len(pos_dict[i]) for i in pos_dict if chrom_dict[i] == chrom])

        print("Getting SNPs...")
        current_CDSs = {i: CDSs[i] for i in CDSs if CDSs[i][0][0][0] == chrom}
        if current_CDSs and total_neutral >= 150:
            #this is to check whether the file you've provided has all the SNPs or whether there's a separate file for each chromosome
            if SNP_file_name_prefix[-4] == ".":
                SNP_file_name = SNP_file_name_prefix
            else:
                SNP_file_name = "{0}_{1}.bed".format(SNP_file_name_prefix, chrom)
            if CDS_SNP_file_name_prefix[-4] == ".":
                CDS_SNP_file_name = CDS_SNP_file_name_prefix
            else:
                CDS_SNP_file_name = "{0}_{1}.bed".format(CDS_SNP_file_name_prefix, chrom)
            
            #only get new SNPs if the file isn't already there
            if os.path.isfile(SNP_file_name):
                current_new_SNPs = False
            else:
                current_new_SNPs = True
            if os.path.isfile(CDS_SNP_file_name):
                current_new_parse = False
            else:
                current_new_parse = True
            #get SNP positions in relative CDS coordinates
            current_SNPs, current_SNPs_to_remove = conservation.get_relative_SNPs(current_CDSs, SNP_file_name, CDS_SNP_file_name, seqs, names, genome, get_new_SNPs = current_new_SNPs, parse_SNPs = current_new_parse, remove_GT = remove_GT)
            #remove positions that overlap the stop codon
            pos_dict = {i: [j for j in pos_dict[i] if (lengths[i] - j) > 3] for i in pos_dict}
            hit_dict = {i: [j for j in hit_dict[i] if (lengths[i] - j) > 3] for i in hit_dict}
            combined_dict = {i: pos_dict[i] + hit_dict[i] for i in hit_dict}

            #filter the hits/controls to remove sites that overlap SNP positions that were filtered out, so they
            #wouldn't count as monomorphic
            pos_dict, hit_dict, combined_dict = filter_on_SNPs({i: current_CDSs[i] for i in current_CDSs if i in hit_dict}, current_SNPs, current_SNPs_to_remove, combined_dict, pos_dict, hit_dict)

            #get MSA from Ensembl Compara database
            coords_file = "temp_data/coords_file{0}.txt".format(random.random())
            MSA_file = "{0}_{1}.txt".format(MSA_file_name_prefix, chrom)
            if new_MSA:
                if not os.path.isfile(MSA_file):
                    print("Obtaining MSA...")
                    eo.get_MSA_gene_list(current_CDSs, coords_file, "EPO", "primates", 85, "homo_sapiens", MSA_file)
                    os.remove(coords_file)
                    eo.flush_tables("localhost", "mysql", "fackel")
            MSA_raw = eo.parse_MSA_output(MSA_file)
            phy_file = "{0}/{1}_{2}_phylip{3}.phy".format(output_folder, dataset, suffix, chrom)
            phy_file_neutral = "{0}/{1}_{2}_phylip{3}_neutral.phy".format(output_folder, dataset, suffix, chrom)
            #if you're ignoring degenerate sites
            if degen_hits:
                degen = {}
                for trans in degen_hits:
                    degen[trans] = degen_hits[trans].copy()
                    degen[trans].update(degen_controls[trans])
            else:
                degen = None
            #based on the MSA and the hit/control positions, create a phylip alignment suitable for phyloFit
            current_CDSs, separate_to_concat_mapping, current_data, current_data_neutral = positions_to_phylip(MSA_raw, CDSs, phylip_data, combined_dict, pos_dict, lengths, phy_file, phy_file_neutral, clean_names, current_SNPs, degen = degen)
    ##        os.remove(MSA_file)
            current_neutral_length = len(current_data_neutral[0])
            #filter out chromosomes where there isn't enough information
            if current_neutral_length >= 150:
                chroms_to_keep.append(chrom)
                #get the tuple for each position in the alignment (a string where each characater corresponds to the base
                #present in a particular species in the alignment)
                tuples_mapping = get_tuples_mapping(current_data)

                print("Calculating prior probabilities...")
                outroot = "{0}/{1}_{2}_phylip{3}".format(output_folder, dataset, suffix, chrom)
                neutral_outroot = "{0}_neutral".format(outroot)

                #get lambda from control sites (scaling factor on tree)
                lambda_b = get_lambda(neutral_outroot, phy_file_neutral, "JC69")

                #get prior probabilities for each of the bases at each site
                model_file = "{0}.mod".format(neutral_outroot)
                pp_final = get_pp(outroot, subst_model, phy_file, model_file, separate_to_concat_mapping, combined_dict, tuples_mapping)

                #format output data appropriately for INSIGHT
                print("Storing information for writing to file later on...")
                control_counts[chrom] = np.sum([len(pos_dict[i]) for i in pos_dict])
                hit_counts[chrom] = np.sum([len(hit_dict[i]) for i in hit_dict])
                current_hit_output, current_neutral_output = get_output_string(pos_dict, hit_dict, pp_final, current_CDSs, current_SNPs, seqs, names, chrom, a, lambda_b, freq_threshold, degen = degen)
                hit_output.append(current_hit_output)
                neutral_output.append(current_neutral_output)
    return(hit_output, neutral_output, chroms_to_keep, hit_counts, control_counts)

def get_output_string(pos_dict, hit_dict, pp_dict, current_CDSs, current_SNPs, seqs, names, chrom, a, lambda_b, freq_threshold, degen = None):
    '''
    Format input string for INSIGHT.
    '''
    current_neutral_output = []
    current_hit_output = []
    #filter out hits/controls for which there is no information on prior probabilities
    combined_dict = {i: sorted(pos_dict[i] + hit_dict[i]) for i in pp_dict}
    #convert relative CDS coordinates to absolute chrom. coordinates
    mapping_dict = map_relative_to_feature(combined_dict, current_CDSs)
    K = 0
    site_count = 0
    for trans in sorted(pp_dict.keys()):
        #get SNPs for current transcript
        if trans in current_SNPs:
            trans_SNPs = current_SNPs[trans]
            #filter out degenerate SNPs
            if degen:
                trans_SNPs = {i: trans_SNPs[i] for i in trans_SNPs if (i not in degen[trans]) or (trans_SNPs[i][0] in degen[trans][i])}              
        else:
            trans_SNPs = []
        for position in combined_dict[trans]:
            #format SNP data for INSIGHT
            if position in mapping_dict[trans]:
                current_string = "\nsite\tchr{0}:{1}".format(chrom, str(mapping_dict[trans][position]))
                try:
                    major_pp = pp_dict[trans][position][base_order.index(seqs[names.index(trans)][position])]
                except KeyError:
                    raise KeyError
                if position not in trans_SNPs:
                    #monomorphic
                    site_type = "M"
                    current_string = current_string + "\t{0}\t{1}".format(site_type, major_pp)
                else:
                    minor_pp = pp_dict[trans][position][base_order.index(trans_SNPs[position][0])]
                    if current_SNPs[trans][position][1] > freq_threshold:
                        #high-frequency SNP
                        site_type = "H"
                    else:
                        #low-frequency SNP
                        site_type = "L"
                    current_string = current_string + "\t{0}\t{1}\t{2}".format(site_type, major_pp, minor_pp)
                if position in pos_dict[trans]:
                    current_neutral_output.append(current_string)
                    site_count = site_count + 1
                    if site_type != "M":
                        #count polymorphic sites
                        K = K + 1
                else:
                    current_hit_output.append(current_string)
    theta = (1/a) * (K/site_count)
    block_header = ["block\tchr{0}:1-300000000\ttheta\t{1}\tlambda\t{2}".format(chrom, theta, lambda_b)]
    current_neutral_output = "".join(block_header + current_neutral_output)
    current_hit_output = "".join(block_header + current_hit_output)
    return(current_hit_output, current_neutral_output)
     
def get_pp(outroot, subst_model, phy_file, model_file, separate_to_concat_mapping, combined_dict, tuples_mapping, min_inf = None, parse_output = True):
    '''
    Get prior probabilities for all the bases at the different sites in an MSA. Note that for phyloFit these
    are posterior probabilities but theyare priors for INSIGHT.
    '''
    #you don't want to compute a tree, just get the posterior probabilities for an existing tree
    #hence all the flags from --post_probs onwards
    arguments = ["phyloFit", "--init-model", model_file, "--out-root", outroot, "--subst-mod", subst_model,
                           "--msa-format", "PHYLIP", "--post-probs", "--scale-only", "--no-rates", "--no-freqs", phy_file]
    if min_inf:
        arguments.extend(["-I", min_inf])
    results = run_process(arguments)
    #parse into convenient dictionary
    if parse_output:
        pp_file = "{0}.postprob".format(outroot)
        pp = rw.read_many_fields(pp_file, " ")
        pp = [[j for j in i if j] for i in pp]
        #the outgroup nodes are labelled from the inside out, starting from 1
        pp = {i[1]: i[-4:] for i in pp}
        pp_final = {}
        #map from coordinates in the concatenated alignment to positions in individual CDSs
        for trans in separate_to_concat_mapping:
            pp_final[trans] = {}
            for position in combined_dict[trans]:
                pp_final[trans][position] = pp[tuples_mapping[separate_to_concat_mapping[trans][position]]]
        return(pp_final)

def get_tuples_mapping(sequences):
    '''
    For each position in an alignment, get a tuple that concatenates the bases from all the different species.
    '''
    tuples = ["".join([j[i] for j in sequences]) for i in range(len(sequences[0]))]
    tuples = [re.sub("N", "*", i) for i in tuples]
    tuples = [re.sub("-", "*", i) for i in tuples]
    tuples_mapping = {pos: i for pos, i in enumerate(tuples)}
    return(tuples_mapping)

def parse_basinhoppin_pos(file):
    '''
    Parse hit/control positions.
    '''
    #not used in main() in the present script but imported into other scripts
    #this is ugly, this function needs to move
    positions = rw.read_many_fields(file, "\t")
    positions = [[i[0], [int(j) for j in i[1].split(",") if j]] for i in positions]
    positions = list_to_dict(positions, 0, 1)
    return(positions)

def parse_degen(file_name):
    '''
    Parse a degenracy file into a nice dictionary with transcript IDs as keys.
    '''
    degen = rw.read_many_fields(file_name, "\t")
    degen = list_to_dict(degen, 0, 1)
    degen = {i: degen[i].split(",") for i in degen}
    for trans in degen:
        separate = [i.split(":") for i in degen[trans]]
        separate = [i for i in separate if len(i) == 2]
        degen[trans] = {int(i[0]): i[1].split("|") for i in separate}
    return(degen)

def parse_INSIGHT_output(INSIGHT_output):
    '''
    Parse the output from INSIGHT.
    '''
    with open(INSIGHT_output) as file:
        output = file.read()
        output_filtered = re.findall("rho[\d\w\.\- \n\t,\[\]: \*]+Post", output)
        estimates = re.findall("(?<=Estimates: )[\d\.\t\- ]+", output_filtered[0])[0]
        estimates = estimates.split(" ")
        try:
            SE = re.findall("(?<=StndrdErr: )[\d\.\t\- ]+", output_filtered[0])[0]
            SE = SE.split(" ")[:3]
        except IndexError:
            SE = ["NA", "NA", "NA"]
        lls = re.findall("(?<=EM status:)[\d \t\-\.e\+]+(?=converged)", output)
        lls = [[j for j in i.split(" ") if j] for i in lls[1:]]
        lls = [abs(float(i[1])) for i in lls]
        lls = [lls[1] - lls[0], lls[2] - lls[0], lls[3] - lls[0]]
        lls = [i * 2 for i in lls]
    output = {"estimates": estimates, "SEs": SE, "chi_sq": lls}
    return(output)

def parse_pos_file(chrom_dict, chrom, pos_file):
    '''
    Parse hits/controls for a particular chromosome.
    '''
    out_dict = {}
    with open(pos_file) as input_file:
        counter = 0
        for line in input_file:
            line = line.split("\t")
            trans = line[0]
            if trans in chrom_dict:
                if chrom_dict[trans] == chrom:
                    out_dict[trans] = [int(i) for i in line[1].split(",")]
            counter = counter + 1
    return(out_dict)

def positions_to_phylip(MSA_raw, CDSs, phylip_data, combined_dict, pos_dict, lengths, phy_file, phy_file_neutral, clean_names, SNP_dict, degen = None):
    '''
    Based on an MSA and a set of hit/control positions, concatenate the alignment at the relevant position
    sites into phylip format suitable for phyloFit.
    '''
    current_CDSs = {}
    #maps from CDS coordinates to coordinates in the final, concatenated alignment
    separate_to_concat_mapping = {}
    cumul_length = 0
    #phylip_data will have everything, phylip_data_neutral will only have control sites
    phylip_data_neutral = copy.deepcopy(phylip_data)
    for trans in sorted(pos_dict.keys()):
        if trans in MSA_raw:
            if trans not in SNP_dict:
                SNP_dict[trans] = []
            #check that you have data from all species
            all_found = []
            for species in phylip_data:
                all_found.extend([species in i["species"] for i in MSA_raw[trans]])
            if False not in all_found:
                #store the human sequence length based on genome annotation data so that you could later
                #check that the sequence is complete (the same length as the human sequence)
                #for all species
                expected_length = lengths[trans]
                human_seq = "".join([i["species"]["homo_sapiens"]["seq"] for i in MSA_raw[trans]])
                human_seq_clean = [i for i in human_seq if i != "-"]
                all_positions = combined_dict[trans]
                if all_positions:
                    #check that teh human sequence you got from the alignment is the expected length
                    #+3 because of the stop
                    if (len(human_seq_clean) + 3) == expected_length:
                        #convert from CDS coordinates to alignment coordinates
                        #you need just the neutral positions for calculating lambda and getting the sequence evolution model
                        #you need all the positions to get Zi values
                        aligned_neutral_positions = conservation.get_aligned_positions([list(human_seq)], pos_dict[trans])
                        aligned_positions = conservation.get_aligned_positions([list(human_seq)], all_positions)
                        all_species_intact = []
                        #revise the expected length to the human alignment length (rather than the length
                        #of the unaligned sequence)
                        revised_expected_length = len(human_seq)               
                        species_temp = {}
                        species_temp_neutral = {}
                        #concatenate the sequence from the different exons
                        #check that it's the expetced length
                        for species in phylip_data:
                            current_seq = "".join([i["species"][species]["seq"] for i in MSA_raw[trans]])
                            current_intact = len(current_seq) == revised_expected_length
                            all_species_intact.append(current_intact)
                            #extract the alignment positions that correspond to hits/controls
                            if current_intact:
                                species_temp[species] = [current_seq[i] for i in aligned_positions]
                                species_temp_neutral[species] = [current_seq[i] for pos, i in enumerate(aligned_neutral_positions) if pos_dict[trans][pos] not in SNP_dict[trans]]
                        #if all species have a complete sequence
                        if False not in all_species_intact:
                            #sanity check
                            if human_seq_clean[:3] != ["A", "T", "G"]:
                                print(trans)
                                print(human_seq_clean)
                                print("Sequence lacks start!")
                                raise Exception
                            if degen:
                                #artificially recode sites as non-divergent if the divergence isn't motif-disrupting
                                species_temp = reencode_disruption(species_temp, all_positions, degen[trans])
                                species_temp_neutral = reencode_disruption(species_temp_neutral, pos_dict[trans], degen[trans], SNP_dict_keys = SNP_dict[trans])
                            #the prior probabilities must be estimated without knowledge of the human sequence
                            #so reencode that as Ns
                            species_temp["homo_sapiens"] = "".join(["N" for i in species_temp["homo_sapiens"]])
                            #concatenate the alignment from this transcript to the full alignment for this chromosome
                            [phylip_data[i].extend(species_temp[i]) for i in phylip_data]
                            [phylip_data_neutral[i].extend(species_temp_neutral[i]) for i in phylip_data_neutral]
                            #store the mappings between CDS coordinates and the cooridnates in the full concatenated alignment
                            separate_to_concat_mapping[trans] = {i: pos + cumul_length for pos, i in enumerate(all_positions)}
                            cumul_length = cumul_length + len(all_positions)
                            current_CDSs[trans] = CDSs[trans]
    #write the final alignment to file
    current_data = ["".join(phylip_data[i]) for i in sorted(phylip_data.keys())]
    current_data_for_writing = MultipleSeqAlignment([SeqRecord(Seq(current_data[pos], IUPAC.ambiguous_dna), id = sorted(clean_names)[pos]) for pos in range(len(current_data))])
    AlignIO.write(current_data_for_writing, phy_file,"phylip-sequential")
    current_data_neutral = ["".join(phylip_data_neutral[i]) for i in sorted(phylip_data_neutral.keys())]   
    current_data_for_writing_neutral = MultipleSeqAlignment([SeqRecord(Seq(current_data_neutral[pos], IUPAC.ambiguous_dna), id = sorted(clean_names)[pos]) for pos in range(len(current_data_neutral))])
    AlignIO.write(current_data_for_writing_neutral, phy_file_neutral,"phylip-sequential")
    return(current_CDSs, separate_to_concat_mapping, current_data, current_data_neutral)

def reduce_dict(pos_dict, fraction):
    '''
    Reduce the sizes of the lists in a hits/controls dictionary to a defined fraction of the initial size by sampling
    bases randomly.
    '''
    output_dict = {i: np.random.choice(pos_dict[i], size = round(fraction * len(pos_dict[i])), replace = False) for i in pos_dict}
    return(output_dict)

def reencode_disruption(species_temp, all_positions, degen, SNP_dict_keys = None):
    '''
    Reencode positions in the alignment to appear as non-divergent if the divergence has been recorded as degenerate.
    '''
    #for the phyloFit stage, you don't want sites that overlap with SNPs
    if SNP_dict_keys:
        all_positions = [i for i in all_positions if i not in SNP_dict_keys]
    #loop over positions
    for pos, site in enumerate(all_positions):
        human_base = species_temp["homo_sapiens"][pos]
        #loop over species in the alignment
        for species in species_temp:
            alt_base = species_temp[species][pos]
            #if the site is divergent
            if alt_base != "-" and alt_base != human_base:
                #if the divergence isn't motif disrupting,
                #reencod the alternative base as the human base
                if alt_base not in degen[site]:
                    species_temp[species][pos] = human_base
    return(species_temp)

def shuffle_dictionaries(dict1, dict2):
    '''
    For negative controls: shuffle the hit/control positions for each transcript,
    maintaining the sizes of the hit/control lists.
    '''
    #combine all sites
    combined_dict = {i: dict1[i] + dict2[i] for i in dict1}
    dict1_out = {}
    dict2_out = {}
    for i in combined_dict:
        current_pos = combined_dict[i].copy()
        #shuffle positions
        random.shuffle(current_pos)
        #divide shuffled positions into pseudo-hits and pseudo-controls
        dict1_out[i] = sorted(current_pos[: len(dict1[i])])
        dict2_out[i] = sorted(current_pos[len(dict1[i]):])
        #sanity check
        if (len(dict1_out[i]) != len(dict1[i])) or (len(dict2_out[i]) != len(dict2[i])):
            print("Dictionary lengths don't match!")
            print(len(dict1[i]))
            print(len(dict1_out[i]))
            print(len(dict2[i]))
            print(len(dict2_out[i]))
    return(dict1_out, dict2_out)

def write_output_file(file_name, output, n):
    '''
    Write INSIGHT input string to file.
    '''
    with open(file_name, "w") as file:
        file.write("samples\t{0}\n".format(n))
        file.write("\n".join(output))
                                      
def main():
    description = "Run INSIGHT on a set of sequences and a set of sites."
    args = parse_arguments(description, ["fasta", "genome", "features_file", "families_file", "suffix", "dataset", "output_folder", "freq_threshold", "n", "hit_file", "control_file", "SNP_file_name_prefix", "CDS_SNP_file_name_prefix", "MSA_file_name_prefix", "trial_file", "trials", "hit_degen_file", "control_degen_file", "hit_reduce", "control_reduce", "new_SNPs", "new_MSA", "shuffle", "nonsyn_hits", "remove_GT", "big_tree"], floats = [7, 18, 19], ints = [8, 15], flags = [20, 21, 22, 23, 24, 25])
    fasta, genome, features_file, families_file, suffix, dataset, general_output_folder, freq_threshold, n, hit_file, control_file, SNP_file_name_prefix, CDS_SNP_file_name_prefix, MSA_file_name_prefix, trial_file, trials, hit_degen_file, control_degen_file, hit_reduce, control_reduce, new_SNPs, new_MSA, shuffle, nonsyn_hits, remove_GT, big_tree = args.fasta, args.genome, args.features_file, args.families_file, args.suffix, args.dataset, args.output_folder, args.freq_threshold, args.n, args.hit_file, args.control_file, args.SNP_file_name_prefix, args.CDS_SNP_file_name_prefix, args.MSA_file_name_prefix, args.trial_file, args.trials, args.hit_degen_file, args.control_degen_file, args.hit_reduce, args.control_reduce, args.new_SNPs, args.new_MSA, args.shuffle, args.nonsyn_hits, args.remove_GT, args.big_tree
    output_folder = "{0}/{1}_{2}".format(general_output_folder, dataset, suffix)

    names, seqs = rw.read_fasta(fasta)

    #prepare feature set and family information
    fs = Feature_Set(features_file, genome)
    fs.set_dataset(dataset)
    if families_file == "None":
        conservation.find_families(fasta, "general/{0}".format(dataset))
        families_file = "general/{0}_families.txt".format(dataset)
    families = rw.read_families(families_file)
    fs.add_families(families)

    make_dir(output_folder)

    general_folder = "DFE/for_everybody"
    make_dir(general_folder)
    if MSA_file_name_prefix == "None":
        MSA_file_name_prefix = "{0}/{1}_MSA".format(general_folder, dataset)

    #read in degeneracy information
    if hit_degen_file != "None":
        degen_hits = parse_degen(hit_degen_file)
        degen_controls = parse_degen(control_degen_file)
    else:
        degen_hits = None
        degen_controls = None

    #get relevant genome features
    transcripts = fs.get_transcripts()
    CDSs = fs.get_CDS()
    lengths = fs.get_lengths(CDSs, CDS = True)
    #filter out sex chromosomes from the analysis
    sex_chromosomes = ["X", "Y"]
    chrom_dict = {i: transcripts[i][0] for i in transcripts if transcripts[i][0] not in sex_chromosomes}
    chroms = list(set(list(chrom_dict.values())))

    clean_names = ["homo", "pan", "pongo", "macaca"]

    #if you're running several trials
    #if just one, it'll still make a single trial file
    if trial_file == "None":
        trial_file = "{0}_{1}_{2}.txt".format(trial_file, suffix, trials)
        

    with open(trial_file, "w") as o_file:
        print(suffix)
        #output file header
        o_file.write("rho\teta\tgamma\tDp\tPw\talpha\ttau\trhose\tetase\tgammase\trholl\tetall\tgammall\n")
        for trial in range(trials):
            print("==========TRIAL {0}==========\n".format(trial))


            #get INSIGHT input data as a string based on divergence and SNP data
            hit_output, neutral_output, chroms_to_keep, hit_counts, control_counts = get_MSA(chroms, chrom_dict, control_file, hit_file, CDSs, lengths, names, seqs, clean_names, freq_threshold, dataset, suffix, genome, output_folder, general_folder, n, SNP_file_name_prefix, CDS_SNP_file_name_prefix, MSA_file_name_prefix, new_SNPs, new_MSA, shuffle, remove_GT, big_tree, hit_reduce = hit_reduce, control_reduce = control_reduce,  degen_hits = degen_hits, degen_controls = degen_controls)

            print("Writing output files...")
            neutral_output_file = "{0}/{1}_{2}_{3}_neutral_input.txt".format(output_folder, dataset, suffix, trial)
            hit_output_file = "{0}/{1}_{2}_{3}_hit_input.txt".format(output_folder, dataset, suffix, trial)
            write_output_file(neutral_output_file, neutral_output, n)
            write_output_file(hit_output_file, hit_output, n)

            print("Running INSIGHT...")
            conservation.INSIGHT(neutral_output_file, hit_output_file, freq_threshold, "../Software/INSIGHT", "{0}_{1}".format(dataset, suffix))

            print("Counting positions on chromosomes...")
            with open("{0}/{1}_{2}_pos_per_chrom.csv".format(output_folder, dataset, suffix), "w") as file:
                file.write("chrom\thits\tcontrols\n")
                for chrom in sorted(chroms_to_keep):
                    file.write("{0}\t{1}\t{2}\n".format(chrom, hit_counts[chrom], control_counts[chrom]))

            INSIGHT_output = "../Software/INSIGHT/{0}_{1}.ins.log".format(dataset, suffix)
            #parse the INSIGHT output and do simple significance testing
            try:
                parsed_output = parse_INSIGHT_output(INSIGHT_output)
                estimates = parsed_output["estimates"]
                SE = parsed_output["SEs"]
                lls = parsed_output["chi_sq"]

                print("\n")
                print("Chisq statistics: {0}".format(" ".join([str(i) for i in lls])))
                rho_pL = scipy.stats.chi2.sf(lls[0], 3)
                print("pL(rho): {0}".format(rho_pL))
                eta_pL = scipy.stats.chi2.sf(lls[1], 1)
                print("pL(eta): {0}".format(eta_pL))
                gamma_pL = scipy.stats.chi2.sf(lls[2], 1)
                print("pL(gamma): {0}".format(gamma_pL))
                
                lls = "\t".join([str(i) for i in lls])
                estimates = "\t".join(estimates)
                SE = "\t".join(SE)
                o_file.write(estimates)
                o_file.write("\t")
                o_file.write(SE)
                o_file.write("\t")
                o_file.write(lls)
                o_file.write("\n")
            #skip trials where INSIGHT failed to produce a full output
            except IndexError:
                print("Skipping...")
                pass

if __name__ == "__main__":
    main()
