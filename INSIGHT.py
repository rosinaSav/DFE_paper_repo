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

base_order = sorted(nc._canon_bases_)
lambda_regex = re.compile("(?<=\(homo:)[\d\.\-e]*(?=\,pan)")

def filter_on_SNPs(current_CDSs, current_SNPs, current_SNPs_to_remove, combined_dict, pos_dict, hit_dict):
    for trans in current_CDSs:
        if trans in current_SNPs:
            trans_SNPs = list(current_SNPs[trans].keys())
        else:
            trans_SNPs = []
        temp_pos = []
        temp_hits = []
        for position in combined_dict[trans]:
            if (trans not in current_SNPs_to_remove) or (position not in current_SNPs_to_remove[trans]) or (position in trans_SNPs):
                if position in pos_dict[trans]:
                    temp_pos.append(position)
                else:
                    temp_hits.append(position)
        pos_dict[trans] = sorted(temp_pos)
        hit_dict[trans] = sorted(temp_hits)
        combined_dict[trans] = sorted(temp_pos + temp_hits)
    return(pos_dict, hit_dict, combined_dict)

def get_lambda(lambda_file_outroot, phy_file, subst_model, min_inf = None):
    lambda_file = "{0}.mod".format(lambda_file_outroot)
    remove_file(lambda_file)
    tree_file = "DFE/UCSC_model.mod"
    arguments = ["phyloFit", "--init-model", tree_file, "--out-root", lambda_file_outroot, "--subst-mod", subst_model,
                           "--msa-format", "PHYLIP", "--scale-only", phy_file]
    if min_inf:
        arguments.extend(["-I", min_inf])
    results = run_process(arguments)
    with open(lambda_file) as file:
        lambda_b = file.read()
    lambda_b = re.findall(lambda_regex, lambda_b)[0]
    return(lambda_b)

def get_MSA(chroms, chrom_dict, control_file, hit_file, CDSs, lengths, names, seqs, clean_names, freq_threshold, dataset, suffix, genome, output_folder, general_folder, n, SNP_file_name_prefix, CDS_SNP_file_name_prefix, MSA_file_name_prefix, new_SNPs, new_MSA, shuffle, remove_GT, big_tree, degen_hits = None, degen_controls = None, hit_reduce = 0, control_reduce = 0):
    subst_model = "JC69"
    a = np.sum([(1/i) for i in range(1, n)])

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
        
        phylip_data = {"homo_sapiens": [], "pongo_abelii": [], "macaca_mulatta": [], "pan_troglodytes": []}
        if big_tree:
            phylip_data = {"gorilla_gorilla": [], "callithrix_jacchus": [], "papio_anubis": [], "chlorocebus_sabaeus": [], "homo_sapiens": [], "pongo_abelii": [], "macaca_mulatta": [], "pan_troglodytes": []}

        print(chrom)

        hit_dict = parse_pos_file(chrom_dict, chrom, hit_file)
        pos_dict = parse_pos_file(chrom_dict, chrom, control_file)

        if shuffle:
            print("Shuffling dictionaries...")
            hit_dict, pos_dict = shuffle_dictionaries(hit_dict, pos_dict)

        if hit_reduce > 0:
            hit_dict = reduce_dict(hit_dict, hit_reduce)
            pos_dict = reduce_dict(pos_dict, control_reduce)

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
            
            if os.path.isfile(SNP_file_name):
                current_new_SNPs = False
            else:
                current_new_SNPs = True
            if os.path.isfile(CDS_SNP_file_name):
                current_new_parse = False
            else:
                current_new_parse = True
            current_SNPs, current_SNPs_to_remove = conservation.get_relative_SNPs(current_CDSs, SNP_file_name, CDS_SNP_file_name, seqs, names, genome, get_new_SNPs = current_new_SNPs, parse_SNPs = current_new_parse, remove_GT = remove_GT)
            #remove positions that overlap the stop codon
            pos_dict = {i: [j for j in pos_dict[i] if (lengths[i] - j) > 3] for i in pos_dict}
            hit_dict = {i: [j for j in hit_dict[i] if (lengths[i] - j) > 3] for i in hit_dict}
            combined_dict = {i: pos_dict[i] + hit_dict[i] for i in hit_dict}
            
            pos_dict, hit_dict, combined_dict = filter_on_SNPs({i: current_CDSs[i] for i in current_CDSs if i in hit_dict}, current_SNPs, current_SNPs_to_remove, combined_dict, pos_dict, hit_dict)

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
            if degen_hits:
                degen = {}
                for trans in degen_hits:
                    degen[trans] = degen_hits[trans].copy()
                    degen[trans].update(degen_controls[trans])
            else:
                degen = None
            current_CDSs, separate_to_concat_mapping, current_data, current_data_neutral = positions_to_phylip(MSA_raw, CDSs, phylip_data, combined_dict, pos_dict, lengths, phy_file, phy_file_neutral, clean_names, current_SNPs, degen = degen)
    ##        os.remove(MSA_file)
            current_neutral_length = len(current_data_neutral[0])
            if current_neutral_length >= 150:
                chroms_to_keep.append(chrom)
                tuples_mapping = get_tuples_mapping(current_data)

                print("Calculating prior probabilities...")
                outroot = "{0}/{1}_{2}_phylip{3}".format(output_folder, dataset, suffix, chrom)
                neutral_outroot = "{0}_neutral".format(outroot)

                lambda_b = get_lambda(neutral_outroot, phy_file_neutral, "JC69")

                model_file = "{0}.mod".format(neutral_outroot)
                pp_final = get_pp(outroot, subst_model, phy_file, model_file, separate_to_concat_mapping, combined_dict, tuples_mapping)
                
                print("Storing information for writing to file later on...")
                control_counts[chrom] = np.sum([len(pos_dict[i]) for i in pos_dict])
                hit_counts[chrom] = np.sum([len(hit_dict[i]) for i in hit_dict])
                current_hit_output, current_neutral_output = get_output_string(pos_dict, hit_dict, pp_final, current_CDSs, current_SNPs, seqs, names, chrom, a, lambda_b, freq_threshold, degen = degen)
                hit_output.append(current_hit_output)
                neutral_output.append(current_neutral_output)
    return(hit_output, neutral_output, chroms_to_keep, hit_counts, control_counts)

def get_output_string(pos_dict, hit_dict, pp_dict, current_CDSs, current_SNPs, seqs, names, chrom, a, lambda_b, freq_threshold, degen = None):
    current_neutral_output = []
    current_hit_output = []
    combined_dict = {i: sorted(pos_dict[i] + hit_dict[i]) for i in pp_dict}
    mapping_dict = map_relative_to_feature(combined_dict, current_CDSs)
    K = 0
    site_count = 0
    for trans in sorted(pp_dict.keys()):
        if trans in current_SNPs:
            trans_SNPs = current_SNPs[trans]
            if degen:
                trans_SNPs = {i: trans_SNPs[i] for i in trans_SNPs if (i not in degen[trans]) or (trans_SNPs[i][0] in degen[trans][i])}              
        else:
            trans_SNPs = []
        for position in combined_dict[trans]:
            if position in mapping_dict[trans]:
                current_string = "\nsite\tchr{0}:{1}".format(chrom, str(mapping_dict[trans][position]))
                try:
                    major_pp = pp_dict[trans][position][base_order.index(seqs[names.index(trans)][position])]
                except KeyError:
                    raise KeyError
                if position not in trans_SNPs:
                    site_type = "M"
                    current_string = current_string + "\t{0}\t{1}".format(site_type, major_pp)
                else:
                    minor_pp = pp_dict[trans][position][base_order.index(trans_SNPs[position][0])]
                    if current_SNPs[trans][position][1] > freq_threshold:
                        site_type = "H"
                    else:
                        site_type = "L"
                    current_string = current_string + "\t{0}\t{1}\t{2}".format(site_type, major_pp, minor_pp)
                if position in pos_dict[trans]:
                    current_neutral_output.append(current_string)
                    site_count = site_count + 1
                    if site_type != "M":
                        K = K + 1
                else:
                    current_hit_output.append(current_string)
    theta = (1/a) * (K/site_count)
    block_header = ["block\tchr{0}:1-300000000\ttheta\t{1}\tlambda\t{2}".format(chrom, theta, lambda_b)]
    current_neutral_output = "".join(block_header + current_neutral_output)
    current_hit_output = "".join(block_header + current_hit_output)
    return(current_hit_output, current_neutral_output)
     
def get_pp(outroot, subst_model, phy_file, model_file, separate_to_concat_mapping, combined_dict, tuples_mapping, min_inf = None, parse_output = True):
    arguments = ["phyloFit", "--init-model", model_file, "--out-root", outroot, "--subst-mod", subst_model,
                           "--msa-format", "PHYLIP", "--post-probs", "--scale-only", "--no-rates", "--no-freqs", phy_file]
    if min_inf:
        arguments.extend(["-I", min_inf])
    results = run_process(arguments)
    if parse_output:
        pp_file = "{0}.postprob".format(outroot)
        pp = rw.read_many_fields(pp_file, " ")
        pp = [[j for j in i if j] for i in pp]
        #the outgroup nodes are labelled from the inside out, starting from 1
        pp = {i[1]: i[-4:] for i in pp}
        pp_final = {}
        for trans in separate_to_concat_mapping:
            pp_final[trans] = {}
            for position in combined_dict[trans]:
                pp_final[trans][position] = pp[tuples_mapping[separate_to_concat_mapping[trans][position]]]
        return(pp_final)

def get_tuples_mapping(sequences):
    tuples = ["".join([j[i] for j in sequences]) for i in range(len(sequences[0]))]
    tuples = [re.sub("N", "*", i) for i in tuples]
    tuples = [re.sub("-", "*", i) for i in tuples]
    tuples_mapping = {pos: i for pos, i in enumerate(tuples)}
    return(tuples_mapping)

def parse_basinhoppin_pos(file):
    #not used in main() in the present script but imported into other scripts
    positions = rw.read_many_fields(file, "\t")
    positions = [[i[0], [int(j) for j in i[1].split(",") if j]] for i in positions]
    positions = list_to_dict(positions, 0, 1)
    return(positions)

def parse_degen(file_name):
    degen = rw.read_many_fields(file_name, "\t")
    degen = list_to_dict(degen, 0, 1)
    degen = {i: degen[i].split(",") for i in degen}
    for trans in degen:
        separate = [i.split(":") for i in degen[trans]]
        separate = [i for i in separate if len(i) == 2]
        degen[trans] = {int(i[0]): i[1].split("|") for i in separate}
    return(degen)

def parse_INSIGHT_output(INSIGHT_output):
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
    current_CDSs = {}
    separate_to_concat_mapping = {}
    cumul_length = 0
    phylip_data_neutral = copy.deepcopy(phylip_data)
    for trans in sorted(pos_dict.keys()):
        if trans in MSA_raw:
            if trans not in SNP_dict:
                SNP_dict[trans] = []
            all_found = []
            for species in phylip_data:
                all_found.extend([species in i["species"] for i in MSA_raw[trans]])
            if False not in all_found:
                expected_length = lengths[trans]
                human_seq = "".join([i["species"]["homo_sapiens"]["seq"] for i in MSA_raw[trans]])
                human_seq_clean = [i for i in human_seq if i != "-"]
                all_positions = combined_dict[trans]
                if all_positions:
                    if (len(human_seq_clean) + 3) == expected_length:
                        #you need just the neutral positions for calculating lambda and getting the sequence evolution model
                        #you need all the positions to get Zi values
                        aligned_neutral_positions = conservation.get_aligned_positions([list(human_seq)], pos_dict[trans])
                        aligned_positions = conservation.get_aligned_positions([list(human_seq)], all_positions)
                        all_species_intact = []
                        revised_expected_length = len(human_seq)               
                        species_temp = {}
                        species_temp_neutral = {}
                        for species in phylip_data:
                            current_seq = "".join([i["species"][species]["seq"] for i in MSA_raw[trans]])
                            current_intact = len(current_seq) == revised_expected_length
                            all_species_intact.append(current_intact)
                            if current_intact:
                                species_temp[species] = [current_seq[i] for i in aligned_positions]
                                species_temp_neutral[species] = [current_seq[i] for pos, i in enumerate(aligned_neutral_positions) if pos_dict[trans][pos] not in SNP_dict[trans]]
                        if False not in all_species_intact:
                            if human_seq_clean[:3] != ["A", "T", "G"]:
                                print(trans)
                                print(human_seq_clean)
                                print("Sequence lacks start!")
                                raise Exception
                            if degen:
                                #artificially recode sites as non-divergent if the divergence isn't motif-disrupting
                                species_temp = reencode_disruption(species_temp, all_positions, degen[trans])
                                species_temp_neutral = reencode_disruption(species_temp_neutral, pos_dict[trans], degen[trans], SNP_dict_keys = SNP_dict[trans])
                            species_temp["homo_sapiens"] = "".join(["N" for i in species_temp["homo_sapiens"]])
                            [phylip_data[i].extend(species_temp[i]) for i in phylip_data]
                            [phylip_data_neutral[i].extend(species_temp_neutral[i]) for i in phylip_data_neutral]
                            separate_to_concat_mapping[trans] = {i: pos + cumul_length for pos, i in enumerate(all_positions)}
                            cumul_length = cumul_length + len(all_positions)
                            current_CDSs[trans] = CDSs[trans]
    current_data = ["".join(phylip_data[i]) for i in sorted(phylip_data.keys())]
    current_data_for_writing = MultipleSeqAlignment([SeqRecord(Seq(current_data[pos], IUPAC.ambiguous_dna), id = sorted(clean_names)[pos]) for pos in range(len(current_data))])
    AlignIO.write(current_data_for_writing, phy_file,"phylip-sequential")
    current_data_neutral = ["".join(phylip_data_neutral[i]) for i in sorted(phylip_data_neutral.keys())]   
    current_data_for_writing_neutral = MultipleSeqAlignment([SeqRecord(Seq(current_data_neutral[pos], IUPAC.ambiguous_dna), id = sorted(clean_names)[pos]) for pos in range(len(current_data_neutral))])
    AlignIO.write(current_data_for_writing_neutral, phy_file_neutral,"phylip-sequential")
    return(current_CDSs, separate_to_concat_mapping, current_data, current_data_neutral)

def reduce_dict(pos_dict, fraction):
    output_dict = {i: np.random.choice(pos_dict[i], size = round(fraction * len(pos_dict[i])), replace = False) for i in pos_dict}
    return(output_dict)

def reencode_disruption(species_temp, all_positions, degen, SNP_dict_keys = None):
    if SNP_dict_keys:
        all_positions = [i for i in all_positions if i not in SNP_dict_keys]
    for pos, site in enumerate(all_positions):
        human_base = species_temp["homo_sapiens"][pos]
        for species in species_temp:
            alt_base = species_temp[species][pos]
            if alt_base != "-" and alt_base != human_base:
                if alt_base not in degen[site]:
                    species_temp[species][pos] = human_base
    return(species_temp)

def shuffle_dictionaries(dict1, dict2):
    combined_dict = {i: dict1[i] + dict2[i] for i in dict1}
    dict1_out = {}
    dict2_out = {}
    for i in combined_dict:
        current_pos = combined_dict[i].copy()
        random.shuffle(current_pos)
        dict1_out[i] = sorted(current_pos[: len(dict1[i])])
        dict2_out[i] = sorted(current_pos[len(dict1[i]):])
        if (len(dict1_out[i]) != len(dict1[i])) or (len(dict2_out[i]) != len(dict2[i])):
            print("Dictionary lengths don't match!")
            print(len(dict1[i]))
            print(len(dict1_out[i]))
            print(len(dict2[i]))
            print(len(dict2_out[i]))
    return(dict1_out, dict2_out)

def write_output_file(file_name, output, n):
    with open(file_name, "w") as file:
        file.write("samples\t{0}\n".format(n))
        file.write("\n".join(output))
                                      
def main():
    description = "Run INSIGHT on a set of sequences and a set of sites."
    args = parse_arguments(description, ["fasta", "genome", "features_file", "families_file", "suffix", "dataset", "output_folder", "freq_threshold", "n", "hit_file", "control_file", "SNP_file_name_prefix", "CDS_SNP_file_name_prefix", "MSA_file_name_prefix", "trial_file", "trials", "hit_degen_file", "control_degen_file", "hit_reduce", "control_reduce", "new_SNPs", "new_MSA", "shuffle", "nonsyn_hits", "remove_GT", "big_tree"], floats = [7, 18, 19], ints = [8, 15], flags = [20, 21, 22, 23, 24, 25])
    fasta, genome, features_file, families_file, suffix, dataset, general_output_folder, freq_threshold, n, hit_file, control_file, SNP_file_name_prefix, CDS_SNP_file_name_prefix, MSA_file_name_prefix, trial_file, trials, hit_degen_file, control_degen_file, hit_reduce, control_reduce, new_SNPs, new_MSA, shuffle, nonsyn_hits, remove_GT, big_tree = args.fasta, args.genome, args.features_file, args.families_file, args.suffix, args.dataset, args.output_folder, args.freq_threshold, args.n, args.hit_file, args.control_file, args.SNP_file_name_prefix, args.CDS_SNP_file_name_prefix, args.MSA_file_name_prefix, args.trial_file, args.trials, args.hit_degen_file, args.control_degen_file, args.hit_reduce, args.control_reduce, args.new_SNPs, args.new_MSA, args.shuffle, args.nonsyn_hits, args.remove_GT, args.big_tree
    output_folder = "{0}/{1}_{2}".format(general_output_folder, dataset, suffix)

    names, seqs = rw.read_fasta(fasta)

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

    if hit_degen_file != "None":
        degen_hits = parse_degen(hit_degen_file)
        degen_controls = parse_degen(control_degen_file)
    else:
        degen_hits = None
        degen_controls = None

    transcripts = fs.get_transcripts()
    CDSs = fs.get_CDS()
    lengths = fs.get_lengths(CDSs, CDS = True)
    sex_chromosomes = ["X", "Y"]
    chrom_dict = {i: transcripts[i][0] for i in transcripts if transcripts[i][0] not in sex_chromosomes}
    chroms = list(set(list(chrom_dict.values())))

    clean_names = ["homo", "pan", "pongo", "macaca"]

    if trial_file == "None":
        trial_file = "{0}_{1}_{2}.txt".format(trial_file, suffix, trials)
        

    with open(trial_file, "w") as o_file:
        print(suffix)
        o_file.write("rho\teta\tgamma\tDp\tPw\talpha\ttau\trhose\tetase\tgammase\trholl\tetall\tgammall\n")
        for trial in range(trials):
            print("==========TRIAL {0}==========\n".format(trial))


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
            except IndexError:
                print("Skipping...")
                pass

if __name__ == "__main__":
    main()
