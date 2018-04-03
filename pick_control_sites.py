'''
Author: Rosina Savisaar.
Pick roughly nucleotide-matched control sites for a set of motif hits.
'''

from bedtools_games import Feature_Set
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
import conservation
import ensembl_ops as eo
from housekeeping import flatten, list_to_dict, make_dir, parse_arguments, pause_script, print_elements, remove_file, run_in_parallel, run_process
import nucleotide_comp as nc
import os
import random
import re
import read_and_write as rw

def get_ancestral_CG(outroot, subst_model, phy_files, model_file, tuples_mapping_dict, anc_CG_file_name, high_CG = None, min_inf = None, macaque = False, comprehensive = False, from_model = False):
    '''
    Get a dictionary that says for each transcript which positions were ancestrally CpG/GpC.
    '''
    #if a file name hasn't been supplied or if the file with the supplied name doesn't exist, determine
    #CpG positions again, otherwise just read them in from the file
    if not anc_CG_file_name or anc_CG_file_name == "None" or not os.path.exists(anc_CG_file_name):
        #you need several in case you have a high_CG dictionary
        pps = []
        for phy_file in phy_files:
            if subst_model == "JC69" or from_model:
                #use an existing substitution model
                arguments = ["phyloFit", "--init-model", model_file, "--out-root", outroot, "--subst-mod", subst_model,
                                       "--msa-format", "PHYLIP", "--post-probs", "--scale-only", phy_file]
            else:
                #estimate a new model
                arguments = ["phyloFit", "--out-root", outroot, "--subst-mod", subst_model,
                                       "--msa-format", "PHYLIP", "--tree", "DFE/full_tree.tree", "--post-probs", phy_file]
                
            if subst_model == "JC69":
                block_size = 4
                tuple_pos_lim = 2
                shift_in_tuple = 0
            else:
                #for dinucleotide models
                block_size = 16
                tuple_pos_lim = 3
                shift_in_tuple = 9

            #turn off when testing                        
            if min_inf:
                arguments.extend(["-I", min_inf])
            results = run_process(arguments)
            #read in posterior probabilities of having various nucelotides ancestrally
            pp_file = "{0}.postprob".format(outroot)
            pp = rw.read_many_fields(pp_file, " ")
            pp = [[j for j in i if j] for i in pp]
            pp = pp[2:]
            #the posterior probability that you had a CpG at a position has to be greater
            #than threshold for a position to be counted as ancestrally CpG
            threshold = 0.5
            #will be over-written if you're doing big tree
            human_pos = 0
            #the outgroup nodes are labelled from the outside in, starting from 1
            if macaque:
                #it's to know whether we're doing big tree or little tree
                if len(pp[0]) == 14:
                    #little tree, mononucleotide
                    pp = {"_".join(i[1:tuple_pos_lim]): [i[len(i) - (3 * block_size): len(i) - (2 * block_size)]] for i in pp}
                elif len(pp[0]) > 14:
                    #big tree/dinucleotide (i.e. it'll give you nonsense if you're trying to do context with the little tree)
                    #the shift_in_tuple is to do with the fact that if you're doing U2S, you want the second tuple and not the first
                    human_pos = 3 + shift_in_tuple
                    if comprehensive:
                        #you want to get all nodes except for node 0, which is the outgroup-ingroup ancestor
                        pp = {"_".join(i[1:tuple_pos_lim]): [i[len(i) - (j * block_size): len(i) - ((j - 1) * block_size)] for j in range(1, 7)] for i in pp}
                    else:
                        pp = {"_".join(i[1:tuple_pos_lim]): [i[len(i) - (6 * block_size): len(i) - (5 * block_size)]] for i in pp}
                else:
                    #for tests etc. where you might only have, say, two species
                    pp = {"_".join(i[1:tuple_pos_lim]): [i[-block_size:]] for i in pp}
            else:
                pp = {"_".join(i[1:tuple_pos_lim]): [i[-block_size:]] for i in pp}
            pps.append(pp)
        anc_CG = {}
        #just to get the length
        example_pp = pps[0][list(pps[0].keys())[0]]
        for trans in tuples_mapping_dict:
            #tuples_mapping_dict has the alignment tuple corresponding to each position
            #because the phyloFit output is organized by tuples, not by positions
            anc_CG[trans] = []
            for node_pos in range(len(example_pp)):
                #if you're using dinucleotides
                if subst_model != "JC69":
                    for pos in sorted(tuples_mapping_dict[trans].keys())[1:]:
                        try:
                            pp_number = 0
                            #if you're gonna produce different output dictionaries for high and low GC regions
                            if high_CG:
                                if pos in high_CG[trans]:
                                    pp_number = 1
                            current_tuple = tuples_mapping_dict[trans][pos]
                            #don't consider positions where there is an alignment gap for human
                            if current_tuple[human_pos] != "*":
##                                print(current_tuple)
##                                print(pps[pp_number])
##                                print("\n")
                                if current_tuple in pps[pp_number]:
                                    current_pp = pps[pp_number][current_tuple][node_pos]
                                else:
                                    current_pp = pps[abs(pp_number - 1)][current_tuple][node_pos]
                                #because it can be either GC or CG, hence 6 or 9
                                if float(current_pp[6]) > threshold or float(current_pp[9]) > threshold:
                                    #you're always testing the second member in the dinucleotide
                                    anc_CG[trans].append(pos - 1)
                                    anc_CG[trans].append(pos)
                        except KeyError:
                            if pos % 100 == 0:
                                pass
                            else:
                                raise KeyError
                else:
                    #if you're using mononucleotides, you have to keep track of what the previous neuclotide was
                    C_prev = False
                    G_prev = False
                    for pos in sorted(tuples_mapping_dict[trans].keys()):
                        pp_number = 0
                        if high_CG:
                            if pos in high_CG[trans]:
                                pp_number = 1
                        current_C = False
                        current_G = False
                        current_tuple = tuples_mapping_dict[trans][pos]
                        if current_tuple[human_pos] != "*":
                            current_pp = pps[pp_number][current_tuple][node_pos]
                            #if current is C and previous was G
                            if float(current_pp[1]) > threshold:
                                if G_prev:
                                    anc_CG[trans].append(G_pos)
                                    anc_CG[trans].append(pos)
                                current_C = True
                            #if current is G and previous was C
                            if float(current_pp[2]) > threshold:
                                if C_prev:
                                    anc_CG[trans].append(C_pos)
                                    anc_CG[trans].append(pos)
                                current_G = True
                            C_prev = False
                            G_prev = False
                            if current_C:
                                C_prev = True
                                #you need to specify the position explicitly because it's not necessarily
                                #the last one if there were dashes
                                C_pos = pos
                            if current_G:
                                G_prev = True
                                G_pos = pos
            anc_CG[trans] = sorted(list(set(anc_CG[trans])))
        remove_file(pp_file)
        if anc_CG_file_name and anc_CG_file_name != "None":
            with open(anc_CG_file_name, "w") as file:
                for trans in anc_CG:
                    to_write = "\t".join([trans, ",".join([str(i) for i in anc_CG[trans]])])
                    file.write(to_write)
                    file.write("\n")
    else:
        #parse
        anc_CG = rw.read_many_fields(anc_CG_file_name, "\t")
        anc_CG = [i for i in anc_CG if len(i) == 2]
        anc_CG = list_to_dict(anc_CG, 0, 1)
        anc_CG = {i: [int(i) for i in anc_CG[i].split(",") if i != ""] for i in anc_CG}
    return(anc_CG)

def get_CpG_dicts(CDSs, chroms, MSA_file_name_prefix, lengths, clean_names, phylip_data, fasta, anc_CG_file_name, high_CG_file_name, fs, macaque_anc = False, pseudoCG = False, comprehensive = False, subst_model = None, return_tuples = False, regions = False):
    '''
    Get two dictionaries, one that says for each transcript which positions are CpG/GpC in macaque
    and one which positions were likely CpG/GpC in the human-macaque ancestor.
    '''
    names, seqs = rw.read_fasta(fasta)
    #if you're gonna determine ancestral CpG positions from scratch rather than reading them in from an existing file
    #if you want to have the name of the file determined automatically
    if (not anc_CG_file_name) or (anc_CG_file_name == "None"):
        new_CG = True
        phy_file = "temp_data/temp_anc_CG{0}.txt".format(random.random())
    #if you want to give the file a name yourself
    elif not os.path.exists(anc_CG_file_name):
        new_CG = True
    else:
        new_CG = False

    if new_CG:
        print("Will get new CpG data...")
        if len(phylip_data) < 8 and comprehensive:
            print("Comprehensive CpG filtering only in big tree mode!")
            raise Exception
        #if you want to pretend some other dinucleotide are CpG
        if pseudoCG:
            CG_kmers = ["C[\-]*T", "A[\-]*G"]
        #the hyphens are there in case the two nucleotides are separated by an indel
        else:
            CG_kmers = ["C[\-]*G", "G[\-]*C"]
        CG_kmers = [re.compile(i) for i in CG_kmers]
        macaque_CG_dict = {}

        anc_CG_concat_full = [[[""]], [[""]]]
        tuples_mapping_dict_full = {}

        for chrom in chroms:

            print(chrom)

            #only leave those CDSs that are on the current chromosome
            current_CDSs = {i: CDSs[i] for i in CDSs if CDSs[i][0][0][0] == chrom}
            coords_file = "temp_data/coords_file{0}.txt".format(random.random())

            #check if the MSA is already at the specified location, otherwise retrieve it
            MSA_file = "{0}_{1}.txt".format(MSA_file_name_prefix, chrom)
            if not os.path.isfile(MSA_file):
                print("Obtaining MSA...")
                eo.get_MSA_gene_list(current_CDSs, coords_file, "EPO", "primates", 85, "homo_sapiens", MSA_file)
                os.remove(coords_file)
                eo.flush_tables("localhost", "mysql", "fackel")
            MSA_raw = eo.parse_MSA_output(MSA_file)
            if high_CG_file_name != "None":
                high_CG = rw.read_many_fields(high_CG_file_name, "\t")
                high_CG = {i[0]: [int(j) for j in i[1:]] for i in high_CG}
            else:
                high_CG = None
            #get concatenated sequences (for determining ancestral CpG positions) and macaque CpG information for this chromosome
            anc_CG_concat, macaque_CG_dict, tuples_mapping_dict = get_CpG_dicts_core(MSA_raw, lengths, phylip_data, CG_kmers, macaque_anc, macaque_CG_dict, high_CG, comprehensive = comprehensive, subst_model = subst_model)
            remove_file(coords_file)
            #add that information to the global dictionaries
            anc_CG_concat_full, tuples_mapping_dict_full = update_anc_CG(anc_CG_concat_full, anc_CG_concat, tuples_mapping_dict_full, tuples_mapping_dict)
            
        phy_files = write_anc_CG(anc_CG_concat_full, anc_CG_file_name, clean_names, macaque_CG_dict)
        pp_file = anc_CG_file_name

    else:
        print("Will read in existing CpG data...")
        pp_file = None
        phy_files = "None"
        high_CG = None
        tuples_mapping_dict_full = None
        macaque_CG_file_name = "{0}_macaque.txt".format(anc_CG_file_name[:-4])
        macaque_CG_dict = rw.read_many_fields(macaque_CG_file_name, "\t")
        macaque_CG_dict = [i for i in macaque_CG_dict if len(i) == 2]
        macaque_CG_dict = list_to_dict(macaque_CG_dict, 0, 1)
        macaque_CG_dict = {i: [int(i) for i in macaque_CG_dict[i].split(",") if i != ""] for i in macaque_CG_dict}
    anc_CG_dict = get_ancestral_CG(pp_file, subst_model, phy_files, "DFE/UCSC_model.mod", tuples_mapping_dict_full, anc_CG_file_name, high_CG = high_CG, macaque = macaque_anc, comprehensive = comprehensive)
    [remove_file(i) for i in phy_files]
    #if you're looking at exon cores/flanks rather than full CDSs
    if regions:
        #you need to have matching bed/fasta files for this to work (with the records in the same order)
        bed = fasta.replace("fasta", "bed")
        transcripts = fs.get_transcripts()
        #for each flank/core, figure out what positions it covers in the full CDS
        mapping_dict = conservation.map_regions_to_CDS(fasta, bed, fs, transcripts, CDSs, trans_ids = True)
        anc_CG_dict = region_CpG(mapping_dict, anc_CG_dict)
    if return_tuples:
        return(anc_CG_dict, macaque_CG_dict, tuples_mapping_dict_full)
    else:
        return(anc_CG_dict, macaque_CG_dict)

def get_CpG_dicts_core(MSA_raw, lengths, phylip_data, CG_kmers, macaque_anc, macaque_CG_dict, high_CG, comprehensive = False, subst_model = None):
    '''
    Core for get_CpG_dicts above.
    '''
    if not subst_model:
        subst_model = "JC69"
    #JC69 is mononucleotide-based, U2S is dinucleotide-based
    if subst_model == "JC69":       
        dinucleotide = False
    else:
        dinucleotide = True
    if high_CG:
        sequence_concat = {0: [], 1: []}
    else:
        sequence_concat = {0: []}
    tuples_mapping_dict = {}
    for trans in MSA_raw:
        #check that the human sequence is present in all the exons
        human_found = ["homo_sapiens" in i["species"] for i in MSA_raw[trans]]
        if False not in human_found:
            expected_length = lengths[trans]
            human_seq = "".join([i["species"]["homo_sapiens"]["seq"] for i in MSA_raw[trans]])
            human_seq_clean = [i for i in human_seq if i != "-"]
            if (len(human_seq_clean) + 3) == expected_length:
                all_species_intact = []
                revised_expected_length = len(human_seq)
                current_seqs = {i: [] for i in phylip_data}
                for species in phylip_data:
                    try:
                        current_seq = "".join([i["species"][species]["seq"] for i in MSA_raw[trans]])
                    except KeyError:
                        current_seq = ""
                    current_seqs[species] = current_seq
                    current_intact = len(current_seq) == revised_expected_length
                    all_species_intact.append(current_intact)
                if False not in all_species_intact:
                    if high_CG:
                        if len(high_CG[trans]) > 0:
                            current_high_CG = realign_high_CG(high_CG[trans], current_seqs["homo_sapiens"])
                        else:
                            current_high_CG = []
                    macaque_sequence = current_seqs["macaca_mulatta"]
                    macaque_CG_pos = [re.finditer(i, macaque_sequence) for i in CG_kmers]
                    macaque_CG_pos = flatten([[(j.start(), j.end()) for j in i] for i in macaque_CG_pos])
                    macaque_CG_pos = flatten([list(range(i[0], i[1])) for i in macaque_CG_pos])
                    current_data = ["".join(current_seqs[i]) for i in sorted(phylip_data.keys())]
                    if high_CG:
                        current_high = ["".join([i for pos, i in enumerate(j) if pos in current_high_CG]) for j in current_data]
                        current_low = ["".join([i for pos, i in enumerate(j) if pos not in current_high_CG]) for j in current_data]
                        sequence_concat[0].append(current_high)
                        sequence_concat[1].append(current_low)
                    else:
                        sequence_concat[0].append(current_data)
                    tuples_mapping_for_CG = get_tuples_mapping(current_data, dinucleotide = dinucleotide)
                    human_aligned_dict = conservation.get_unaligned_pos_dict(current_seqs["homo_sapiens"])
                    tuples_mapping_for_CG = {human_aligned_dict[i]: tuples_mapping_for_CG[i] for i in tuples_mapping_for_CG if human_aligned_dict[i] != None}
                    macaque_CG_pos = sorted(list(set([human_aligned_dict[i] for i in macaque_CG_pos if human_aligned_dict[i] != None])))
                    macaque_CG_dict[trans] = macaque_CG_pos
                    tuples_mapping_dict[trans] = tuples_mapping_for_CG
    return(sequence_concat, macaque_CG_dict, tuples_mapping_dict)

def get_tuples_mapping(sequences, dinucleotide = False):
    '''
    For each position in an alignment, return a string that has the nucleotide at that position in each of the species in the aligment.
    '''
    tuples = ["".join([j[i] for j in sequences]) for i in range(len(sequences[0]))]
    tuples = [re.sub("N", "*", i) for i in tuples]
    tuples = [re.sub("-", "*", i) for i in tuples]
    if dinucleotide:
        #separate the two tuples with an underscore
        old_tuples = tuples.copy()
        tuples = ["_".join([tuples[i - 1], tuples[i]]) for i in range(1, len(tuples))]
        tuples = ["_".join(["".join("*" for i in range(len(old_tuples[0]))), old_tuples[0]])] + tuples
    tuples_mapping = {pos: i for pos, i in enumerate(tuples)}
    return(tuples_mapping)

def realign_high_CG(high_CG, human_seq):
    '''
    Convert high GC positions from coordinates in the unaligned sequence to coordinates in the
    aligned sequence.
    '''
    #the business with blocks is so you'd also include hyphens interleaved with the high CG positions
    CG_blocks = nc.contiguous_blocks(high_CG)
    current_high_CG_temp = conservation.get_aligned_positions([list(human_seq), ""], high_CG)
    current_high_CG = []
    start = 0
    for block in CG_blocks:
        block_length = block[1] - block[0]
        end = start + (block_length) - 1
        if end >= len(current_high_CG_temp):
            end = len(current_high_CG_temp) - 1
        current_high_CG.extend(list(range(current_high_CG_temp[start], current_high_CG_temp[end] + 1)))
        start = start + block_length
    current_high_CG = sorted(current_high_CG)
    return(current_high_CG)

def region_CpG(mapping_dict, anc_CG_dict):
    '''
    Based on mapping_dict, which has the indices where the core/flank maps in the full CDS and
    anc_CG_dict, which has the indices that were ancestrally CpG/GpC in the full CDS,
    determine which positions were ancestrally CpG/GpC in the flank/core.
    '''
    anc_CG_dict_regions = {}
    #loop over the cores/flanks
    for region in mapping_dict:
        #check what CDS it corresponds to
        trans_name = mapping_dict[region]["idn"]
        if trans_name in anc_CG_dict:
            #check which of the positions in the flank overlap with ancestral CpG positions (you check using the full CDS indices but then use enumerate to get back flank/core indices)
            current_CpG = [pos for pos, i in enumerate(mapping_dict[region]["flank indices in CDS"]) if i in anc_CG_dict[trans_name]]
            anc_CG_dict_regions[region] = current_CpG
    return(anc_CG_dict_regions)

def update_anc_CG(anc_CG_concat_full, anc_CG_concat, tuples_mapping_dict_full, tuples_mapping_dict):
    '''
    Update the full concatenated string (for determining ancestral CpG positions) with that from a single chromosome.
    '''
    for pos, key in enumerate(list(anc_CG_concat.keys())):
        [anc_CG_concat_full[pos].append(i) for i in anc_CG_concat[key]]
    for trans in tuples_mapping_dict:
        tuples_mapping_dict_full[trans] = tuples_mapping_dict[trans]
    return(anc_CG_concat_full, tuples_mapping_dict_full)

def write_anc_CG(anc_CG_concat_full, anc_CG_file_name, clean_names, macaque_CG_dict, test_only = False):
    '''
    Write alignment for determining ancestral positions to file.
    '''
    #remove the dummy when not using a high_CG_file
    anc_CG_concat_full = [i for i in anc_CG_concat_full if i != [['']]]
    phy_files = ["temp_data/temp_phy{0}.phy".format(random.random()) for i in anc_CG_concat_full]
    pp_file = anc_CG_file_name

    #this is to get rid of the empty dummy at the very start of each value
    anc_CG_concat_full = [[j for j in i if len(j) > 1] for i in anc_CG_concat_full]
    anc_CG_concat_full = [["".join([k[j] for k in i]) for j in range(len(i[0]))] for i in anc_CG_concat_full]
    if test_only:
        return(anc_CG_concat_full)
    anc_CG_concat_full = [MultipleSeqAlignment([SeqRecord(Seq(i[pos], IUPAC.ambiguous_dna), id = sorted(clean_names)[pos]) for pos in range(len(i))]) for i in anc_CG_concat_full]
    [AlignIO.write(anc_CG_concat_full[i], phy_files[i],"phylip-sequential") for i in range(len(anc_CG_concat_full))]
    if anc_CG_file_name:
        macaque_CG_file_name = "{0}_macaque.txt".format(anc_CG_file_name[:-4])
        with open(macaque_CG_file_name, "a") as file:
            for trans in macaque_CG_dict:
                to_write = "\t".join([trans, ",".join([str(i) for i in macaque_CG_dict[trans]])])
                file.write(to_write)
                file.write("\n")
    return(phy_files)
    
def main():
        description = "Pick roughly nucleotide-matched control sites for a set of motif hits."
        args = parse_arguments(description, ["fasta", "genome", "features_file", "families_file", "dataset", "motifs_file", "run_number", "hit_file", "niter", "stepsize", "control_file", "error_file", "MSA_file_name_prefix", "anc_CG_file_name", "high_CG_file_name", "exclude_file", "brute_mapping", "verbose", "old_motif_format", "nonsyn_hits", "top_set_only", "remove_GT", "leave_CG", "remove_ancestral_CpG", "replacement_control", "macaque_anc", "remove_macaque_CpG", "big_tree", "pseudoCG", "comprehensive", "context", "prone_sites", "CG_gene_filter", "match_size", "raw", "regions"], ints = [6, 8, 9], flags = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35])
        fasta, genome, features_file, families_file, dataset, motifs_file, run_number, hit_file, niter, stepsize, control_file, error_file, MSA_file_name_prefix, anc_CG_file_name, high_CG_file_name, exclude_file, brute_mapping, verbose, old_motif_format, nonsyn_hits, top_set_only, remove_GT, leave_CG, remove_ancestral_CpG, replacement_control, macaque_anc, remove_macaque_CpG, big_tree, pseudoCG, comprehensive, context, prone_sites, CG_gene_filter, match_size, raw, regions = args.fasta, args.genome, args.features_file, args.families_file, args.dataset, args.motifs_file, args.run_number, args.hit_file, args.niter, args.stepsize, args.control_file, args.error_file, args.MSA_file_name_prefix, args.anc_CG_file_name, args.high_CG_file_name, args.exclude_file, args.brute_mapping, args.verbose, args.old_motif_format, args.nonsyn_hits, args.top_set_only, args.remove_GT, args.leave_CG, args.remove_ancestral_CpG, args.replacement_control, args.macaque_anc, args.remove_macaque_CpG, args.big_tree, args.pseudoCG, args.comprehensive, args.context, args.prone_sites, args.CG_gene_filter, args.match_size, args.raw, args.regions

        #argparse can't do booleans
        if anc_CG_file_name == "None":
            anc_CG_file_name = None

        #I store motif data in one of two formats
        if old_motif_format:
            motifs = rw.read_names(motifs_file)[1:]
        else:
            motifs = rw.read_motifs(motifs_file)
            #if you're doing RBP motifs and only want motifs that were found to be enriched in Savisaar and Hurst 2017
            if top_set_only:
                summary_data = rw.read_many_fields("RBP/RBP_hg38_introncontaining_new.txt", "\t")

                summary_dict = list_to_dict(summary_data, 0, 4, floatify = True)

                motifs = {RBP: motifs[RBP] for RBP in motifs if (summary_dict[RBP] < 0.1)}
            motifs = list(set(flatten(list(motifs.values()))))

        #create an instance of a Feature_Set object and associate a structure of paralogous families to it, unless if you've said to ignore that (used when analyzing exon flanks/cores)
        fs = Feature_Set(features_file, genome)
        fs.set_dataset(dataset)
        if families_file == "None":
            conservation.find_families(fasta, "general/{0}".format(dataset))
            families_file = "general/{0}_families.txt".format(dataset)

        if families_file != "ignore":
            families = rw.read_families(families_file)
            fs.add_families(families)

        general_folder = "DFE/for_everybody"
        make_dir(general_folder)
        #if you've already retrieved MSAs from ensembl
        if MSA_file_name_prefix == "None":
            MSA_file_name_prefix = "{0}/{1}_MSA".format(general_folder, dataset)

        #admin
        transcripts = fs.get_transcripts()
        CDSs = fs.get_CDS()
        lengths = fs.get_lengths(CDSs, CDS = True)
        #only consider genes that are not on the sex chromosomes
        sex_chromosomes = ["X", "Y"]
        chrom_dict = {i: transcripts[i][0] for i in transcripts if transcripts[i][0] not in sex_chromosomes}
        chroms = list(set(list(chrom_dict.values())))

        #U2S is a dinucleotide-based substitution model, JC69 is mononucleotide-based
        if context:
            subst_model = "U2S"
        else:
            subst_model = "JC69"

        #names used in the MSA (there's a character restriction in the phylip files so you can't use the full name)
        clean_names = ["homo", "pan", "pongo", "macaca"]
        phylip_data = {"homo_sapiens": [], "pongo_abelii": [], "macaca_mulatta": [], "pan_troglodytes": []}
        if big_tree:
            clean_names = ["calli", "chloro", "gorilla", "homo", "macaca", "pan", "papio", "pongo"]
            phylip_data = {"gorilla_gorilla": [], "callithrix_jacchus": [], "papio_anubis": [], "chlorocebus_sabaeus": [], "homo_sapiens": [], "pongo_abelii": [], "macaca_mulatta": [], "pan_troglodytes": []}

        if remove_ancestral_CpG or remove_macaque_CpG or CG_gene_filter:
            anc_CG_dict, macaque_CG_dict = get_CpG_dicts(CDSs, chroms, MSA_file_name_prefix, lengths, clean_names, phylip_data, fasta, anc_CG_file_name, high_CG_file_name, fs, macaque_anc = macaque_anc, pseudoCG = pseudoCG, comprehensive = comprehensive, subst_model = subst_model, regions = regions)
                                
        else:
            anc_CG_dict = None
            macaque_CG_dict = None
        
        if replacement_control:
            nc.fit_control_pos_to_hits_replacement(fasta, motifs, run_number, hit_file, control_file, anc_CG_dict, macaque_CG_dict, family_seed = 5, CG_gene_filter = CG_gene_filter, niter = niter, verbose = verbose, brute_mapping = brute_mapping, stepsize = stepsize, write_errors = error_file, fs = fs, nonsyn_hits = nonsyn_hits, leave_CG = leave_CG, remove_ancestral_CpG = remove_ancestral_CpG, remove_macaque_CpG = remove_macaque_CpG, pseudoCG = pseudoCG, prone_sites = prone_sites, match_size = match_size, raw = raw, exclude_file = exclude_file)
        else:
            nc.fit_control_pos_to_hits_wrapper(fasta, motifs, run_number, hit_file, control_file, anc_CG_dict, macaque_CG_dict, family_seed = 5, CG_gene_filter = CG_gene_filter, niter = niter, verbose = verbose, brute_mapping = brute_mapping, stepsize = stepsize, write_errors = error_file, fs = fs, nonsyn_hits = nonsyn_hits, leave_CG = leave_CG, remove_ancestral_CpG = remove_ancestral_CpG, remove_macaque_CpG = remove_macaque_CpG, pseudoCG = pseudoCG, prone_sites = prone_sites, match_size = match_size)


if __name__ == "__main__":
    main()
