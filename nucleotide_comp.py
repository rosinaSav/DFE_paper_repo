'''
Author: Rosina Savisaar.
Module containing functions that have to do with
nucleotide composition, k-mer frquencies etc.
'''

from Bio.SeqUtils.CodonUsage import SynonymousCodons as _genetic_code_
import csv
from cython_func import calc_density_c, calc_density_for_concat_c, calc_density_for_concat_several_c
from housekeeping import flatten, print_elements, remove_file, run_in_parallel, run_process, overlap, update_counter
import my_stats as ms
import numpy as np
import operator
import os
import random
import re
import read_and_write as rw
from scipy.optimize import basinhopping
import scipy.stats
import itertools as it
import time

#global variables
_IUPAC_dict_ = {"AG": "R", "CT": "Y", "CU": "Y", "CG": "S", "AT": "W", "AU": "W", "GT": "K", "GU": "K", "AC": "M",
              "CGT": "B", "CGU": "B", "AGT": "D", "AGU": "D", "ACT": "H", "ACU": "H", "ACG": "V", "ACGT": "N", "ACGU": "N",
              "A": "A", "T": "T", "U": "U", "G": "G", "C": "C"}

_IUPAC_dict_reverse_ = {'B': 'CGT', 'N': 'ACGT', 'K': 'GT', 'Y': 'CT', 'M': 'AC',
                      'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'S': 'CG', 'W': 'AT', 'R': 'AG', "A": "A",
                      "T": "T", "U": "U", "G": "G", "C": "C"}

fourfold_deg_list = fourfold_deg_list()

#over-all composition of hg38 from get_base_comp(bases, "Genomes/hg38/hg38.fa")
_genome_comp_dict_ = {'A': 0.295185715129424, 'G': 0.2047602545425745, 'T': 0.2961350003087717, 'C': 0.2039190300192298}
_canon_bases_ = ["A", "T", "C", "G"]

#mapping from codons to amino acids
_genetic_code_reverse_ = {}
for aa in _genetic_code_:
    for codon in _genetic_code_[aa]:
        _genetic_code_reverse_[codon] = aa

_rev_comp_dict_ = {}
_rev_comp_dict_["A"] = "T"
_rev_comp_dict_["T"] = "A"
_rev_comp_dict_["C"] = "G"
_rev_comp_dict_["G"] = "C"

def calc_nt_comp_error(counts, hit_mononts):
    '''
    Given the nucleotide composition of the true hits and counts the different bases in control sequence,
    calculate the error in mononucleotide composition.
    '''
    counts = [round(i) for i in counts]
    total = np.sum(counts)
    if total:
        sim_mononts = {i: counts[pos]/total for pos, i in enumerate(sorted(list(hit_mononts.keys())))}
        return(np.sum([abs(sim_mononts[i] - hit_mononts[i]) for i in sim_mononts]))
    else:
        #if all the counts are 0
        return(np.inf)

def calc_nt_freqs(motifs, order, return_array = False, no_phase = False, count = False):
    '''
    Calculate the order-nucleotide frequencies of a set of motifs.
    no_phase if you only want to consider the first reading frame (no overlaps between kmers)
    returns numpy array if return_array, else returns dictionary
    '''
    all_kmers = generate_all_kmers(order)
    all_kmers_in_motifs = []
    if no_phase:
        phase_range = 1
    else:
        phase_range = order
    #loop over the different possible reading frames
    for phase in range(phase_range):
        for motif in motifs:
            current_strings = [motif[i:i+order] for i in range(phase, len(motif), order)]
            #so it wouldn't split the k-mer if it's halfway through the kmer when it reaches the end of the motif
            current_strings = [i for i in current_strings if (len(i) == order) and ("N" not in i)]
            all_kmers_in_motifs.extend(current_strings)
    kmer_number = len(all_kmers_in_motifs)
    if return_array:
        if count:
            result_array = np.array([all_kmers_in_motifs.count(kmer) for kmer in all_kmers])
        else:
            result_array = np.array([all_kmers_in_motifs.count(kmer)/kmer_number for kmer in all_kmers])
        if abs(np.sum(result_array) - 1) > 0.0001 :
            print(result_array)
            print("Frequencies don't sum up to 1!")
            raise Exception
        return(result_array)
    freqs_dict = {}
    if kmer_number:
        for kmer in all_kmers:
            if count:
                freqs_dict[kmer] = all_kmers_in_motifs.count(kmer)
            else:
                freqs_dict[kmer] = all_kmers_in_motifs.count(kmer)/kmer_number
    else:
        return(None)
    if not count:
        if abs(np.sum(list(freqs_dict.values())) - 1) > 0.0001:
                print(freqs_dict)
                print(motifs)
                print("Frequencies don't sum up to 1!")
                raise Exception        
    return(freqs_dict)

def calc_nt_freqs_sequence(fasta, order):
    '''
    Calculate the order-mer frequencies of a set of sequences.
    '''
    kmers = generate_all_kmers(order)
    names, seqs = rw.read_fasta(fasta)
    print(fasta)
    output = {i: {} for i in names}
    for name in names:
        for kmer in kmers:
            output[name][kmer] = 0
    for pos, seq in enumerate(seqs):
        if pos % 1000 == 0:
            print(pos)
        seq = seq.upper()
        current_name = names[pos]
        total = 0
        for index in range(len(seq)):
            try:
                current_kmer = seq[index: index + order]
                if current_kmer in kmers:
                    output[current_name][current_kmer] = output[current_name][current_kmer] + 1
                    total = total + 1
            except IndexError:
                pass
        if total > 0:
            output[current_name] = {i: output[current_name][i]/total for i in kmers}
    return(output)

def check_amount_of_information(input_dict, lengths_dict, sim_file_folder, sim_file_suffix, fraction, threshold, fs = None):
    '''
    Check that either the true density is greater than _threshold_ or at least _fraction_ of the simulants have a density greater than _threshold_.
    input_dict has RBPs as keys and the raw density of the motifs as values.
    lengths_dict has sequence identifiers as keys and the lengths of the corresponding CDSs in bp as values.
    fs is a Feature_Set object from bedtools_games. Necessary if you want to average over families.
    '''
    output_dict = {}
    for protein in input_dict:
        if input_dict[protein] > threshold:
            output_dict[protein] = True
        else:
            input_file_name = "{0}/{1}{2}".format(sim_file_folder, protein, sim_file_suffix)
            current_sim_densities = parse_sim_densities(lengths_dict, input_file_name, fs = fs)
            above_threshold = [i for i in current_sim_densities if i > threshold]
            current_fraction = len(above_threshold)/len(current_sim_densities)
            if current_fraction >= fraction:
                output_dict[protein] = True
            else:
                output_dict[protein] = False
    return(output_dict)

def codon_diversity(CDS):
    '''
    Given a CDS, returns the number of unique codons divided by the total number of codons. A crappy measure because the
    number of unique codons is expected to saturate. You should use windows or something.
    '''
    if len(CDS) % 3 != 0:
        print("The length of the CDS must be a multiple of 3!")
        raise Exception
    all_codons = [CDS[i:i+3] for i in range(0, len(CDS), 3)]
    unique_codons = list(set(all_codons))
    print("\n")
    print(all_codons)
    print(unique_codons)
    return(len(unique_codons)/len(all_codons))

def codon_positions(fasta, out_file, hit_file, human = False, ancestral = None, mono = False, trint = False, exclude_file = None):
    '''
    Makes a list of all the fourfold-degenerate codons in a fasta file and stores their positions.
    '''
    names, seqs = rw.read_fasta(fasta)
    CG_2mers = ["CG", "GC"]
    CG_lengths = [2, 2]
    CG_regex = motif_to_regex(CG_2mers)
    if hit_file:
        hits = rw.read_pos(hit_file)
    else:
        hits = names.copy()
    if ancestral:
        ancestral = rw.read_pos(ancestral)
    if exclude_file:
        exclude = rw.read_pos(exclude_file)
    out_dict = {}
    counter = 0
    for pos, seq in enumerate(seqs):
        counter = update_counter(counter, 1000)
        if names[pos] in hits:
            current_4f = get_4fold_deg(seq)
            if hit_file:
                current_hits = hits[names[pos]]
            else:
                current_hits = []
            if human:
                CG_pos = get_motif_set_density(CG_regex, CG_lengths, seq, concat = True)["positions"]
            else:
                CG_pos = []
            if ancestral:
                anc_pos = ancestral[names[pos]]
            else:
                anc_pos = []
            if exclude_file:
                exclude_pos = exclude[names[pos]]
            else:
                exclude_pos = []
            for codon_pos in current_4f:
                if (codon_pos not in current_hits) and (codon_pos not in CG_pos) and (codon_pos not in anc_pos) and (codon_pos not in exclude_pos):
                    if mono:
                        #if you're doing this with single nucleotides rather than full codons
                        current_codon = seq[codon_pos]
                    elif trint:
                        current_codon = "".join(seq[codon_pos - 1 : codon_pos + 2])                        
                    else:
                        current_codon = "".join(seq[codon_pos - 2 : codon_pos + 1])
                    if current_codon not in out_dict:
                        out_dict[current_codon] = {}
                    if names[pos] not in out_dict[current_codon]:
                        out_dict[current_codon][names[pos]] = []
                    out_dict[current_codon][names[pos]].append(str(codon_pos))
    with open(out_file, "w") as file:
        for codon in sorted(list(out_dict.keys())):
            file.write(codon)
            for name in sorted(list(out_dict[codon].keys())):
                file.write("|{0},".format(name))
                file.write(",".join(out_dict[codon][name]))
            file.write("\n")

def codon_usage_single(sequence, limit, reverse = False):
    '''
    Given a sequence, determine the codon (and thus amino acid) used at each position)
    limit: how far into the sequence should you go (in bp)
    if reverse, start from the end of the sequence
    '''
    if len(sequence) < limit:
        print("The sequence has to be at least {0} bp long!".format(limit))
        raise Exception
    if len(sequence)%3 != 0:
        print("The length of the sequence is not a multiple of 3!")
        raise Exception
    if not reverse:
        sequence = sequence[:limit]
    else:
        sequence = sequence[(len(sequence) - limit): len(sequence)]
    sequence_length = len(sequence)/3
    result = {}
    counter = 0
    for i in range(0, len(sequence), 3):
        codon = sequence[i: i + 3]
        aa = _genetic_code_reverse_[codon]
        result[counter] = [aa, codon]
        counter = counter + 1
    if reverse:
        result = {int(sequence_length - i - 1): result[i] for i in result}
    return(result)

def codon_usage_matrix(sequences, limit, reverse = False):
    '''
    Given a set of sequences, calculate codon usage for each amino acid at each position.
    limit: how far into the sequence should you go (in bp)
    if reverse, start from the end of the sequence
    '''
    sequences = [i for i in sequences if len(i) >= limit]
    if limit%3 != 0:
        print("The required analysis length has to be a multiple of 3!")
        raise Exception
    codon_limit = int(limit/3)
    result = {}
    for site in range(codon_limit):
        result[site] = {}
        for aa in _genetic_code_:
            result[site][aa] = {"number": 0, "codons": {}}
            for codon in _genetic_code_[aa]:
                result[site][aa]["codons"][codon] = 0
    for seq in sequences:
        current_data = codon_usage_single(seq, limit, reverse)
        for site in range(codon_limit):
            result[site][current_data[site][0]]["number"] = result[site][current_data[site][0]]["number"] + 1
            result[site][current_data[site][0]]["codons"][current_data[site][1]] = result[site][current_data[site][0]]["codons"][current_data[site][1]] + 1
    for site in range(codon_limit):
        for aa in _genetic_code_:
            current_number = result[site][aa]["number"]
            if current_number == 0:
                result[site][aa] = None
            else:
                for codon in _genetic_code_[aa]:
                    result[site][aa]["codons"][codon] = result[site][aa]["codons"][codon]/current_number
    return(result)

def consensus_from_PWM(PWM, bases, threshold, transform = False):
    '''
    Take a PWM/PSSM and convert it into a consensus sequence.
    '''
    PWM = np.array(PWM)
    #you want sites in the columns and bases in the rows
    if transform:
        PWM = PWM.transpose()
    consensus = ""
    #loop over the sites
    for site in range(PWM.shape[1]):
        #pick out the bases whose frequency/count/whatever surpasses the threshold
        enriched = np.where(PWM[:,site] >= threshold)
        enriched_bases = bases[enriched]
        #if there's several, get the IUPAC ambiguous nucleotide
        if len(enriched_bases) > 1:
            IUPAC_base = get_ambiguous_nucl(enriched_bases)
            consensus = consensus + IUPAC_base
        else:
            consensus = consensus + enriched_bases[0]
    return(consensus)

def contiguous_blocks(positions):
    '''
    Find contiguous sequences in a list of integers.    
    '''
    first_cluster_start = [positions[0]]
    last_cluster_end = [positions[-1]]
    cluster_starts = [j for i,j in enumerate(positions[1:]) if j - positions[i] != 1]
    first_cluster_start.extend(cluster_starts)
    cluster_starts = first_cluster_start
    cluster_ends = [j for i,j in enumerate(positions[:-1]) if positions[i+1] - j != 1]
    cluster_ends.extend(last_cluster_end)
    #to make Python-style ranges
    cluster_ends = [i + 1 for i in cluster_ends]
    cluster_ranges = zip(cluster_starts, cluster_ends)
    cluster_ranges = [(i[0], i[1]) for i in cluster_ranges]
    return(cluster_ranges)

def count_codons(names, seqs, hits, mono = False, trint = False):
    '''
    Given a set of sequences and a hit file, count how many times each codon appears at each of the hit positions.
    '''
    output_codons = {}
    for pos, seq in enumerate(seqs):
        if names[pos] in hits:
            current_hits = hits[names[pos]]
            for position in current_hits:
                if mono:
                    codon = seq[position]
                elif trint:
                    codon = seq[position - 1: position + 2]
                else:
                    codon = seq[position - 2: position + 1]
                if codon not in output_codons:
                    output_codons[codon] = 0
                output_codons[codon] = output_codons[codon] + 1
    return(output_codons)

def DNA_RNA_conversion(sequence):
    '''
    Convert between DNA and RNA sequences (Us to Ts or vice versa).
    '''
    if "U" in sequence:
        sequence = re.sub("U", "T", sequence)
    elif "T" in sequence:
        sequence = re.sub("T", "U", sequence)
    return(sequence)

def fit_control_pos_to_hits_replacement(fasta, motifs, run_number, hit_file, control_file, anc_CG, macaque_CG, raw = False, niter = None, exclude_file = None, verbose = False, verbose_detailed = False, brute_mapping = False, stepsize = None, fs = None, write_errors = False, alt_anc_CGs = None, nonsyn_hits = False, tuples_mapping = None, family_seed = None, leave_CG = False, remove_ancestral_CpG = False, remove_macaque_CpG = False, CG_gene_filter = False, pseudoCG = False, prone_sites = False, match_size = False):
    '''
    Pick positions from a control sequence so as to match the nucleotide composition to a particular configuration, using sampling with replacement.
    '''
    names, seqs = rw.read_fasta(fasta)
    true_pos_dict, hit_mononts_dict, counts_dict, true_lengths_dict, monont_pos_dict, bounds_dict = prepare_control_pos_for_hits(names, seqs, motifs, anc_CG, macaque_CG, brute_mapping = False, verbose_detailed = verbose_detailed, fs = fs, nonsyn_hits = nonsyn_hits, tuples_mapping = tuples_mapping, alt_anc_CGs = alt_anc_CGs, family_seed = family_seed, leave_CG = leave_CG, remove_ancestral_CpG = remove_ancestral_CpG, CG_gene_filter = CG_gene_filter, remove_macaque_CpG = remove_macaque_CpG, exclude_file = exclude_file, pseudoCG = False, match_size = match_size, prone_sites = prone_sites)

    print("Writing hit positions to file...")

    names_to_keep = []
    counter = 0
    with open(hit_file, "w") as file:
        for name in names:
            if name in counts_dict and counts_dict[name]:
                counter = counter + 1
                names_to_keep.append(name)
                file.write("{0}\t{1}\n".format(name, ",".join([str(i) for i in true_pos_dict[name]])))
    names = names_to_keep

    print("{0} sequences were kept.".format(counter))

    print("Writing control positions to file...")
    with open(control_file, "w") as file:
        for name in names:
            if monont_pos_dict[name] and true_pos_dict[name]:
                #if you want to simply take all possible control positions with no nucleotide composition matching
                if raw:
                    #mash together and sort the positions across the four nucleotides
                    control_pos = sorted(flatten([monont_pos_dict[name][i] for i in monont_pos_dict[name]]))
                else:
                    #you need to know how many of each base you have at the hit positions
                    #but in hit_mononts_dict you have frequencies rather than counts
                    #so you multiply each frequency by the total number of hit positions and convert to integer to get counts
                    hits_per_base = {i: int(hit_mononts_dict[name][i] * len(true_pos_dict[name])) for i in hit_mononts_dict[name]}
                    #sample with repacement to get as many control positions of each base as there are hit psoitions
                    #unless if there is 0 of a particular base, in which case you obviously sample none of that base
                    control_pos = {i: np.random.choice(monont_pos_dict[name][i], hits_per_base[i], replace = True) if monont_pos_dict[name][i] else [] for i in hits_per_base}
                    #merge the list across the different bases but NB! don't uniquify!
                    control_pos = sorted(flatten(list(control_pos.values())))
                file.write("{0}\t{1}\n".format(name, ",".join([str(i) for i in control_pos])))

def fit_control_pos_to_hits(seq_names, counts_dict, hit_mononts_dict, bounds_dict, run, true_lengths_dict, options_supplied):
    '''
    Pick positions from a control sequence so as to match the nucleotide composition to a particular configuration.
    '''
    lengths = []
    counts = []
    errors = []
    result_names = []
    #defaults
    options = {"stepsize": 10, "verbose": False, "niter": 500}
    #override defaults
    for key in options_supplied:
        options[key] = options_supplied[key]

    counter = 0
    for seq_name in seq_names:
        initial_counts = counts_dict[seq_name]
        if initial_counts:
            hit_mononts = hit_mononts_dict[seq_name]
            bounds = bounds_dict[seq_name]
            true_length = true_lengths_dict[seq_name]
            counter = counter + 1
            np.random.seed()
            minim_result = basinhopping(calc_nt_comp_error, initial_counts, minimizer_kwargs = {"bounds": bounds, "args": (hit_mononts,)}, stepsize = options["stepsize"], niter = options["niter"])
            try:
                best_counts = [round(i) for i in minim_result.x]
            except ValueError:
                best_counts = minim_result.x
            counts.append(best_counts)
            length = np.sum(best_counts)
            relative_length = length/true_length
            lengths.append(relative_length)
            frequencies = {i: best_counts[pos]/length for pos, i in enumerate(sorted(list(hit_mononts.keys())))}
            error = [hit_mononts[i] - frequencies[i] for i in sorted(list(hit_mononts.keys()))]
            errors.append(error)
            result_names.append(seq_name)
            if counter % 100 == 0:
                print("{0}\\\\\\{1}".format(run, counter))
            if options["verbose"]:
                print(best_counts)
                print("C: {0} vs {1}.".format(hit_mononts["C"], frequencies["C"]))
                print("A: {0} vs {1}.".format(hit_mononts["A"], frequencies["A"]))
                print("G: {0} vs {1}.".format(hit_mononts["G"], frequencies["G"]))
                print("T: {0} vs {1}.".format(hit_mononts["T"], frequencies["T"]))
                print("ERROR:")
                print(np.mean([abs(i) for i in error]))
                print("LENGTH:")
                print(relative_length)
                print("\n")
    return(result_names, counts, lengths, errors)

def fit_control_pos_to_hits_wrapper(fasta, motifs, run_number, hit_file, control_file, anc_CG, macaque_CG, niter = None, verbose = False, verbose_detailed = False, brute_mapping = False, stepsize = None, fs = None, write_errors = False, alt_anc_CGs = None, nonsyn_hits = False, tuples_mapping = None, family_seed = None, leave_CG = False, remove_ancestral_CpG = False, remove_macaque_CpG = False, CG_gene_filter = False, pseudoCG = False, prone_sites = False, match_size = False):
    '''
    Pick control positions to match hit positions in nucleotide composition.
    '''

    names, seqs = rw.read_fasta(fasta)

    #prepare data: where the true hits are, how many of each mononucleotide they contain,
    #what nucleotides you have available at potential control sites and what the relevant positions are
    true_pos_dict, hit_mononts_dict, counts_dict, true_lengths_dict, monont_pos_dict, bounds_dict = prepare_control_pos_for_hits(names, seqs, motifs, anc_CG, macaque_CG, brute_mapping = False, verbose_detailed = verbose_detailed, fs = fs, nonsyn_hits = nonsyn_hits, tuples_mapping = tuples_mapping, alt_anc_CGs = alt_anc_CGs, family_seed = family_seed, leave_CG = leave_CG, remove_ancestral_CpG = remove_ancestral_CpG, CG_gene_filter = CG_gene_filter, remove_macaque_CpG = remove_macaque_CpG, pseudoCG = False, match_size = match_size, prone_sites = prone_sites)

    print("Writing hit positions to file...")

    names_to_keep = []
    counter = 0
    with open(hit_file, "w") as file:
        for name in names:
            if name in counts_dict and counts_dict[name]:
                counter = counter + 1
                names_to_keep.append(name)
                file.write("{0}\t{1}\n".format(name, ",".join([str(i) for i in true_pos_dict[name]])))
    names = names_to_keep

    print("{0} sequences were kept.".format(counter))

    run_dict = {}
    skews = []

    #options for optimization
    options = {"verbose": verbose, "niter": niter, "stepsize": stepsize}

    #we're gonna run the optimization run_number times and pick the run that led to the least skewed divergences from
    #the nucleotide composition of the hits
    #(the optimization algorithm minimizes absolute error but in this second step, we minimize skew)
    for run in range(run_number):
        run_dict[run] = {"names": [], "counts": [], "lengths": [], "errors": []}

        #run optimization once (separately for each transcript, allowing you to parallelize
        #across transcripts)
        results = run_in_parallel(names, ["foo", counts_dict, hit_mononts_dict, bounds_dict, run, true_lengths_dict, options], fit_control_pos_to_hits, workers = 20)
        for result in results:
            current_result = result.get()
            run_dict[run]["names"].extend(current_result[0])
            run_dict[run]["counts"].extend(current_result[1])
            run_dict[run]["lengths"].extend(current_result[2])
            run_dict[run]["errors"].extend(current_result[3])

        #skew, defined as the number of transcripts that have positive error - the number of
        #transcripts that have negative error for the current base, summed across the four bases
        skew = np.sum([abs(len([i[j] for i in run_dict[run]["errors"] if i[j] > 0]) - len([i[j] for i in run_dict[run]["errors"] if i[j] < 0])) for j in range(4)])
        skews.append(skew)

    #pick run with lowest skew
    best_run = skews.index(np.min(skews))
    print(skews)
    print(best_run)
    print("\n")

    if write_errors:
        write_minimizer_error(run_dict[best_run]["errors"], run_dict[best_run]["lengths"], run_dict[best_run]["names"], write_errors)

    with open(control_file, "w") as file:
        for pos, name in enumerate(run_dict[best_run]["names"]):
            current_sim_pos = pick_control_pos(run_dict[best_run]["counts"][pos], monont_pos_dict[name])
            file.write("{0}\t{1}\n".format(name, ",".join([str(i) for i in current_sim_pos])))
    
def fourfold_deg_list():
    '''
    Get a list of all the fourfold degenerate codons.
    '''
    codons = ["TGA","TAG","TAA","TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"]
    aas =  ["STOP","STOP","STOP","Phe", "Phe", "Leu", "Leu", "Leu", "Leu", "Leu", "Leu", "Ile", "Ile", "Ile", "Met", "Val", "Val", "Val", "Val", "Ser", "Ser", "Ser", "Ser", "Pro", "Pro", "Pro", "Pro", "Thr", "Thr", "Thr", "Thr", "Ala", "Ala", "Ala", "Ala", "Tyr", "Tyr", "His", "His", "Gln", "Gln", "Asn", "Asn", "Lys", "Lys", "Asp", "Asp", "Glu", "Glu", "Cys", "Cys", "Trp", "Arg", "Arg", "Arg", "Arg", "Ser", "Ser", "Arg", "Arg", "Gly", "Gly", "Gly", "Gly"]
    bases = ["A","T","C","G"]
    #create a list for storing all codons that have a fourfold degenerate site in their third position
    #loop over all the codons and store the first two bases as a string.
    #then loop over the four bases and append each of them in turn to the string containing the first two bases in the current codon.
    #check whether the amino-acid encoded for by the new codon is the same as that encoded for by the new codon you just made.
    #if this is the case for all four bases, add the current codon to the list you created in the beginning of this section.
    fourfold_deg_codons = []
    for i in range(len(codons)):
        first_two = "".join(codons[i][0:2])
        codon_change = False
        for j in bases:
            codon_comp = "".join([first_two,j])
            if aas[i] != aas[codons.index(codon_comp)]:
                codon_change = True
        if codon_change == False:
            fourfold_deg_codons.append(codons[i])
    return(fourfold_deg_codons)

def GC_along_sequence(sequence, length):
    '''
    Which sites in this sequence have GC?
    length: how far into the sequence should you go?
    '''
    GC = ["G", "C"]
    bases = ["A", "T", "C", "G"]
    sequence = sequence[:length]
    result = np.array([i for i in range(length)], dtype = "float")
    for pos, site in enumerate(sequence):
        if site in bases:
            if site in GC:
                result[pos] = 1
            else:
                result[pos] = 0
        else:
            result[pos] = np.nan
    return(result)

def GC4_along_sequence(sequence, length, remove_start = False, bases = None):
    '''
    Which fourfold degenerate sites in this sequence have GC?
    '''
    if not bases:
        bases = ["G", "C"]
    if length % 3 != 0:
        print("The required length has to be a multiple of 3!")
        raise Exception
    sequence = sequence[:length]
    if remove_start:
        sequence = sequence[3:]
    fourfold = get_4fold_deg(sequence)
    result = np.array([i for i in range(int(len(sequence)/3))], dtype = "float")
    third_pos = [i for i in range(2, len(sequence), 3)]
    for pos, site in enumerate(third_pos):
        if site in fourfold:
            if sequence[site] in bases:
                result[pos] = 1
            else:
                result[pos] = 0
        else:
            result[pos] = np.nan
    return(result)


def GC_along_sequence_set(names, sequences, length, fs = None):
    '''
    Calculate GC at the different sites in a set of sequences. Return also the number of data points.
    '''
    names = [names[pos] for pos, i in enumerate(names) if len(sequences[pos]) >= length]
    sequences = [i for i in sequences if len(i) >= length]
    GC_values = {}
    for i in range(len(names)):
        if i % 100 == 0:
            print(i)
        GC_values[names[i]] = GC_along_sequence(sequences[i], length)
    print(len(GC_values))
    if fs:
        GC_values = fs.average_over_families_2d(GC_values, remove_nans = True, nan_for_None = True)
    GC_values = np.array(list(GC_values.values()))
    result = np.nanmean(GC_values, axis = 0)
    n = GC_values.shape[0]
    return(result, n)

def GC4_along_sequence_set(names, sequences, length, remove_start = False, fs = None, require_a_member = True, bases = None):
    '''
    Calculate GC4 at the different sites in a set of sequences.
    '''
    if not bases:
        bases = ["G", "C"]
    names = [names[pos] for pos, i in enumerate(names) if len(sequences[pos]) >= length]
    sequences = [i for i in sequences if len(i) >= length]
    GC4_values = {}
    for i in range(len(names)):
        if i % 100 == 0:
            print(i)
        GC4_values[names[i]] = GC4_along_sequence(sequences[i], length, remove_start = remove_start, bases = bases)
    print(len(GC4_values))
    if fs:
        GC4_values = fs.average_over_families_2d(GC4_values, remove_nans = True, nan_for_None = True, require_a_member = require_a_member, remove_empty = True)
    GC4_values = np.array(list(GC4_values.values()))
    result = np.nanmean(GC4_values, axis = 0)
    n = GC4_values.shape[0]
    return(result, n)

def GC_average_along_sequence_set(names, sequences, length, fs = None):
    '''
    Calculate GC at the different sites in a set of sequences and average within codons. Return also the number of data points.
    '''
    if length % 3 != 0:
        print("The requested length has to be a multiple of 3!")
        raise Exception
    names = [names[pos] for pos, i in enumerate(names) if len(sequences[pos]) >= length]
    sequences = [i for i in sequences if len(i) >= length]
    GC_values = {}
    for i in range(len(names)):
        if i % 100 == 0:
            print(i)
        GC_values[names[i]] = GC_along_sequence(sequences[i], length)
    print(len(GC_values))
    if fs:
        GC_values = fs.average_over_families_2d(GC_values, remove_nans = True, nan_for_None = True)
    GC_values = np.array(list(GC_values.values()))
    result_temp = np.nanmean(GC_values, axis = 0)
    result = []
    for pos in range(0, len(result_temp), 3):
        current_mean = np.mean(result_temp[pos:pos + 3])
        result.append(current_mean)
    n = GC_values.shape[0]
    return(np.array(result), n)

def GC_corr_nonsyn_syn(names, sequences, length, fs = None, positions = None):
    '''
    Scan a bunch of CDSs and calculate the correlation between GC content at synonymous and GC content at non-synonymous sites.
    '''
    corrs = []
    if length % 3 != 0:
        print("Length has to be a multiple of three!")
        raise Exception
    map_to_binary = {"A": 0, "T": 0, "G": 1, "C": 1}
    for position in range(0, length, 3):
        print(position)
        current_bases = {}
        for site in range(3):
            current_bases[site] = {names[i]: map_to_binary[sequences[i][position + site]] if sequences[i][position + site] in map_to_binary else np.nan for i in range(len(sequences))}
            if fs:
                current_bases[site] = fs.average_over_families(current_bases[site], use_nan = True)
        if positions:
            ns = [current_bases[positions[0]][i] for i in sorted(current_bases[0].keys())]
            s = [current_bases[positions[1]][i] for i in sorted(current_bases[0].keys())]
        else:
            ns = [np.mean([current_bases[0][i], current_bases[1][i]]) for i in sorted(current_bases[0].keys())]
            s = [current_bases[2][i] for i in sorted(current_bases[0].keys())]
        ns_copy = ns.copy()
        ns = [ns[i] for i in range(len(ns)) if (not np.isnan(ns[i]) and not np.isnan(s[i]))]
        s = [s[i] for i in range(len(s)) if (not np.isnan(ns_copy[i]) and not np.isnan(s[i]))]
        corrs.append(scipy.stats.spearmanr(ns, s)[0])
    datapoint_number = len(current_bases[0])
    return(corrs, datapoint_number)
        

def generate_all_kmers(k):
    '''
    Generate all possible DNA k-mers of a particular k.
    '''
    bases = ["A", "T", "C", "G"]
    kmers = [''.join(i) for i in it.product(bases, repeat=k)]
    return(kmers)

def generate_random_sequence(length, bases, seed = None):
    '''
    Generate a random sequence of a specified length.
    '''
    if seed:
        random.seed(seed)
    sequence = [random.choice(bases) for i in range(length)]
    sequence = "".join(sequence)
    return(sequence)

def get_4fold_deg(sequence):
    '''
    Determine where the fourfold degenerate positions in a sequence are.
    '''
    fourfold_pos = []
    for i in range(0,len(sequence),3):
        if sequence[i:i+3] in fourfold_deg_list:
            fourfold_pos.append(i+2)
    return(fourfold_pos)

def get_ambiguous_nucl(nucl_list):
    '''
    Get the ambiguous IUPAC base corresponding to a list of DNA bases.
    '''
    query = sorted(nucl_list)
    query = "".join(query)
    ambiguous_nucl = _IUPAC_dict_[query]
    return(ambiguous_nucl)

def get_base_comp(bases, fasta, count = False):
    '''
    Get the over-all base frequencies in the sequences in a fasta file.
    '''
    global_counter = 0
    base_counters = {base: 0 for base in bases}
    with open(fasta) as file:
        for line in file:
            if line[0] == ">":
                print(line)
            else:
                for base in bases:
                    #sum up upper-case and lower-case occurrences
                    current_count_upper = line.count(base)
                    current_count_lower = line.count(base.lower())
                    current_count = current_count_lower + current_count_upper
                    base_counters[base] = base_counters[base] + current_count
                    #keeps track of the total number of bases so you could later convert counts to frequencies
                    #you can't just sum the lengths of the sequences because
                    #then you would also be counting N-bases
                    global_counter = global_counter + current_count
    if not count:
        for i in base_counters:
            base_counters[i] = base_counters[i]/global_counter
    return(base_counters)

def get_base_comp_motifs(bases, motifs):
    '''
    Get the over-all base composition of a set of motifs. Presumes that the motifs only contain bases supplied in the list 'bases'.
    '''
    motifs = "".join(motifs)
    length = len(motifs)
    base_dict = {base: (motifs.count(base))/length for base in bases}
    return(base_dict) 

def get_GC(sequence, alternative_bases = None):
    '''
    Get the GC content of a sequence. Can also be used for other sets of bases than GC.
    '''
    sequence = sequence.upper()
    sequence = [i for i in sequence if i!= "N"]
    if alternative_bases:
        bases = alternative_bases
    else:
        bases = ["G","C"]
    hits = 0
    for i in sequence:
        if i in bases:
            hits = hits + 1
    GC = hits/len(sequence)
    return(GC)

def get_GC4(sequence, phase, alternative_bases = None):
    '''
    Get the GC content at the fourfold degenerate bases of a sequence.
    '''
    #make sure the sequence starts and ends with full codons.
    sequence = trim_sequence(sequence, phase)
    #figure out where the 4-fold degenerate sites are.
    fourfold_deg_pos = get_4fold_deg(sequence)
    if not fourfold_deg_pos:
        return None
    if alternative_bases:
        bases = alternative_bases
    else:
        bases = ["G","C"]
    hits = 0
    for i in range(2,len(sequence),3):
        if sequence[i] in bases:
            if i in fourfold_deg_pos:
                hits = hits + 1
    GC4 = hits/len(fourfold_deg_pos)
    return(GC4)

def get_longest_run(motifs):
    '''
    Given a set of motifs, calculate how long is the longest mononucleotide
    run that you observe with each of the four DNA bases.
    '''
    max_motif_length = max([len(i) for i in motifs])
    bases_dict = {i: [] for i in ["A", "T", "C", "G"]}
    for base in bases_dict:
        for length in range(0, max_motif_length + 1):
            current_run = "".join([base for i in range(length)])
            for motif in motifs:
                if current_run in motif:
                    bases_dict[base].append(length)
        bases_dict[base] = max(bases_dict[base])
    return(bases_dict)

def get_motif_enrichment(motifs, fasta, sim_file_prefix, n_sim):
    '''
    Determine which motifs out of a set of supplied motifs are enriched
    in a particular fasta file when compared to control fastas.
    Concatenate the sequences within the fasta.
    Return a dictionary of Z-scores.
    '''
    names, sequences = rw.read_fasta(fasta)
    concat_sequence = "".join(sequences)
    concat_sequence = concat_sequence.upper()
    #compile a lookahead regex for each motif
    motif_regex = [re.compile("".join([i[0],"(?=",i[1:],")"])) for i in motifs]
    motif_lengths = [len(i) for i in motifs]
    real_densities = [calc_density_c(motif_regex[i], concat_sequence, motif_lengths[i])[0] for i in range(len(motifs))]
    sim_densities = np.zeros((len(motifs), n_sim))
    for i in range(n_sim):
        current_names, current_sequences = rw.read_fasta("".join([sim_file_prefix, str(i), ".fasta"]))
        current_concat_sequence = "".join(current_sequences)
        current_concat_sequence = current_concat_sequence.upper()
        sim_densities[:,i] = [calc_density_c(motif_regex[j], current_concat_sequence, motif_lengths[j])[0] for j in range(len(motifs))]
    Zscore_dict = {}
    for pos, i in enumerate(motifs):
        Zscore_dict[i] = ms.calc_Zscore(real_densities[pos], sim_densities[pos, :].copy())
    return(Zscore_dict)

def get_motif_enrichment_CLIP(motifs, fasta, sim_file_prefix, n_sim):
    '''
    Determine which motifs out of a set of supplied motifs are enriched
    in a particular fasta file when compared to control fastas.
    ***Don't concatenate*** the sequences within the fasta.
    Return a dictionary of enrichment p-values.
    '''
    #read in the sequences of interest
    names, sequences = rw.read_fasta(fasta)
    sequences = [i.upper() for i in sequences]
    #create a numpy array of booleans that will have one row per sequence in the fasta and one column per motif
    real_values = np.zeros((len(sequences), len(motifs)), dtype=bool)
    #put in a True or a False in each cell depending on whether that motif appears in that sequence
    for pos, i in enumerate(motifs):
        real_values[:, pos] = [i in j for j in sequences]
    sim_counts = np.zeros((len(motifs), n_sim))
    #do the same for each control sequence
    for i in range(n_sim):
        current_names, current_sequences = rw.read_fasta("".join([sim_file_prefix, str(i), ".fasta"]))
        current_sequences = [sequence.upper() for sequence in current_sequences]
        for pos, j in enumerate(motifs):
            current_sim_values = [j in k for k in current_sequences]
            #count the 'True's
            sim_counts[pos, i] = len([l for l in current_sim_values if l == True])
    enrichment_p_dict = {}
    for pos, i in enumerate(motifs):
        current_real_values = real_values[:, pos]
        #count the 'True's
        real_count = len(current_real_values[current_real_values])
        enrichment_p_dict[i] = ms.calc_eff_p(real_count, sim_counts[pos,], greater = True)
    return(enrichment_p_dict)

def get_motif_set_density(motifs, motif_lengths, sequence, n_sim = False, simulants = False, concat = False, trim = False, raw = False, dont_expand = False):
    '''
    Calculate the combined density of a set of motifs in a sequence.
    '''
    #determine where the motifs match in the sequence
    matches = calc_density_for_concat_several_c(motifs, sequence, motif_lengths)
    if raw:
        return(matches)
    matches = [j.flatten() for j in matches if len(j) > 0]
    seq_length = len(sequence)
    #uniquify the matches across all the different motifs and store either the count or the density
    if matches:
        matches = np.concatenate(tuple(matches))
        matches = np.unique(matches)
        count = len(matches)
        if not concat:
            density = count/seq_length
    else:
        count = 0
        density = 0
    #if you want to determine an enrichment p-value
    if n_sim:
        sim_densities = [0 for i in range(n_sim)]
        sim_matches = [[] for i in range(n_sim)]
        #if you want to return a density
        if not concat:
            #if you're shuffling codons rather than motifs
            if type(simulants) != list:
                for i in range(n_sim):
                    #calculate the density for each simulant
                    current_sequence = shuffle_sequence(sequence, units = simulants, trim = trim)
                    sim_matches[i] = np.unique(np.array(flatten([calc_density_c(motifs[j], current_sequence, motif_lengths[j], return_count = True)[1] for j in range(len(motifs))])))
                    sim_densities[i] = np.divide(len(sim_matches[i]), seq_length)  
            else:
                for i in range(n_sim):
                    current_motifs = simulants[i]
                    sim_matches[i] = np.unique(np.array(flatten([calc_density_c(current_motifs[j], sequence, motif_lengths[j], return_count = True)[1] for j in range(len(current_motifs))])))
                    sim_densities[i] = np.divide(len(sim_matches[i]), seq_length)  
            ND = ms.normalize(density, sim_densities)
            return({"density": density, "positions": matches, "simulated densities": sim_densities,
                    "simulated positions": sim_matches, "ND": ND})
        #if you want to return the number of bases that overlap
        else:
            if type(simulants) != list:
                for i in range(n_sim):
                    #calculate the count for each simulant
                    current_sequence = shuffle_sequence(sequence, units = simulants, trim = trim)
                    current_matches = calc_density_for_concat_several_c(motifs, current_sequence, motif_lengths)
                    current_matches = [j.flatten() for j in current_matches if len(j) > 0]
                    if current_matches:
                        current_matches = np.concatenate(tuple(current_matches))
                        sim_matches[i] = np.unique(current_matches)
                        sim_densities[i] = len(sim_matches[i])
                    else:
                        sim_densities[i] = 0
                        sim_matches[i] = current_matches

            else:
                for i in range(n_sim):
                    current_motifs = simulants[i]
                    current_matches = calc_density_for_concat_several_c(current_motifs, sequence, motif_lengths)
                    current_matches = [j.flatten() for j in current_matches if len(j) > 0]
                    if current_matches:
                        current_matches = np.concatenate(tuple(current_matches))
                        sim_matches[i] = np.unique(current_matches)
                        sim_densities[i] = len(sim_matches[i])
                    else:
                        sim_densities[i] = 0
                        sim_matches[i] = current_matches
            return({"count": count, "positions": matches, "simulated counts": sim_densities,
                    "simulated positions": sim_matches})
    else:
        if concat:
            return({"count": count, "positions": matches})
        else:
            return({"density": density, "positions": matches})

def get_neighbours(motifs, individual = False):
    '''
    For a set of motifs, generate all the unique motifs that are a single base substitution away.
    Don't allow motifs that are in the original set.
    individual: rather than generate all neighbouring motifs, randomly pick one neighbour for each of the
    motifs in the original set
    '''
    results = []
    bases = ["A", "T", "C", "G"]
    #because strings are annoying
    motifs_list = [list(i) for i in motifs]
    for motif in motifs_list:
        #this list will only be used in individual mode
        temp_results = []
        for site in range(len(motif)):
            #make everything now, including the motif itself, clean up later
            for base in bases:
                new_motif = motif.copy()
                new_motif[site] = base
                if not individual:
                    results.append("".join(new_motif))
                else:
                    temp_results.append("".join(new_motif))
        if individual:
            found = False
            while not found:
                picked_motif = "".join(random.choice(temp_results))
                if picked_motif not in results and picked_motif not in motifs:
                    found = True
                    results.append(picked_motif)
    #uniquify
    results = list(set(results))
    results = [i for i in results if i not in motifs]
    return(results)

def get_sequence_set_density(fasta_name, bed_name, motifs, simulants, n_sim, densities_output_file_name, sim_densities_output_file_name, positions_output_file_name, sim_positions_output_folder_name, feature_set = None, concat = False, positions = True, verbose = False, workers = None, trim = False, two_seqs = False):
    '''
    Determine the density of a set of motifs in the sequences in a fasta file.
    '''
    densities = {}
    sim_densities = {}
    NDs = {}
    lengths = {}
    fasta_names, fasta_seq = rw.read_fasta(fasta_name)
    if two_seqs:
        fasta_seq = [i.split("|") for i in fasta_seq]
        fasta_seq = [i[0] for i in fasta_seq]
        fasta_seq = ["".join([j for j in i if j in _canon_bases_]) for i in fasta_seq]
    bed_dict = {}
    #if, for instance, the sequence names in your fasta file refer to coordinates and not to transcript names,
    #you can supply a bed file that has one record per sequence in the fasta file (in the same order!) and
    #python will use the names in the bed file rather than in the fasta file to refer to the sequences
    #(can be necessary, for example, when averaging across families)
    if bed_name:
        bed = rw.read_many_fields(bed_name, "\t")
        for i in range(len(bed)):
            bed_dict[bed[i][3]] = bed[i]
            fasta_names[i] = bed[i][3]
    fasta_records = zip(fasta_names, fasta_seq)
    fasta_records = [list(i) for i in fasta_records]
    #calculate the real and simulated densities for all the different sequences in the fasta file (can also return counts instead of sequences when "concat" is True)
    trim_arg = {"trim": trim}
    results = run_in_parallel(fasta_records, ["foo", motifs, n_sim, simulants, feature_set, concat, positions, trim_arg, verbose], get_sequence_set_density_core, workers = None)
    with open(densities_output_file_name, "w") as dens_output, open(sim_densities_output_file_name, "w") as sim_dens_output, open(positions_output_file_name, "w") as pos_output:
        #store the real densities/counts, simulated densities/counts and the positions of the real hits in the various output files
        dens_output_writer = csv.writer(dens_output, delimiter = "\t")
        sim_dens_output_writer = csv.writer(sim_dens_output, delimiter = "\t")
        pos_output_writer = csv.writer(pos_output, delimiter = "\t")       
        for i in results:
            curr_densities, curr_sim_densities, curr_NDs_or_lengths, curr_real_pos_dict, curr_pos_by_simulation = i.get()
            for j in list(curr_densities.keys()):
                if n_sim:
                    dens_output_writer.writerow([j, curr_densities[j], curr_NDs_or_lengths[j]])
                    sim_dens_output_writer.writerow([j] + curr_sim_densities[j])
                try:
                    pos_output_writer.writerow(curr_real_pos_dict[j])
                except KeyError:
                    pass
                #also store the data in lists
                densities[j] = curr_densities[j]
                if n_sim:
                    sim_densities[j] = curr_sim_densities[j]
                if concat:
                    lengths[j] = curr_NDs_or_lengths[j]
                else:
                    if n_sim:
                        NDs[j] = curr_NDs_or_lengths[j]
    #average the statistics over families
    if feature_set:
        densities = feature_set.average_over_families(densities)
        NDs = feature_set.average_over_families(NDs)
        sim_densities = feature_set.average_over_families_2d(sim_densities)
        if concat:
            lengths = feature_set.average_over_families(lengths)
    #if you're meant to be returning a single point estimate for the whole set of sequences
    if concat:
        total_length = np.sum(np.array(list(lengths.values())))
        density = np.sum(np.array(list(densities.values())))/total_length
        if n_sim:
            sim_densities = [np.array(i) for i in list(sim_densities.values())]
            #sum the counts for each simulate and divide by the total length of the sequences
            densities_per_sim = [np.array([i[j] for i in sim_densities]) for j in range(n_sim)]
            densities_per_sim = [np.sum([j for j in i if j != None])/total_length for i in densities_per_sim]
            mean_sim = np.mean(densities_per_sim)
            ND = (density - mean_sim)/mean_sim
            sd_sim = np.std(densities_per_sim, ddof = 1)
            #probability that a density this high could have been produced by chance
            p = ms.calc_eff_p(density, densities_per_sim)
            #probability that a density this low could have been produced by chance
            depletion_p = ms.calc_eff_p(density, densities_per_sim, greater = False)
            Z = ms.calc_Zscore(density, densities_per_sim)
            output_dict = {"density": density, "ND": ND, "effective p": p, "depletion p": depletion_p, "simulated densities": densities_per_sim, "Z": Z, "simulant sd": sd_sim}
        else:
            output_dict = {"density": density}           
    else:
        median_density = np.median(np.array(list(densities.values())))
        if n_sim:
            NDs = [i for i in list(NDs.values()) if i!= None]
            median_ND = np.median(np.array(NDs))
            sim_densities = list(sim_densities.values())
            #take the median density obtained in each simulation
            medians_per_sim = [np.array([i[j] for i in sim_densities]) for j in range(n_sim)]
            medians_per_sim = [np.median(i) for i in medians_per_sim]
            sd_sim = np.std(medians_per_sim, ddof = 1)
            p = ms.calc_eff_p(median_density, medians_per_sim)
            depletion_p = ms.calc_eff_p(median_density, medians_per_sim, greater = False)
            Z = ms.calc_Zscore(median_density, medians_per_sim)
            output_dict = {"median density": median_density, "median ND": median_ND, "simulated densities": medians_per_sim, "effective p": p, "depletion p": depletion_p, "Z": Z, "simulant sd": sd_sim}
        else:
            output_dict = {"median density": median_density}           
    return(output_dict)

def get_sequence_set_density_core(fasta_records, motifs, n_sim, simulants, feature_set, concat, positions, trim = False, verbose = False):
    '''
    The bit that is parallelized in get_sequence_set_density.
    '''
    #to avoid circular import
    from bedtools_games import matches_to_bed
    motif_regex = [re.compile("".join([i[0],"(?=",i[1:],")"])) for i in motifs]
    motif_lengths = [len(i) for i in motifs]
    #check that you're not meant to be using shuffled sequences as control
    if type(simulants) == list:
        simulants = [[re.compile("".join([i[0],"(?=",i[1:],")"])) for i in j] for j in simulants]
    densities = {}
    sim_densities = {}
    NDs_or_lengths = {}
    real_pos_dict = {}
    if n_sim:
        pos_by_simulation = [[] for i in range(n_sim)]
    else:
        pos_by_simulation = []
    for pos, record in enumerate(fasta_records):
        if pos%100 == 0:
            if verbose:
                print(pos)
        current_name = record[0]
        sequence = record[1]
        if feature_set:
            #this is in case there is, say, "_full_CDS" tagged onto the sequence name, which would mess up the averaging over families
            current_name = current_name.split("_")
            current_name = current_name[0]
        if concat:
            current_dict = get_motif_set_density(motif_regex, motif_lengths, sequence, n_sim = n_sim, simulants = simulants, concat = True, trim = trim)
            curr_count = current_dict["count"]
            densities[current_name] = curr_count
            if n_sim:
                curr_sim_counts = current_dict["simulated counts"]
                sim_densities[current_name] = curr_sim_counts
            NDs_or_lengths[current_name] = len(sequence)
            if positions:
                if len(current_dict["positions"]) > 0:
                    real_positions = [int(i) for i in current_dict["positions"]]
                    real_pos_dict[current_name] = [current_name] + real_positions
                if n_sim:
                    curr_sim_positions = [[int(j) for j in i] for i in current_dict["simulated positions"]]
                    curr_sim_positions = [[current_name] + i if len(i) > 0 else None for i in curr_sim_positions]
                    for i in range(n_sim):
                        if curr_sim_positions[i] != None:
                            pos_by_simulation[i].append(curr_sim_positions[i])
        else:
            current_dict = get_motif_set_density(motif_regex, motif_lengths, sequence, n_sim = n_sim, simulants = simulants, trim = trim)
            curr_density = current_dict["density"]
            densities[current_name] = curr_density
            if n_sim:
                curr_ND = current_dict["ND"]
                NDs_or_lengths[current_name] = curr_ND
                curr_sim_densities = current_dict["simulated densities"]
                sim_densities[current_name] = curr_sim_densities
                temp_sim_densities = [current_name]
                temp_sim_densities.extend(curr_sim_densities)
            if positions:
                if len(current_dict["positions"]) > 0:
                    real_positions = [int(i) for i in current_dict["positions"]]
                    real_pos_dict[current_name] = [current_name] + real_positions
                if n_sim:
                    curr_sim_positions = [[int(j) for j in i] for i in current_dict["simulated positions"]]
                    curr_sim_positions = [[current_name] + i if len(i) > 0 else None for i in curr_sim_positions]
                    for i in range(n_sim):
                        if curr_sim_positions[i] != None:
                            pos_by_simulation[i].append(curr_sim_positions[i])
    return(densities, sim_densities, NDs_or_lengths, real_pos_dict, pos_by_simulation)

def index_to_codon(index):
    '''
    Given a 0-based index, return the 0-based indices of all the bases in the codon.
    '''
    if index%3 == 0:
        return(list(range(index, index+3)))
    elif index%3 == 1:
        return(list(range(index - 1, index+2)))
    else:
        return(list(range(index - 2, index + 1)))

def Ke_delta(sequence, Ke_dict, mutation, sim = False):
    '''
    Given a sequence and a mutation, calculate the impact of the mutation on Ke score.
    '''
    position = mutation[0]
    true_base = mutation[1]
    var_base = mutation[2]
    seq_list = list(sequence)
    delta = Ke_delta_core(sequence, seq_list, true_base, var_base, position, Ke_dict)
    if sim:
        sim_deltas = []
        base_pos = [i for i in range(len(sequence)) if sequence[i] == sequence[position]]
        base_pos = [i for i in base_pos if i != position] 
        for sim_pos in base_pos:
            sim_delta = Ke_delta_core(sequence, seq_list, true_base, var_base, sim_pos, Ke_dict)
            sim_deltas.append(sim_delta)
        output = {}
        output["delta"] = delta
        output["Z"] = ms.calc_Zscore(delta, sim_deltas)
        output["p"] = ms.calc_eff_p(delta, sim_deltas, greater = False)
        output["sim. number"] = len(sim_deltas)
    else:
        output = delta
    return(output)

def Ke_delta_core(sequence, seq_list, true_base, var_base, position, Ke_dict):
    '''
    Core for Ke_delta.
    '''
    true_score = Ke_score(sequence, Ke_dict, position)
    if np.isnan(true_score):
        return(None)
    mut_seq = seq_list.copy()
    if "X" in var_base:
        del_length = len(var_base)
        danger_pos = [i for i in range(position, position + del_length)]
        mut_seq = [i for pos, i in enumerate(mut_seq) if pos not in danger_pos]
    elif true_base == "I":
        mut_seq = [mut_seq[: (position + 1)]] + list(var_base) + [mut_seq[(position + 1): ]]
        mut_seq = flatten(mut_seq)
    else:
        mut_seq[position] = var_base
    mut_seq = "".join(mut_seq)
    mut_score = Ke_score(mut_seq, Ke_dict, position)
    if np.isnan(mut_score):
        return(None)
    delta = mut_score - true_score
    return(delta)

def Ke_score(sequence, Ke_dict, position):
    '''
    Given a sequence, a dictionary of motif scores and a position in the sequence,
    calculate the average Ke motif score at that position.
    '''
    fragment = sequence[(position - 6 + 1): (position + 6)]
    if not fragment:
        fragment = sequence[: (position + 6)]
    scores = []
    for index in range(0, (len(fragment) - 6 + 1)):
        current_motif = fragment[index: index + 6]
        if current_motif in Ke_dict:
            scores.append(Ke_dict[current_motif])
    if not scores:
        return(0)
    result = np.nanmean(scores)
    return(result)

def kmers_from_nc(k, n, comp_dict = None, genome_comp = False, return_freqs = False, seed = None, verbose = True):
    '''
    Generate n k-mers that would be sampled from either the hg38 mononucleotide
    composition or from the base composition in comp_dict.
    '''
    if seed:
        np.random.seed(seed = seed)
    if genome_comp:
        comp_dict = _genome_comp_dict_.copy()
    #numpy random.choice wants the probabilities in an array
    probs_array = np.array([comp_dict[i] for i in _canon_bases_])
    motifs = []
    for i in range(n):
        current_motif = np.random.choice(_canon_bases_, size = k, replace = True, p = probs_array)
        motifs.append("".join(current_motif))
    #check the mononucleotide composition of the motifs you generated
    flat_motifs = "".join(motifs)
    obtained_freqs = [(flat_motifs.count(i))/len(flat_motifs) for i in _canon_bases_]
    obtained_dict = {}
    for pos, i in enumerate(_canon_bases_):
        obtained_dict[i] = obtained_freqs[pos]
    if verbose:
        print("Generated {0} {1}-mers.\nInput frequencies:\n{2}\nOutput frequencies:\n{3}".format(len(motifs), len(motifs[0]), comp_dict, obtained_dict))
    if return_freqs:
        return(motifs, obtained_dict)
    else:
        return(motifs)
    
def make_simulants(motifs, n_sim, output_file_name = None, retro = False, mono = False, no_duplicates = False, remove_stops = False, remove_existing = False, cap_runs = False, exclude = None, seed = None, concat = True):
    '''
    Given a set of motifs, generate n_sim sets of simulants with the same over-all dinucleotide frequencies.
    motifs: a list of strings (the motifs for which you want to generate simulants)
    n_sim: number of simulant sets required
    mono: use mono- rather than dinucleotide frequencies
    no_duplicates: don't allow duplicates within simulant sets
    remove_stops: don't allow simulants that contain stop codons
    remove_existing: don't allow simulants to coincide with motifs in the true set
    cap_runs: don't allow mononucleotide runs whose length exceeds that of the longst run of that base in the true motifs
    exclude: don't allow any of the motifs within this list
    seed: supply a seed for the PRNG
    concat: concatenate the true motifs before extracting dinucleotides
    '''
    if cap_runs:
        longest_runs_numbers = get_longest_run(motifs)
        longest_runs_strings = ["".join([base for i in range(longest_runs_numbers[base] + 1)]) for base in longest_runs_numbers]
    motifs = [list(i) for i in motifs]
    nts = flatten(motifs)
    if mono:
        dints = flatten(motifs)
    else:
        if concat:
            dints = []
            #grab all the dinucleotides in the two reading frames
            for i in range(0, len(nts) - 1, 2):
                dints.append([nts[i], nts[i+1]])
            for i in range(1, len(nts) - 1, 2):
                dints.append([nts[i], nts[i+1]])
        else:
            dints = []
            for motif in motifs:
                for i in range(0, len(motif) - 1, 2):
                    dints.append([motif[i], motif[i+1]])
                for i in range(1, len(motif) - 1, 2):
                    dints.append([motif[i], motif[i+1]])
            
    #right now you have a list of lists. Turn that into a list of strings where each string is one dinucleotide.
    dints = ["".join(i) for i in dints]
    simulants = [["" for j in motifs] for i in range(n_sim)]
    for i, j in enumerate(motifs):
        if mono:
            dint_number = len(j)
        else:
            dint_number = len(j) // 2
        if len(j) % 2 == 0:
            even = True
        else:
            even = False
        for k in range(n_sim):
            found = False
            while not found:
                problem = False
                sim_motif = []
                for l in range(dint_number):
                    if seed:
                        random.seed(a = seed)
                        seed = seed + 1
                    sim_motif.append(random.choice(dints))
                #if the length of the motif is not an even number, add on a random mononucleotide from a bag made previously.
                if (not even) and (not mono):
                    sim_motif.append(random.choice(nts))
                sim_motif = "".join(sim_motif)
                if remove_stops:
                    if "TAA" in sim_motif or "TAG" in sim_motif or "TGA" in sim_motif:
                        problem = True
                if remove_existing:
                    if list(sim_motif) in motifs:
                        problem = True
                if cap_runs:
                    for run in longest_runs_strings:
                        if run in sim_motif:
                            problem = True
                if no_duplicates:
                    if sim_motif in simulants[k]:
                        problem = True
                if exclude:
                    if sim_motif in exclude:
                        problem = True
                if not problem:
                    found = True
                    simulants[k][i] = sim_motif
    if output_file_name:
        if not retro:
            with open(output_file_name, "w") as file:
                for i in range(n_sim):
                    file.write("|".join(simulants[i]))
                    file.write("\n")
        else:
            for i in range(n_sim):
                with open("{0}/fake_ESEs_{1}.txt".format(output_file_name, i), "w") as file:
                    file.write("\n".join(simulants[i]))
                    file.write("\n")
    return(simulants)

def make_simulants_markov(motifs, n_sim, output_file_name = None, remove_stops = False, remove_existing = False, noise = False):
    '''
    Given a set of motifs, make n_sim sets of simulants with the simulants sampled from the
    dinucleotide composition of the original set, using a Markov model to generate the simulants.
    '''
    #get mononucleotide frequencies
    monont_freqs = calc_nt_freqs(motifs, 1)
    monont_freqs_bases = sorted(list(monont_freqs.keys()))
    #convert the mononucleotide frequencies into a list
    monont_freqs_values = [monont_freqs[i] for i in monont_freqs_bases]
    #get dinucleotide frequencies
    dint_freqs = calc_nt_freqs(motifs, 2)
    #make a dictionary where the keys are bases and the values are lists of transition probabilities.
    #the order within each list is the same as that within the mononucleotide frequencies list
    transition_dict = nt_freqs_to_transitions(dint_freqs)
    transition_dict_list_form = {base: [transition_dict[base][base2] for base2 in monont_freqs_bases] for base in transition_dict}
    if noise:
        transition_dict_list_form =  {base: transition_dict_list_form[base] + np.random.uniform(0, 0.01, 4) for base in transition_dict_list_form}
        transition_dict_list_form =  {base: transition_dict_list_form[base]/np.sum(transition_dict_list_form[base]) for base in transition_dict_list_form}   
    simulants = []
    longest_runs_numbers = get_longest_run(motifs)
    longest_runs_strings = ["".join([base for i in range(longest_runs_numbers[base] + 1)]) for base in longest_runs_numbers]
    for sim in range(n_sim):
        if sim%100 == 0:
            print(sim)
        current_simulants = []
        for motif in motifs:
            found = False
            counter = 0
            while not found:
                counter = counter + 1
                problem = False
                #sample the first base based on mononucleotide frequencies
                current_motif = []
                try:
                    current_motif.append(np.random.choice(monont_freqs_bases, size = 1, p = monont_freqs_values)[0])
                except ValueError:
                    print(monont_freqs_values)
                    raise Exception
                #then, for all the other bases in the motif:
                for base in range(1, len(motif)):
                    previous_base = current_motif[-1]
                    #np.random.choice will return an array of bases so you have to get the first one
                    current_motif.append(np.random.choice(monont_freqs_bases, size = 1, p = transition_dict_list_form[previous_base])[0])
                current_motif = "".join(current_motif)
                if remove_stops:
                    #if there is a stop in the motif you just made, leave _found_ as False and go for another round
                    if "TAA" in current_motif or "TGA" in current_motif or "TAG" in current_motif:
                        problem = True
                if remove_existing:
                    if current_motif in motifs:
                        problem = True
                for run in longest_runs_strings:
                    if run in current_motif:
                        problem = True
                if (not problem) or (counter == 100):
                    found = True
                    current_simulants.append(current_motif)
        simulants.append(current_simulants)
    simulants_pooled = flatten(simulants)
    new_monont_freqs = calc_nt_freqs(simulants_pooled, 1)
    new_dint_freqs = calc_nt_freqs(simulants_pooled, 2)
    print("Mononucleotide frequencies in original set:")
    print(monont_freqs)
    print("Mononucleotide frequencies in simulants:")
    print(new_monont_freqs)
    print("Dinucleotide frequencies in original set:")
    print(dint_freqs)
    print("Dinucleotide frequencies in simulants:")
    print(new_dint_freqs)
    return(simulants)

def make_simulants_within(motifs, n_sim):
    '''
    Given a set of motifs, make n_sim sets of simulants with each simulant motif sampled from the
    dinucleotide composition of a particular real motif, not that of the whole set.
    '''
    even = False
    simulants = [[] for i in range(n_sim)]
    for motif in motifs:
        if len(motif) % 2 == 0:
            even = True
        motif_length = len(motif)
        #get the dinucleotides in the motif
        dints_phase0 = [motif[i:i+2] for i in range(0, len(motif) - 1, 2)]
        dints_phase1 = [motif[i:i+2] for i in range(1, len(motif) - 1, 2)]
        dints = dints_phase0 + dints_phase1
        for sim in range(n_sim):
            found = False
            while not found:
                new_motif = [random.choice(dints) for i in range(int(motif_length/2))]
                if not even:
                    new_motif.append(random.choice(motif))
                new_motif = "".join(new_motif)
                problem = False
                if not problem:
                    found = True                      
                    simulants[sim].append(new_motif)
    return(simulants)

def map_ranks_to_counts(hit_base_ranks, monont_pos):
    '''
    Given the ranks of bases in true motif hits and the numbers
    of different bases available in the rest of the sequence,
    construct an initial guess for the base counts to use for control.
    '''
    counts = []
    for key in sorted(list(monont_pos.keys())):
        rank = hit_base_ranks[key]
        if rank == 3:
            counts.append(len(monont_pos[key]))
        elif rank == 0:
            counts.append(0)
        elif rank == 2:
            counts.append(round((len(monont_pos[key])/4 * 3)))
        elif rank == 1:
            counts.append(round(len(monont_pos[key])/4))
        else:
            print("Invalid rank!")
            print(hit_base_ranks)
            raise Exception
    return(counts)

def map_ranks_to_counts_brute(hit_base_ranks, monont_pos, hit_mononts):
    '''
    Given the ranks of bases in true motif hits, the numbers
    of different bases available in the rest of the sequence and
    the frequencies of bases in true hits,
    construct the base counts to use for control.
    '''
    counts = []
    top_base = [i for i in hit_base_ranks if hit_base_ranks[i] == 3][0]
    max_freq = hit_mononts[top_base]
    max_number = len(monont_pos[top_base])
    for key in sorted(list(monont_pos.keys())):
        rank = hit_base_ranks[key]
        current_freq = hit_mononts[key]
        current_number = len(monont_pos[key])
        if rank == 3:
            counts.append(current_number)
        else:
            expected = round((current_freq * max_number)/max_freq)
            if expected > current_number:
                counts.append(current_number)
            else:
                counts.append(expected)
    return(counts)

def mask_motifs(fasta, motifs, output_fasta):
    '''
    Given a fasta and a set of motifs, replace all the hits with X's.
    '''
    motif_lengths = [len(i) for i in motifs]
    motifs = motif_to_regex(motifs)
    names, seqs = rw.read_fasta(fasta)
    new_seqs = []
    for pos, seq in enumerate(seqs):
        if pos % 1000 == 0:
            print(pos)
        #1st element of tuple because you want just the positions, not the length of the overlap
        matches = get_motif_set_density(motifs, motif_lengths, seq, concat = True)["positions"]
        new_seq = ["X" if pos in matches else i for pos, i in enumerate(seq)]
        new_seq = "".join(new_seq)
        new_seqs.append(new_seq)
    rw.write_to_fasta(names, new_seqs, output_fasta)
    print("Wrote masked sequence to file.")

def maxent(sequence, site):
    '''
    Calculate the MaxEntScan score of a sequence.
    '''
    temp_file_name = "temp_data/temp_file{0}.txt".format(random.random())
    with open(temp_file_name, "w") as file:
        file.write(sequence)
    results = run_process(["perl", "score{0}.pl".format(site), temp_file_name])
    results = results.split("\t")
    score = float(results[1])
    remove_file(temp_file_name)
    return(score)

def mono_runs(motifs):
    '''
    Calculate what fraction of the motifs have mononucleotide runs of three or more bases.
    '''
    counter = 0
    for motif in motifs:
        if "AAA" in motif or "TTT" in motif or "CCC" in motif or "GGG" in motif:
            counter = counter + 1
    return(counter/len(motifs))

def mono_runs_core(sequence, regex):
    '''
    Core for the mono_runs_sequence function.
    '''
    matches = []
    all_hits = [re.finditer(i, sequence) for i in regex]
    for base_hit in all_hits:
        hit_pos = [[i for i in range(hit.start(), hit.end())] for hit in base_hit]
        [matches.extend(i) for i in hit_pos]
    try:
        current_fraction = len(matches)/len([i for i in sequence if i in _canon_bases_])
    except ZeroDivisionError:
        current_fraction = None
    return(current_fraction)
    

def mono_runs_sequence(names, sequences, min_run, n_sim = False, fs = None, output_file = None):
    '''
    For each sequence in _sequences_, calculate the fraction composed of mononucleotide runs
    of _min_run_ or more bases. If _n_sim_, shuffle each sequence _n_sim_ times and get an
    empirical distribution.
    '''
    results = {}
    if n_sim:
        results["real fractions"] = {}
        results["sim. fractions"] = {}
        results["norm. fractions"] = {}
    min_run = str(min_run)
    regex = [re.compile(base + "{" + min_run + ",}") for base in _canon_bases_]
    for pos, name in enumerate(names):
        if pos % 100 == 0:
            print(pos)
        current_sequence = sequences[pos]
        current_fraction = mono_runs_core(current_sequence, regex)
        if not n_sim:
            results[name] = current_fraction
        else:
            results["real fractions"][name] = current_fraction
            results["sim. fractions"][name] = []
            sequence_list = list(current_sequence)
            for sim in range(n_sim):
                np.random.seed(sim)
                new_seq = sequence_list.copy()
                np.random.shuffle(new_seq)
                new_seq = "".join(new_seq)
                current_sim_fraction = mono_runs_core(new_seq, regex)
                results["sim. fractions"][name].append(current_sim_fraction)
            current_norm_fraction = ms.normalize(current_fraction, results["sim. fractions"][name])
            results["norm. fractions"][name] = current_norm_fraction
    if not n_sim:
        if fs:
            results = fs.average_over_families(results)
        if output_file:
            with open(output_file, "w") as file:
                file.write("gene\treal fraction\n")
                for gene in sorted(list(results.keys())):
                    file.write("{0}\t{1}\n".format(gene, results[gene]))
            results_output = np.median(list(results.values()))
            results = results_output
    else:
        if fs:
            results["real fractions"] = fs.average_over_families(results["real fractions"])
            results["sim. fractions"] = fs.average_over_families_2d(results["sim. fractions"])
            results["norm. fractions"] = fs.average_over_families(results["norm. fractions"])
        if output_file:
            with open(output_file, "w") as file:
                file.write("gene\treal fraction\tmean sim. fraction\tnorm. fraction\n")
                sim_fractions = [[] for i in range(n_sim)]
                for gene in sorted(list(results["real fractions"].keys())):
                    mean_sim_fraction = np.mean(results["sim. fractions"][gene])
                    norm_fraction = results["norm. fractions"][gene]
                    file.write("{0}\t{1}\t{2}\t{3}\n".format(gene, results["real fractions"][gene], mean_sim_fraction, norm_fraction))
                    for sim in range(n_sim):
                        sim_fractions[sim].append(results["sim. fractions"][gene][sim])
                sim_fractions = [np.median(i) for i in sim_fractions]
                results_output = {}
                results_output["real median"] = np.median(list(results["real fractions"].values()))
                results_output["mean simulated median"] = np.mean(sim_fractions)
                results_output["normalized fraction"] = ms.normalize(results_output["real median"], sim_fractions)
                results_output["p"] = ms.calc_eff_p(results_output["real median"], sim_fractions, greater = False)
                results = results_output
    return(results)

def motif_kmers(motifs, k):
    '''
    Given  a set of motifs, return all the unique k-mers in the set.
    '''
    freqs_dict = calc_nt_freqs(motifs, k)
    kmers = []
    for key in freqs_dict:
        if freqs_dict[key] > 0:
            kmers.append(key)
    return(kmers)

def motif_phase(motif, sequence):
    '''
    Given a motif and a sequence, count how many of the hits are in each phase.
    '''
    positions = [hit.start() for hit in motif.finditer(sequence)]
    phase_dict = {}
    for phase in range(3):
        hits = [i for i in positions if i%3 == phase]
        phase_dict[phase] = len(hits)
    return(phase_dict)

def motif_phase_deviation(motifs, sequences, return_chi = False):
    '''
    Given a set of motifs and a set of sequences, calculate to what extent the phase distribution of the hits deviates from a uniform distribution.
    '''
    phase_counts = motif_phase_sequence_set(motifs, sequences)
    output = []
    chis = []
    for motif in phase_counts:
        total = sum(list(phase_counts[motif].values()))
        expected = total/3
        deviations = [abs(i - expected)/expected if expected else np.nan for i in list(phase_counts[motif].values())]
        output.append(np.nanmean(deviations))
        if return_chi:
            if expected:
                chi = sum([((i - expected)**2)/expected for i in list(phase_counts[motif].values())])
            else:
                chi = np.nan
            chis.append(chi)
    if return_chi:
        return([np.nanmean(output), np.nanmean(chis)])
    return(np.nanmean(output))

def motif_phase_motif_set(motifs, motif_regex, sequence):
    '''
    Given a set of motifs and a sequence, count how many of the hits are in each phase for each motif.
    '''
    output = {}
    for pos, regex in enumerate(motif_regex):
        output[motifs[pos]] = motif_phase(regex, sequence)
    return(output)

def motif_phase_sequence_set(motifs, sequences):
    '''
    Given a set of motifs and a set of sequences, count how many of the hits are in each phase for each motif, summing across sequences.
    '''
    final_dict = {j: {i: 0 for i in range(3)} for j in motifs}
    motif_regex = motif_to_regex(motifs)
    for sequence in sequences:
        current_dict = motif_phase_motif_set(motifs, motif_regex, sequence)
        for motif in motifs:
            for phase in range(3):
                final_dict[motif][phase] = final_dict[motif][phase] + current_dict[motif][phase]
    return(final_dict)

def motif_to_regex(motifs):
    '''
    Convert a string into a lookahead regex where only the first base
    is matched and the rest is in the lookahead.
    '''
    regex = [re.compile("".join([i[0],"(?=",i[1:],")"])) for i in motifs]
    return(regex)

def motifs_from_peaks(bed, genome, k, n, output_file = None, fasta = None, exclude_overlaps = False):
    '''
    Given a bed file, use RSAT to determine the n k-mers that are most over-represented.
    '''
    from bedtools_games import fasta_from_intervals
    if not fasta:
        fasta = "temp_data/temp_fasta_{0}.fasta".format(random.random())
        fasta_from_intervals(bed, fasta, genome, force_strand = True)
    if not output_file:
        output_file = "temp_data/rsat_out.txt"
    arguments = ["oligo-analysis", "-i", fasta, "-format", "fasta", "-seqtype", "dna", "-o", output_file, "-l",
                k, "-markov", 1, "-return", "occ,proba,rank", "-sort", "-1str", "-nogrouprc"]
    if exclude_overlaps:
        arguments.append("-noov")
    run_process(arguments)
    data = rw.read_many_fields(output_file, "\t")
    data = data[1: (n + 1)]
    data = [i[0].upper() for i in data]
    if not fasta:
        os.remove(fasta)
    if not output_file:
        os.remove(output_file)
    return(data)

def motifs_to_file_retro(motifs, motifs_file_name):
    '''
    Write a list of motifs into an old school motifs file.
    '''
    with open("{0}.csv".format(motifs_file_name), "w") as file:
        file.write("{0}\n".format(motifs_file_name))
        for i in motifs:
            file.write("{0}\n".format(i))

def nt_freqs_to_transitions(frequencies):
    '''
    Convert a dictionary of dinucleotide frequencies into a dictionary of transition probabilities from one base to the other.
    '''
    transition_dict = {}
    for base in _canon_bases_:
        transition_dict[base] = {}
        temp_dict = {}
        for base2 in _canon_bases_:
            current_freq = frequencies["".join([base, base2])]
            temp_dict[base2] = current_freq
        total_fraction = sum(list(temp_dict.values()))
        if total_fraction == 0:
            transition_dict[base] = {base2: 0.25 for base2 in _canon_bases_}
        else:
            transition_dict[base] = {base2: temp_dict[base2]/total_fraction for base2 in _canon_bases_}
    return(transition_dict)

def nucleotide_comp_skew(motifs):
    '''
    Given a set of motifs, calculate to what extent their base composition deviates from a uniform distribution.
    '''
    skew_sum = 0
    freqs = calc_nt_freqs(motifs, 1)
    for base in _canon_bases_:
        current_skew = abs(freqs[base] - 0.25)
        skew_sum = current_skew + skew_sum
    skew = skew_sum/4
    return(skew)

def nucleotide_comp_skew_sequence(fasta, window_size = None):
    '''
    Given a set of sequences, calculate to what extent their base composition deviates from a uniform distribution.
    '''
    output_dict = {}
    if not window_size:
        freqs = calc_nt_freqs_sequence(fasta, 1)
        for name in freqs:
            skew_sum = 0
            for base in _canon_bases_:
                current_skew = abs(freqs[name][base] - 0.25)
                skew_sum = current_skew + skew_sum
            output_dict[name] = skew_sum/4
    else:
        names, seqs = rw.read_fasta(fasta)
        for pos, name in enumerate(names):
            if pos%1000 == 0:
                print(pos)
            total_skew = []
            current_seq = seqs[names.index(name)]
            for index in range(len(current_seq)):
                current_kmer = current_seq[index: index + window_size]
                if (len(current_kmer) == window_size) and ("X" not in current_kmer):
                    skew_sum = 0
                    current_freqs = calc_nt_freqs([current_kmer], 1)
                    if current_freqs:
                        for base in _canon_bases_:
                            current_skew = abs(current_freqs[base] - 0.25)
                            skew_sum = current_skew + skew_sum
                        total_skew.append(skew_sum/4)
            if len(total_skew) > 0:
                output_dict[name] = np.mean(total_skew)
            else:
                output_dict[name] = None
    return(output_dict)

def overlaps_within_motif_set(motifs):
    '''
    Calculate the extent to which there is overlaps within the set of motifs.
    '''
    starts = list(set([motif[:3] for motif in motifs]))
    ends = list(set([motif[-3:] for motif in motifs]))
    overlaps = len(overlap(starts, ends))
    return(overlaps/len(motifs))

def palindromes(motifs):
    '''
    Only counting the first 3 and the last 3 bases, what fraction of the motifs in the set are palindromes?
    '''
    counter = 0
    for motif in motifs:
        start = motif[:3]
        end = motif[-3:]
        if start == end[::-1]:
            counter = counter + 1
    return(counter/len(motifs))
        
def parse_retro_motif_density_output(dataset_name, output_suffix):
    '''
    Parse the output from find_ESEs_wrapper.py into a dictionary.
    '''
    folder_name = "{0}_folder".format(dataset_name)
    sim_densities = rw.read_numbers("{0}/{1}_ESE_density_per_sim{2}.csv".format(folder_name, dataset_name, output_suffix))
    sim_densities = [float(i) for i in sim_densities]
    densities = rw.read_many_fields("{0}/{1}_ESE_density{2}.csv".format(folder_name, dataset_name, output_suffix), ",")
    real_densities = [float(i[1]) for i in densities]
    sim_dens_per_gene = [float(i[2]) for i in densities]
    median_density = np.median(real_densities)
    NDs = [(real_densities[i] - sim_dens_per_gene[i])/sim_dens_per_gene[i] for i in range(len(real_densities))]
    median_ND = np.median(NDs)
    p = ms.calc_eff_p(median_density, sim_densities)
    output = {}
    output["median density"] = median_density
    output["median ND"] = median_ND
    output["p"] = p
    return(output)

def parse_sim_densities(lengths_dict, input_file_name, fs = None):
    '''
    Parses an RBP_density.py style sim densities file from a run with concatenation.
    '''
    data = rw.read_many_fields(input_file_name, "\t")
    data = {i[0]: [int(j) for j in i[1:]] for i in data}
    if fs:
        data = fs.average_over_families_2d(data)
        lengths_dict = fs.average_over_families(lengths_dict)
    sorted_keys = sorted(list(lengths_dict.keys()))
    data = np.array([data[i] for i in sorted_keys])
    data_colsums = np.sum(data, axis = 0)
    total_length = np.sum(np.array([lengths_dict[i] for i in sorted_keys]))
    densities = np.divide(data_colsums, total_length)
    return(densities)

def pick_control_pos(counts, monont_pos):
    '''
    Pick control positions from among a set of positions.
    '''
    counts = [round(i) for i in counts]
    result = {}
    new_sim_hits = []
    for pos, monont in enumerate(sorted(monont_pos.keys())):
        if len(monont_pos[monont]) > 0:
            new_pos = np.random.choice(monont_pos[monont], size = counts[pos], replace = False)
            new_sim_hits.extend(new_pos)
    new_sim_hits = sorted(list(set(new_sim_hits)))        
    return(new_sim_hits)

def pick_control_pos_codons(hit_file, fasta, control_file, human = False, ancestral = None, mono = False, trint = False, exclude_file = None):
    '''
    Given a hit file, pick control positions so as to match the codon composition (across all the different CDSs).
    '''
    codon_pos_file = "temp_data/codon_pos{0}.txt".format(random.random())
    codon_positions(fasta, codon_pos_file, hit_file, human = human, ancestral = ancestral, mono = mono, trint = trint, exclude_file = exclude_file)
    codon_pos_dict = rw.parse_codon_positions(codon_pos_file)
    remove_file(codon_pos_file)
    names, seqs = rw.read_fasta(fasta)
    hits = rw.read_pos(hit_file)
    output_codons = count_codons(names, seqs, hits, mono = mono, trint = trint)
    pos_list = []
    for codon in output_codons:
        #because np.random.choice doesn't like that it's a list of tuples
        dummy = [i for i in range(0, len(codon_pos_dict[codon]))]
        pos_list.extend([codon_pos_dict[codon][i] for i in np.random.choice(dummy, output_codons[codon], replace = True)])
    rw.write_codon_control_pos(pos_list, control_file, hits)

def prepare_control_pos_for_hits(names, seqs, motifs, anc_CG, macaque_CG, tuples_mapping = None, brute_mapping = False, CG_gene_filter = False, family_seed = None, exclude_file = None, alt_anc_CGs = None, fs = None, nonsyn_hits = False, leave_CG = False, verbose_detailed = False, remove_ancestral_CpG = False, remove_macaque_CpG = False, pseudoCG = False, match_size = False, prone_sites = False):
    motif_lengths = [len(i) for i in motifs]
    motif_regex = motif_to_regex(motifs)

    if fs:
        try:
            picked = fs.pick_random_members(family_seed)
        except AttributeError:
            picked = names.copy()
    else:
        picked = names.copy()

    if pseudoCG:
        CG_2mers = ["CT", "AG"]
    else:
        CG_2mers = ["CG", "GC"]
    CG_lengths = [2, 2]
    CG_regex = motif_to_regex(CG_2mers)

    if prone_sites:
        prone_2mers = [".G", "C."]
        prone_lengths = [2, 2]
        prone_regex = motif_to_regex(prone_2mers)

    print("Preparing sequence data...")

    hit_mononts_dict = {}
    counts_dict = {}
    true_lengths_dict = {}
    monont_pos_dict = {}
    bounds_dict = {}
    true_pos_dict = {}

    if remove_ancestral_CpG or verbose_detailed or CG_gene_filter:
        picked = [i for i in picked if i in anc_CG]
    if remove_macaque_CpG or verbose_detailed:
        picked = [i for i in picked if i in macaque_CG]

    if exclude_file != "None" and exclude_file != None:
        exclude_pos = rw.read_pos(exclude_file)
    else:
        exclude_pos = {i: [] for i in picked}

    counter = 0
    for pos, seq in enumerate(seqs):
        name = names[pos]
        if name in picked:
            true_pos_dict[name] = None
            hit_mononts_dict[name] = None
            counts_dict[name] = None
            true_lengths_dict[name] = None
            monont_pos_dict[name] = None
            bounds_dict[name] = None
            fourfold_pos = get_4fold_deg(seq)
            CG_pos = get_motif_set_density(CG_regex, CG_lengths, seq, concat = True)["positions"]
            if prone_sites:
                prone_pos = get_motif_set_density(prone_regex, prone_lengths, seq, concat = True)["positions"]
            counter = counter + 1
            if counter % 100 == 0:
                print(counter)
            true_hits = get_motif_set_density(motif_regex, motif_lengths, seq, concat = True)["positions"]
            
            true_hits_copy = true_hits.copy()
            if not nonsyn_hits:
                true_hits = [i for i in true_hits if i in fourfold_pos]
            else:
                true_hits = [i for i in true_hits if i not in fourfold_pos]
            if verbose_detailed:
                for hit in true_hits:
                    try:
                        if ((hit in CG_pos) or (hit in macaque_CG)) and (hit not in alt_anc_CGs[1][name]):
                            print("\n")
                            order = [3, 5, 2, 7, 1, 6, 4, 0] 
                            prev = [tuples_mapping[name][hit - 1][i] for i in order]
                            cur = [tuples_mapping[name][hit][i] for i in order]
                            subs = [tuples_mapping[name][hit + 1][i] for i in order]
                            for base in [prev, cur, subs]:
                                print(" ".join((base[:4] + ["||"] + base[4:-1] + ["||"] + [base[-1]])))
                            print("***")
                            print("Human + macaque: {0}.".format(str((hit in CG_pos) or (hit in macaque_CG))))
                            print("Ancestor (JC69): {0}.".format(str(hit in anc_CG[name])))
                            print("Ancestor (comp): {0}.".format(str(hit in alt_anc_CGs[0][name])))
                            print("Ancestor (U2S): {0}.".format(str(hit in alt_anc_CGs[1][name])))
                    except KeyError:
                        pass
            if not leave_CG:
                true_hits = [i for i in true_hits if i not in CG_pos]
##            temp_seq = "".join([seq[i] for i in anc_CG[name] if i in true_hits])
##            freqs = calc_nt_freqs(temp_seq, 1)
##            print(freqs)
            if prone_sites:
                true_hits = [i for i in true_hits if i not in prone_pos]
            if remove_ancestral_CpG:
                true_hits = [i for i in true_hits if i not in anc_CG[name]]
            if remove_macaque_CpG:
                true_hits = [i for i in true_hits if i not in macaque_CG[name]]
            if true_hits:
                true_pos_dict[name] = true_hits
                hit_seq = "".join([seq[i] for i in true_hits])
                hit_mononts = calc_nt_freqs([hit_seq], 1)
                hit_mononts_dict[name] = hit_mononts
                true_length = len(hit_seq)
                true_lengths_dict[name] = true_length
                nts = [i for i in sorted(list(hit_mononts.keys())) if hit_mononts[i] > 0]
                sim_hits = [[pos for pos, i in enumerate(seq) if i == j] for j in nts]
                monont_pos = {i: [] for i in hit_mononts}
                found = False
                for pos, nt in enumerate(nts):
                    if not leave_CG:
                        monont_pos[nt] = list(set([i for i in sim_hits[pos] if i not in true_hits_copy and i not in CG_pos and i in fourfold_pos]))
                    else:
                        monont_pos[nt] = list(set([i for i in sim_hits[pos] if i not in true_hits_copy and i in fourfold_pos]))
                    if remove_ancestral_CpG:
                        monont_pos[nt] = [i for i in monont_pos[nt] if i not in anc_CG[name]]
                    if remove_macaque_CpG:
                        monont_pos[nt] = [i for i in monont_pos[nt] if i not in macaque_CG[name]]
                    if prone_sites:
                        monont_pos[nt] = [i for i in monont_pos[nt] if i not in prone_pos]
                    #note that exclude_pos only filters the controls and not the hits
                    monont_pos[nt] = [i for i in monont_pos[nt] if i not in exclude_pos[name]]
                    if len(monont_pos[nt]) > 0:
                        found = True
                if found:
                    monont_pos_dict[name] = monont_pos
                    hit_base_ranks = sorted(hit_mononts.items(), key = operator.itemgetter(1))
                    hit_base_ranks = {i[0]: pos for pos, i in enumerate(hit_base_ranks)}
                    if match_size:
                        counts = [hit_mononts[i] for i in sorted(list(hit_mononts.keys()))]
                    else:
                        if brute_mapping:
                            counts = map_ranks_to_counts_brute(hit_base_ranks, monont_pos, hit_mononts)
                        else:
                            counts = map_ranks_to_counts(hit_base_ranks, monont_pos)
                    counts_dict[name] = counts
                    bounds = [(0, len(monont_pos[i])) for i in sorted(list(monont_pos.keys()))]
                    bounds_dict[name] = bounds
    return(true_pos_dict, hit_mononts_dict, counts_dict, true_lengths_dict, monont_pos_dict, bounds_dict)

def reassign_bases(motifs):
    '''
    Swap the bases around in a set of motifs (replacing, for instance,
    all the As with Cs, all the Ts with As...).
    '''
    print("Old motifs:")
    print(motifs)
    bases = ["A", "T", "C", "G"]
    base_dict = {i: i for i in bases}
    new_bases = bases.copy()
    #keep on shuffling around the bases until you get a configuration where none of the bases have been
    #assigned to themselves
    while True:
        random.shuffle(new_bases)
        for pos, i in enumerate(base_dict):
            #you store it in lower case because otherwise the substitution in the motifs
            #won't work correctly. (For instance, imagine that you first turn all the As into Cs and then all the Cs into As.
            #If you didn't go through lower case, you would rewrite the substitutions performed in the first step
            #in the second step.)
            base_dict[i] = new_bases[pos].lower()
        identical = False
        for i in base_dict:
            if base_dict[i] == i.lower():
               identical = True
        if not identical:
            break
    #do the actual substitution in the motifs
    for i in range(len(motifs)):
        for j in base_dict:
            motifs[i] = re.sub(j, base_dict[j], motifs[i])
        motifs[i] = motifs[i].upper()
    print("New motifs:")
    print(motifs)
    print("\n")
    return(motifs)

def reassign_bases_all_possible(motifs):
    '''
    Perform all the reassignments possible for a set of motifs (see reassign_bases).
    '''
    bases = ["A", "T", "C", "G"]
    base_dict = {i: i for i in bases}
    all_new_bases = list(it.permutations(bases))
    final_motifs = []
    for base_set in all_new_bases:
        for pos, i in enumerate(base_dict):
            base_dict[i] = base_set[pos].lower()
        identical = False
        for i in base_dict:
            if base_dict[i] == i.lower():
               identical = True
        if not identical:
            temp_motifs = motifs.copy()
            for i in range(len(temp_motifs)):
                for j in base_dict:
                    temp_motifs[i] = re.sub(j, base_dict[j], temp_motifs[i])
                temp_motifs[i] = temp_motifs[i].upper()
            final_motifs.append(temp_motifs)
    return(final_motifs)

def remove_positions(seq, nts, positions):
    '''
    Given a sequence, a set of kmers with k between 1 and 3, and a set of indices for the sequence,
    remove all indices where the corresponding codon overlaps with a kmer from the set.
    '''
    lengths = [len(i) for i in nts]
    nts = motif_to_regex(nts)
    to_remove = []
    for pos, kmer in enumerate(nts):
        matches = re.finditer(kmer, seq)
        if matches:
            for match in matches:
                offending_pos = [i for i in range(match.start(), match.start() + lengths[pos] + 1)]
                to_remove.extend(offending_pos)
    clean_positions = []
    for codon in range(0, len(positions) - 1, 3):
        if (positions[codon] not in to_remove) and (positions[codon + 1] not in to_remove) and (positions[codon + 2] not in to_remove):
            clean_positions.extend(positions[codon: codon + 3])
    return(clean_positions)

def repeats(fasta, order):
    '''
    Given a set of sequences, what fraction of each is made up of mononucleotide repeats that are order or more long?
    '''
    output = {}
    names, seqs = rw.read_fasta(fasta)
    initial = "".join(["([ATCG])" for i in range(order)])
    expansion = "".join(["\{0}".format(i) for i in range(1, order + 1)])
    regex = initial + "(" + expansion + "){1,}"
    regex = re.compile(regex)
    for pos, seq in enumerate(seqs):
        current_hits = re.finditer(regex, seq)
        current_hits = [i.group(0) for i in current_hits]
        hit_length = np.sum([len(i) for i in current_hits])
        fraction = hit_length/len([i for i in seq if i in _canon_bases_])
        output[names[pos]] = fraction
    return(output)

def repeats_within_motif(motifs, order, immediate = False):
    '''
    In a set of motifs, what fraction of motifs contain repeated order-nucleotides?
    '''
    counter = 0
    for motif in motifs:
        found = False
        for base in range(len(motif)):
            current_string = motif[base: base + order]
            if len(current_string) == order:
                if immediate:
                    if motif[base + order: base + (order * 2)] == current_string:
                        found = True
                else:
                    remaining_sequence = motif[base + order:]
                    if current_string in remaining_sequence:
                        found = True
        if found:
            counter = counter + 1
    return(counter/len(motifs))

def rev_comp(base):
    '''
    Reverse complement a base.
    '''
    return(_rev_comp_dict_[base])

def shuffle_sequence(sequence, units, trim = False, seed = None):
    '''
    Shuffle a sequence by chunks of size units.
    '''
    if seed:
        random.seed(seed)
    remainder = len(sequence) % units
    if remainder != 0:
        if not trim:
            print("In order to shuffle chunks of size {0}, the length of the sequence must be a multiple of {0}!\n Set trim to True if sequences should be trimmed from the end to comply!".format(units))
            raise Exception
        else:
            sequence = sequence[:-remainder]            
    chunks_in_seq = [sequence[i:i+units] for i in range(0, len(sequence), units)]
    random.shuffle(chunks_in_seq)
    sequence = "".join(chunks_in_seq)
    return(sequence)

def split_clusters(clusters, threshold):
    '''
    Split clusters if their length exceeds a threshold.
    '''
    output_clusters = []
    for cluster in clusters:
        if (cluster[1] - cluster[0]) <= threshold:
            output_clusters.append(cluster)
        else:
            curr_list = [i for i in range(cluster[0], cluster[1])]
            new_ranges = [(curr_list[i], curr_list[i] + threshold) if curr_list[i] + threshold <= cluster[1] else (curr_list[i], cluster[1]) for i in range(0, len(curr_list), threshold)]
            output_clusters.extend(new_ranges)
    return(output_clusters)

def substring_in_motifs(string, motifs):
    '''
    Check whether any of the motifs in a set of motifs contain a given substring.
    '''
    result = False
    string_length = len(string)
    for motif in motifs:
        if string_length <= len(motif):
            if string in motif:
                result = True
                return(result)
    return(result)

def third_position_bases(sequences, limit, reverse = False):
    '''
    Given a set of sequences, determine how frequently each base is used at third positions.
    '''
    sequences = [i for i in sequences if len(i) >= limit]
    if limit%3 != 0:
        print("The required analysis length has to be a multiple of 3!")
        raise Exception
    codon_limit = int(limit/3)
    if not reverse:
        indices = range(codon_limit)
    else:
        indices = range(-1, (-codon_limit - 1), -1)
    result = {}
    for site in range(codon_limit):
        result[site] = {i: 0 for i in _canon_bases_}
    for seq in sequences:
        third_sites_only = [seq[i] for i in range(2, len(seq), 3)]
        for pos, site in enumerate(indices):
            result[pos][third_sites_only[site]] = result[pos][third_sites_only[site]] + 1
    for site in range(codon_limit):
        current_number = sum(result[site].values())
        result[site] = {i: result[site][i]/current_number for i in result[site]}
    return(result)
                        
def trim_sequence(sequence, phase):
    '''
    Trim a sequence so that it would start and end with a full codon.
    Phase 0 means that the sequence starts with the first base of a codon.
    Phase 1: the sequence starts with the second base of a codon (there is one base missing).
    Phase 2: the sequence starts with the third base of a codon (there are two bases missing).
    '''
    #trim the start based on the phase information that was given
    if phase == 0:
        pass
    elif phase == 1:
        sequence = sequence[2:]
    elif phase == 2:
        sequence = sequence[1:]
    else:
        print("Invalid phase information!")
        sys.exit()
    #trim the end based on the length of the remaining sequence.
    if len(sequence)%3 == 0:
        return(sequence)
    elif len(sequence)%3 == 1:
        return(sequence[:-1])
    else:
        return(sequence[:-2])

def trim_sequence_report(sequence, phase):
    '''
    Rather than simply trim a sequence, also return a tuple specifying how much was trimmed.
    '''
    output = [None, None]
    #trim the start based on the phase information that was given
    if phase == 0:
        output[0] = 0
    elif phase == 1:
        output[0] = 2
        sequence = sequence[2:]
    elif phase == 2:
        output[0] = 1
        sequence = sequence[1:]
    else:
        print("Invalid phase information!")
        raise Exception
    #trim the end based on the length of the remaining sequence.
    if len(sequence)%3 == 0:
        output[1] = 0
    elif len(sequence)%3 == 1:
        output[1] = 1
        sequence = sequence[:-1]
    else:
        output[1] = 2
        sequence = sequence[:-2]
    return(output, sequence)

def unravel_consensus(consensus):
    '''
    Return all the motifs that can be generated from a consensus sequence.
    '''
    consensus = list(consensus)
    paradigms = [_IUPAC_dict_reverse_[i] for i in consensus]
    result_tuples = list(it.product(*paradigms))
    sequences = ["".join(i) for i in result_tuples]
    return(sequences)

def within_motif_skew(motifs):
    '''
    Estimate to what extent the individual motifs in the set have more skewed mononucleotide frequencies than the set as a whole.
    '''
    mononts = calc_nt_freqs(motifs, 1, return_array = True)
    freqs_array = np.zeros((len(motifs), 4))
    for pos, motif in enumerate(motifs):
        current_freqs = calc_nt_freqs([motif], 1, return_array = True)
        freqs_array[pos,:] = current_freqs
    freqs_array = abs(freqs_array - mononts)
    skew_sums = np.sum(freqs_array, axis = 1)
    return(np.mean(skew_sums))

def write_codon_usage_matrix(usage_dict, output_file):
    '''
    Write the output from codon_usage_matrix to file.
    '''
    with open(output_file, "w") as file:
        first_row = [str(i) for i in usage_dict]
        first_row = ["amino_acid", "codon"] + first_row
        file.write(",".join(first_row))
        file.write("\n")
        for aa in sorted(list(_genetic_code_.keys())):
            for codon in sorted(_genetic_code_[aa]):
                to_write = [aa, codon]
                for position in usage_dict:
                    if usage_dict[position][aa] == None:
                        to_write.append("NA")
                    else:
                        to_write.append(str(usage_dict[position][aa]["codons"][codon]))
                file.write(",".join(to_write))
                file.write("\n")

def write_minimizer_error(errors, lengths, names, error_file):
    print("Writing error to file {0}.".format(error_file))
    with open(error_file, "w") as file:
        file.write("sequence\t")
        file.write("\t".join(sorted(_canon_bases_)))
        file.write("\tlength")
        file.write("\n")
        for pos, seq in enumerate(errors):
            file.write("{0}\t".format(names[pos]))
            file.write("\t".join([str(i) for i in errors[pos]]))
            file.write("\t")
            file.write(str(lengths[pos]))
            file.write("\n")

def write_codon_usage_numbers(usage_dict, output_file):
    '''
    Write the "numbers" field from a codon usage matrix dictionary to file.
    '''
    with open(output_file, "w") as file:
        first_row = [str(i) for i in usage_dict]
        first_row = ["amino_acid"] + first_row
        file.write(",".join(first_row))
        file.write("\n")
        for aa in sorted(list(_genetic_code_.keys())):
            to_write = [aa]
            for position in usage_dict:
                if usage_dict[position][aa] == None:
                    to_write.append("NA")
                else:
                    to_write.append(str(usage_dict[position][aa]["number"]))
            file.write(",".join(to_write))
            file.write("\n")

def write_third_position_bases(usage_dict, output_file):
    '''
    Write the output dictionary from third_position_bases to file.
    '''
    with open(output_file, "w") as file:
        first_row = [str(i) for i in usage_dict]
        first_row = ["base"] + first_row
        file.write(",".join(first_row))
        file.write("\n")
        for base in _canon_bases_:
            to_write = [base]
            for position in usage_dict:
                to_write.append(str(usage_dict[position][base]))
            file.write(",".join(to_write))
            file.write("\n")
