import csv
import re
import os
import collections
from bedtools_games import bed_to_CDS_indices, CDS_to_bed_mapping, get_sequence, intersect_bed_return_bed, region_indices_to_full_indices, write_features_to_bed
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Phylo.PAML import baseml, codeml, yn00
from housekeeping import flatten, line_count, list_to_dict, make_dir, overlap, pause_script, print_elements, remove_file, run_in_parallel, run_process
import my_stats as ms
import nucleotide_comp as nc
import numpy as np
from operations_on_reads import bound_unbound_pairs, convert2bed, significant_peaks_only
import os
import random
import read_and_write as rw
import shutil
import string
import structure as struct
import tempfile
import time

#global variables
codon_dict = {'GGA': 'G', 'GGC': 'G', 'TTA': 'L', 'CTG': 'L', 'TTC': 'F', 'GAA': 'E', 'CTC': 'L', 'CGC': 'R', 'CAT': 'H', 'GTG': 'V', 'GGG': 'G',
 'GTC': 'V', 'CTA': 'L', 'GCC': 'A', 'TCG': 'S', 'CAC': 'H', 'AAG': 'K', 'AAT': 'N', 'CCG': 'P', 'TGT': 'C', 'CCC': 'P', 'CAA': 'Q',
 'CGT': 'R', 'GCG': 'A', 'ATA': 'I', 'ACA': 'T', 'TGC': 'C', 'CTT': 'L', 'ACT': 'T', 'CGA': 'R', 'TAC': 'Y', 'AGA': 'R', 'ACC': 'T',
 'GTA': 'V', 'ACG': 'T', 'TGG': 'W', 'AAC': 'N', 'CAG': 'Q', 'AGC': 'S', 'TCC': 'S', 'GCA': 'A', 'AGG': 'R', 'ATG': 'M', 'GAC': 'D',
 'TCA': 'S', 'TAT': 'Y', 'ATT': 'I', 'CCA': 'P', 'ATC': 'I', 'AAA': 'K', 'TTG': 'L', 'CCT': 'P', 'GGT': 'G', 'AGT': 'S', 'GCT': 'A',
 'TTT': 'F', 'TCT': 'S', 'GTT': 'V', 'GAG': 'E', 'CGG': 'R', 'GAT': 'D'}

cwd = os.getcwd()
user = re.search("(?<=/Users/)\w*(?=/)", cwd)
#the try...except is in case I'm running this on Watson
try:
    user = user.group(0)
except Exception:
    pass
muscle_exe = "/Users/{0}/Documents/Software/muscle".format(user)

twofold_LeuArg = [["T","T","A"],["T","T","G"],["A","G","A"],["A","G","G"]]
fourfold_LeuArg = [["C","T","T"],["C","T","C"],["C","T","A"],["C","T","G"],["C","G","T"],["C","G","C"],["C","G","A"],["C","G","G"]]

class NoDataException(Exception):
    pass

def accessibility_evolution(RBP, true_bams, mock_bam, CDS, flat_CDS, true_peaks, genome, window_size, n_sim, correspondances_file_name, alignment_folder_name, fs, transcripts, cleanup = False):
    '''
    Given one or more CLIP-seq datasets, as well as a size-matched input dataset, check for selection to sequester unbound target motifs.
    '''
    motif_length = 6
    log = {}
    #you just have to remember the last one is the mock
    bams = true_bams + [mock_bam]
    beds = []

    for pos, bam in enumerate(bams):
        bed = bam[:-4] + "_reads.bed"
        temp_bam = bam[:-4] + "_temp.bam"
        log[bam] = {}
        beds.append(bed)
        #only leave the 2nd read
        run_process(["samtools", "view", "-hb", "-f", "128", bam], file_for_output = temp_bam)
        convert2bed(temp_bam, bed, group_flags = False)
        remove_file(temp_bam)
        if cleanup:
            remove_file(bam)
        record_count = line_count(bed)
        print("R{0}, raw reads: {1}.---{2}".format(pos + 1, record_count, RBP))
        log[bam]["raw count"] = record_count
        #only leave the reads that overlap with CDSs
        intersect_bed_return_bed(bed, flat_CDS, overlap = False, write_both = False, sort = False,
                                 force_strand = True, modify_chr_ids = True, no_dups = True, use_bedops = False, output_file = bed)
        record_count = line_count(bed)
        print("R{0}, overlapping reads: {1}.---{2}\n".format(pos + 1, record_count, RBP))
        log[bam]["overlap count"] = record_count

    all_motifs = []
    for pos, peak_file in enumerate(true_peaks):
        log[peak_file] = {}
        record_count = line_count(peak_file)
        print("R{0}, raw peaks: {1}.---{2}".format(pos + 1, record_count, RBP))
        log[peak_file]["raw count"] = record_count
        new_peak_file = peak_file[:-4] + "_clean.bed"
        #only leave peaks that overlap with CDSs
        intersect_bed_return_bed(peak_file, flat_CDS, overlap = False, write_both = False, sort = False,
                                 force_strand = True, modify_chr_ids = True, no_dups = True, use_bedops = False, output_file = new_peak_file)
        record_count = line_count(new_peak_file)
        print("R{0}, overlapping peaks: {1}.---{2}".format(pos + 1, record_count, RBP))
        log[peak_file]["overlap count"] = record_count
        signif_peak_file = peak_file[:-4] + "_signif.bed"
        #only leave peaks that overlap significantly more reads than the corresponding region in the mock
        significant_peaks_only(new_peak_file, beds[pos], beds[-1], log[bams[pos]]["raw count"], log[bams[-1]]["raw count"], signif_peak_file, threshold = 0.05)
        record_count = line_count(signif_peak_file)
        print("R{0}, significant peaks: {1}.---{2}\n".format(pos + 1, record_count, RBP))
        log[peak_file]["significant count"] = record_count
        #otherwise fasta_from_intervals for RSAT won't work
        modified_signif_peak_file = peak_file[:-4] + "_signif_modified.bed"
        peaks = rw.read_many_fields(signif_peak_file, "\t")
        for position in range(len(peaks)):
            peaks[position][0] = peaks[position][0].lstrip("chr")
        rw.write_to_csv(peaks, modified_signif_peak_file, "\t")
        #use RSAT to see which motifs are over-represented witin the peaks
        motifs = nc.motifs_from_peaks(modified_signif_peak_file, genome, motif_length, 20, output_file = signif_peak_file[:-4] + "_motifs.txt")
        print(motifs[:10])
        print("({0}, R{1})".format(RBP, pos + 1))
        print("\n")
        all_motifs.append(motifs[:10])
        true_peaks[pos] = signif_peak_file
        if cleanup:
            remove_file(new_peak_file)
            remove_file(modified_signif_peak_file)

    for pos, bam in enumerate(true_bams):

        bed = bam[:-4] + "_reads.bed"

        pairs, pairs_coords = bound_unbound_pairs(CDS, bed, true_peaks[pos], genome, motifs[:10], motifs, window_size, bed[:-4] + "_peak_hits.fasta", bed[:-4] + "_desert_hits.fasta")
        pairs_number = len(pairs)
        print("R{0}, peak-desert pairs: {1}.---{2}".format(pos + 1, pairs_number, RBP))
        log[bam]["pairs number"] = pairs_number

        acc_results = effect_on_structure(pairs, pairs_coords, window_size, CDS, motifs[:10], motifs, n_sim, correspondances_file_name, alignment_folder_name, transcripts, fs, genome)
        median_acc_ratio = acc_results["median ratio"]
        print("R{0}, median acc. ratio: {1}.---{2}".format(pos + 1, median_acc_ratio, RBP))
        log[bam]["median acc. ratio"] = median_acc_ratio
        acc_p = acc_results["p"]
        print("R{0}, acc. p: {1}.".format(pos + 1, acc_p))
        log[bam]["acc. p"] = acc_p

        if "mutation score" in acc_results:
            site_count = acc_results["site count"]
            mutation_score = acc_results["mutation score"]
            print("R{0}, site count: {1}.---{2}".format(pos + 1, site_count, RBP))
            print("R{0}, mutation score: {1}---{2}.".format(pos + 1, mutation_score, RBP))
            log[bam]["mutation score"] = mutation_score
            log[bam]["site count"] = site_count

            mutation_p = acc_results["mutation p"]
            print("R{0}, mutation p: {1}.---{2}\n".format(pos + 1, mutation_p, RBP))
            log[bam]["mutation p"] = mutation_p
        else:
            log[bam]["mutation score"] = None
            log[bam]["mutation p"] = None

        if cleanup:
            remove_file(bed)
    if cleanup:
        remove_file(beds[-1])

    motif_overlap = len(overlap(all_motifs[0], all_motifs[1]))/len(motifs)
    print("Motif overlap: {0}.---{1}\n\n".format(motif_overlap, RBP))
    log["motif overlap"] = motif_overlap
           
    return(log)

def blast_all_against_all(db_name, fasta_file_name, output_file_name):
    '''
    Blast all the sequences in a fasta file against each-other.
    '''
    run_process(["makeblastdb", "-in", fasta_file_name, "-out",
                 "/Users/{0}/Documents/Software/ncbi-blast-2.2.30+/db/{1}".format(user, db_name),
                 "-dbtype", "nucl"])
    run_process(["Blastn", "-task", "blastn", "-query", fasta_file_name,
                 "-db", "/Users/{0}/Documents/Software/ncbi-blast-2.2.30+/db/{1}".format(user, db_name),
                 "-out", output_file_name, "-outfmt", "10", "-evalue", "1e-04", "-num_threads", str(int((os.cpu_count()/2)-1))])

def calculate_dS_core(motifs, motif_lengths, input_dict_file_name, alignment_folder_name, check_density = False, map_from_regions = False, region_name = None, method = None, positions_file = None, exclude_positions = None, sequence_file = None, remove_CG = False, only_4f = False):
    '''
    Core function for dS_from_hits below.
    '''
    #if no method is specified for calculating dS, use the Yang and Nielsen (2000) model
    if not method:
        method = "yn"
    temp_file_names = ["temp_data/temp_alignment{0}.phy".format(random.random()) for i in range(4)]
    two_seqs = False
    if method == "baseml":
        two_seqs = True
    if positions_file:
        pos_file = open(positions_file, "w")
    if exclude_positions:
        exclude_file = open(exclude_positions)
    if sequence_file:
        make_dir("temp_data/seq_files_for_ds")
        if sequence_file == "generate":
            sequence_file = "temp_data/seq_files_for_ds/sim_seq_{0}.txt".format(random.random())                    
        seq_file = open(sequence_file, "w")
    with open(temp_file_names[0], "w") as output_file, open(temp_file_names[1], "w") as orth_output_file, open(temp_file_names[2], "w") as header_file, open(input_dict_file_name) as input_file:
        #write a phylip file containing a concatenation of all the bases overlapping motif hits in a set of sequences, aligned to the orthologous sequence in another species
        #most of the input data is read from the file in input_dict_file_name (prepared in input_dict_for_dS below)
        output_file.write("id  ")
        orth_output_file.write("orth_id  ")
        aligning_positions = []
        #each line in the input_dict_file corresponds to one sequence
        for line in input_file:
            #parse current line into a dictionary
            idn, input_dict = parse_input_dict_line(line)
            #if you're not analyzing full CDSs but rather, say, exonic subregions
            if map_from_regions:
                full_CDS = input_dict["full CDS"][region_name]
                #get the positions of motif hits within the sequence
                positions = nc.get_motif_set_density(motifs, motif_lengths, full_CDS, concat = True)["positions"]
                #convert them to positions relative to the corresponding full CDS
                positions = region_indices_to_full_indices(positions, input_dict["flank indices in CDS"][region_name])
                #get the positions of two-fold and four-fold Leucine/Arginine codons
                two_fold = input_dict["two-fold"][region_name]
                four_fold = input_dict["four-fold"][region_name]
            else:
                full_CDS = input_dict["full CDS"]
                positions = nc.get_motif_set_density(motifs, motif_lengths, full_CDS, concat = True)["positions"]
                two_fold = input_dict["two-fold"]
                four_fold = input_dict["four-fold"]
            #expand motif hit positions to get full codons
            if only_4f:
                positions = [i for i in positions if i in input_dict["all_fourfold"]]
            if not two_seqs:
                [positions, to_change_two_fold, to_change_four_fold] = fill_codons(positions, two_fold, four_fold)
            if remove_CG:
                positions = [i for i in positions if i not in input_dict["CG"]]
            if positions_file:
                pos_file.write("{0}\t{1}\n".format(idn, ",".join([str(i) for i in positions])))
            if exclude_positions:
                next_line = exclude_file.readline()
                next_line = next_line.rstrip("\n")
                next_line = next_line.split("\t")
                if next_line[0] != idn:
                    print("IDs don't match!")
                    print(idn)
                    print(next_line)
                    raise Exception
                to_exclude = next_line[1].split(",")
                to_exclude = [int(i) for i in to_exclude if i]
                positions = [i for i in positions if i not in to_exclude]
            if sequence_file:
                current_sequence = "".join([full_CDS[i] for i in positions])
                seq_file.write(current_sequence)
            if positions:
                #the aligned full CDS from either species
                aligned_sequences = input_dict["aligned sequences"]
                #convert positions in the CDS to positions in the aligned sequence (which has dashes for gaps)
                aligned_positions = get_aligned_positions(aligned_sequences, positions)                     
                if not two_seqs and (len(aligned_positions)%3) != 0:
                    print("Problem filling up the codons!")
                    raise Exception
                if not two_seqs:
                    #convert the positions of two- and fourfold degenerate Arginine/Leucine codons to indices in the aligned sequence
                    to_change_two_fold = get_aligned_positions(aligned_sequences, to_change_two_fold)
                    to_change_four_fold = get_aligned_positions(aligned_sequences, to_change_four_fold)
                    #change the Leucine/Arginine codons to either Glutamic acid or Alanine codons depending on the degeneracy of the first base
                    temp_aligned_sequences = changeLeuArg(to_change_two_fold, to_change_four_fold, aligned_sequences)
                else:
                    temp_aligned_sequences = aligned_sequences.copy()
                #extract the bases that correspond to the motif hit positions from the aligned sequence, for both species
                temp_seq = [0 for pos in aligned_positions]
                temp_orth_seq = temp_seq.copy()
                for enum, pos in enumerate(aligned_positions):
                    temp_seq[enum] = temp_aligned_sequences[0][pos]
                    temp_orth_seq[enum] = temp_aligned_sequences[1][pos]
                if not two_seqs:
                    #remove in-frame stop codons from the sequences or codeml will sit quietly and do nothing
                    [temp_seq, temp_orth_seq] = remove_stops(temp_seq, temp_orth_seq)
                #write motif hit sequences to file
                output_file.write("".join(temp_seq))
                orth_output_file.write("".join(temp_orth_seq))
        output_file.write("\n")
        output_file.close()
        seq_length = os.stat(temp_file_names[0]).st_size
        seq_length = seq_length - 5
        if positions_file:
            pos_file.close()
        if exclude_positions:
            exclude_file.close()
        if sequence_file:
            seq_file.close()
        if seq_length > 0:
            #make the first line of the phylip file
            header_file.write(" 2 {0}\n".format(seq_length))
            for i in [orth_output_file, header_file]:
                i.close()
            #concatenate the header line, the sequence from the first species (usually human) and the sequence from the second species (usually macaque)
            run_process(["cat", temp_file_names[2], temp_file_names[0], temp_file_names[1]], file_for_output = temp_file_names[3])
            #use PAML codeml to calculate dS
            cml_input = temp_file_names[3]
            cml_output = "{0}.out".format(temp_file_names[3][:-4])
            ctl_file = "{0}.ctl".format(random.random())
            if two_seqs:
                variable = "tree length"
            else:
                variable = "dS"
            ds = run_codeml(cml_input, cml_output, ctl_file = ctl_file, method = method)[variable]
            #clean up temp files
            for i in temp_file_names:
                os.remove(i)
            os.remove(cml_output)
            os.remove(ctl_file)
        else:
            return(None)
    return(ds)

def calculate_dS_core_no_degeneracy(motifs, motif_lengths, input_dict_file_name, alignment_folder_name, check_density = False):
    '''
    Core function for dS_from_hits below (without allowing for sites where the two species differ yet both are motif).
    '''
    temp_file_names = ["temp_data/temp_alignment{0}.phy".format(random.random()) for i in range(4)]
    with open(temp_file_names[0], "w") as output_file, open(temp_file_names[1], "w") as orth_output_file, open(temp_file_names[2], "w") as header_file, open(input_dict_file_name) as input_file:
        #write a phylip file containing a concatenation of all the bases overlapping motif hits in a set of sequences, aligned to the orthologous sequence in another species
        #most of the input data is read from the file in input_dict_file_name (prepared in input_dict_for_dS below)
        output_file.write("id  ")
        orth_output_file.write("orth_id  ")
        aligning_positions = []
        #each line in the input_dict_file corresponds to one sequence
        for line in input_file:
            #parse current line into a dictionary
            idn, input_dict = parse_input_dict_line(line)
            aligned_sequences = input_dict["aligned sequences"]
            full_CDS = input_dict["full CDS"]
            full_other_CDS = "".join([i for i in aligned_sequences[1] if i != "-"])
            positions = nc.get_motif_set_density(motifs, motif_lengths, full_CDS, concat = True)["positions"]
            other_positions = nc.get_motif_set_density(motifs, motif_lengths, full_other_CDS, concat = True)["positions"]
            two_fold = input_dict["two-fold"]
            four_fold = input_dict["four-fold"]
            #expand motif hit positions to get full codons
            [positions, to_change_two_fold, to_change_four_fold] = fill_codons(positions, two_fold, four_fold)
            #convert positions in the CDS to positions in the aligned sequence (which has dashes for gaps)
            aligned_positions = get_aligned_positions(aligned_sequences, positions)
            other_aligned_positions = get_aligned_positions(list(reversed(aligned_sequences)), other_positions)
            if len(aligned_positions)%3 != 0:
                print("Problem filling up the codons!")
                raise Exception
            #convert the positions of two- and fourfold degenerate Arginine/Leucine codons to indices in the aligned sequence
            to_change_two_fold = get_aligned_positions(aligned_sequences, to_change_two_fold)
            to_change_four_fold = get_aligned_positions(aligned_sequences, to_change_four_fold)
            aligned_positions = remove_degenerate_differences(aligned_sequences, aligned_positions, other_aligned_positions)
            #change the Leucine/Arginine codons to either Glutamic acid or Alanine codons depending on the degeneracy of the first base
            temp_aligned_sequences = changeLeuArg(to_change_two_fold, to_change_four_fold, aligned_sequences)
            #extract the bases that correspond to the motif hit positions from the aligned sequence, for both species
            temp_seq = [0 for pos in aligned_positions]
            temp_orth_seq = temp_seq.copy()
            for enum, pos in enumerate(aligned_positions):
                temp_seq[enum] = temp_aligned_sequences[0][pos]
                temp_orth_seq[enum] = temp_aligned_sequences[1][pos]
            #remove in-frame stop codons from the sequences or codeml will sit quietly and do nothing
            [temp_seq, temp_orth_seq] = remove_stops(temp_seq, temp_orth_seq)
            #write motif hit sequences to file
            output_file.write("".join(temp_seq))
            orth_output_file.write("".join(temp_orth_seq))
        output_file.write("\n")
        output_file.close()
        seq_length = os.stat(temp_file_names[0]).st_size
        seq_length = seq_length - 5
        if seq_length > 0:
            #make the first line of the phylip file
            header_file.write(" 2 {0}\n".format(seq_length))
            for i in [orth_output_file, header_file]:
                i.close()
            #concatenate the header line, the sequence from the first species (usually human) and the sequence from the second species (usually macaque)
            run_process(["cat", temp_file_names[2], temp_file_names[0], temp_file_names[1]], file_for_output = temp_file_names[3])
            #use PAML codeml to calculate dS
            cml_input = temp_file_names[3]
            cml_output = "{0}.out".format(temp_file_names[3][:-4])
            ctl_file = "{0}.ctl".format(random.random())
            variable = "dS"
            ds = run_codeml(cml_input, cml_output, ctl_file = ctl_file, method = "gy")[variable]
            #clean up temp files
            for i in temp_file_names:
                os.remove(i)
            os.remove(cml_output)
            os.remove(ctl_file)
        else:
            return(None)
    return(ds)

def changeLeuArg(to_change_two_fold, to_change_four_fold, aligned_sequences):
    '''
    given the positions of Leucine and Arginine codons in the first sequence of a pair of aligned sequences and the sequences themselves,
    it changes those that are 2-fold degenerate into Glutamic acid codons and those that are 4-fold degenerate
    into Alanine codons (w/o changing the third base of the codon)
    this is done only for Leucine and Arginine codons that were introduced because the third (and maybe second) base of the codon are part of a motif hit
    but where the first base is not part of the motif hit. the goal is to prevent such first sites from being counted as synonymous by PAML
    all while preserving any divergences from macaque at the third site.
    '''
    #Glutamic acid
    two_fold_replacement = ["G","A"]
    #Alanine
    four_fold_replacement = ["G","C"]
    for i in to_change_two_fold:
        for j in range(len(aligned_sequences)):
            #always change the human sequence. don't change the macaque sequence if the codon is an indel (or you'll end up with a codon like "GA-")
            if not (aligned_sequences[j][i] == "-" and j == 1):
                aligned_sequences[j][i:i+2] = two_fold_replacement
    for i in to_change_four_fold:
        for j in range(len(aligned_sequences)):
            if not (aligned_sequences[j][i] == "-" and j == 1):
                aligned_sequences[j][i:i+2] = four_fold_replacement
    return(aligned_sequences)

def check_disruption(pos, var_base, human_seq, motifs, motif_lengths):
    '''
    Given a position in a sequence and the base observed at that position in another species, check whether that substitution decreases the motif content of the local region.
    '''
    max_length = max(motif_lengths)
    start_pos = pos - max_length + 1
    start_trim = 0
    if start_pos < 0:
        start_pos = 0
        start_trim = pos - max_length + 1
    end_pos = pos + max_length
    if end_pos > len(human_seq):
        end_pos = len(human_seq)
    fragment = human_seq[start_pos:end_pos]
    current_dens = nc.get_motif_set_density(motifs, motif_lengths, fragment)["density"]
    fragment = list(fragment)
    fragment[max_length - 1 + start_trim] = var_base
    new_dens = nc.get_motif_set_density(motifs, motif_lengths, "".join(fragment))["density"]
    if new_dens < current_dens:
        return(True)
    else:
        return(False)

def check_ORF_integrity(sequence, PTC_check = True):
    '''
    Given a sequence, check whether it's a nice and clean ORF.
    '''
    #check that the lentgh is a multiple of three
    if len(sequence)%3 != 0:
        return(False, "incorrect length")
    #check start codon
    if sequence[:3] != "ATG":
        return(False, "incorrect start")
    #check stop
    if sequence[-3:] not in ["TAA", "TGA", "TAG"]:
        return(False, "incorrect stop")
    #check that there are no non-canonical bases
    for i in sequence:
        if i not in ["A","T","C","G"]:
            return(False, "invalid base")
    #check that there are no premature termination codons
    if PTC_check:
        seq_length = len(sequence)
        for i in range(0, seq_length, 3):
            current_codon = sequence[i:i+3]
            if i+3 != seq_length and current_codon in ["TAG", "TAA", "TGA"]:
                return(False, "PTC")
    return(True, None)

def cons_by_dinucl(CDS_fasta_file, motifs, correspondances_file, alignment_folder_name, dinucl, picked = None, map_from_regions = None):
    '''
    For each dinucleotide, calculate its conservation within motifs vs elsewhere.
    '''
    #necessary for getting motif positions
    motif_lengths = [len(i) for i in motifs]
    #convert motifs from strings to lookahead regexes (lookahead so you would catch overlaps)
    motifs = nc.motif_to_regex(motifs)
    #prepare output dictionary
    result = {i: {"subst. in motifs": 0, "frequency in motifs": 0, "subst. in non-motifs": 0, "frequency in non-motifs": 0} for i in dinucl}
    #read in input fasta
    names, seqs = rw.read_fasta(CDS_fasta_file)
    #file that has correspondances between tehe gene identifiers of the focal species and of the ortholog
    correspondances = rw.read_many_fields(correspondances_file, ",")
    correspondance_dict = {}
    for i in correspondances:
        correspondance_dict[i[0]] = i[1]
    #if picked is specified (a list of gene identifiers), only analyze those genes, otherwise analyze everything
    if not picked:
        picked = names[:]
    motif_dinucl_sum = 0
    non_motif_dinucl_sum = 0
    #loop over the sequences
    for pos, seq in enumerate(seqs):
        if pos % 1000 == 0:
            print(pos)
        current_name = names[pos]
        if current_name in picked:
            #if analyzing CDS subregions rather than full CDSs
            if map_from_regions:
                #switch from the name used in the fasta to the full CDS identifier
                current_gene_name = map_from_regions[current_name]["idn"]
            else:
                current_gene_name = current_name
            #get the sequence positions that overlap with the motifs
            motif_pos = nc.get_motif_set_density(motifs, motif_lengths, seq, concat = True)["positions"]
            #get positions of fourfold degenerate sites
            fourfold_deg = nc.get_4fold_deg(seq)
            #only consider motif hit positions that are fourfold degenerate
            motif_pos = [i for i in motif_pos if i in fourfold_deg]
            if map_from_regions:
                #map motif hit positions from relative indices in the subregion to indices in the full CDS
                region_pos = fourfold_deg.copy()
                fourfold_deg = region_indices_to_full_indices(fourfold_deg, map_from_regions[current_name]["flank indices in CDS"])
            #get the identifier of the orthologous gene
            orth_idn = correspondance_dict[current_gene_name]
            #stupid historical nonsense
            if "_" not in orth_idn:
                orth_idn = orth_idn + "_0"
            #go to the folder that contains phylip files with alignments between the CDSs of the two species (created in keep_conserved_pc), parse the relevant one
            phy_file_name = "{0}/{1}_{2}.phy".format(alignment_folder_name, current_gene_name, orth_idn)
            aligned_sequences = [list(str(i.seq)) for i in SeqIO.parse(phy_file_name, "phylip-sequential")]
            #convert the fourfold degenerate positions obtained in the unaligned sequence to indices in the aligned sequence (i.e. account for
            #shifts due to indels)
            aligned_fourfold = get_aligned_positions(aligned_sequences, fourfold_deg)
            #loop over the fourfold degenerate positions
            for pos2, site in enumerate(aligned_fourfold):
                in_motif = False
                change = False
                #get position in the unaligned sequence (you want to be extracting the dinucleotides
                #and you can't use the aligned sequence for that because the position might be next to a gap)
                if map_from_regions:
                    raw_pos = region_pos[pos2]
                else:
                    raw_pos = fourfold_deg[pos2]
                #get the two dinucleotides that overlap the position
                current_dinucl = [seq[(raw_pos - 1): (raw_pos + 1)], seq[raw_pos: (raw_pos + 2)]]
                #this is in case the first or the last third position is 4-fold degenerate. This obviously won't happen with full ORFs
                #but might with subregions of the CDS
                current_dinucl = [i for i in current_dinucl if len(i) > 1]
                #if this fourfold degenerate site overlaps with a motif hit
                if raw_pos in motif_pos:
                    in_motif = True
                #if macaque or whoever has the same base as human or whoever
                if aligned_sequences[0][site] != aligned_sequences[1][site]:
                    change = True
                #update the counters based on whether or not the dinucleotide is part of a motif
                #and whether or not it's different in the two species
                for dint in current_dinucl:
                    if in_motif:
                        motif_dinucl_sum = motif_dinucl_sum + 1
                        result[dint]["frequency in motifs"] = result[dint]["frequency in motifs"] + 1
                        if change:
                            result[dint]["subst. in motifs"] = result[dint]["subst. in motifs"] + 1
                    else:
                        non_motif_dinucl_sum = non_motif_dinucl_sum + 1
                        result[dint]["frequency in non-motifs"] = result[dint]["frequency in non-motifs"] + 1
                        if change:
                            result[dint]["subst. in non-motifs"] = result[dint]["subst. in non-motifs"] + 1
    #convert counts into frequencies
    for dint in dinucl:
        if result[dint]["frequency in motifs"] == 0:
            result[dint]["subst. in motifs"] = None
            result[dint]["subst. in non-motifs"] = None
            result[dint]["frequency in non-motifs"] = None
        else:
            result[dint]["subst. in motifs"] = result[dint]["subst. in motifs"]/result[dint]["frequency in motifs"]
            result[dint]["frequency in motifs"] = result[dint]["frequency in motifs"]/motif_dinucl_sum
            if result[dint]["frequency in non-motifs"] == 0:
                result[dint]["subst. in non-motifs"] = None
            else:
                result[dint]["subst. in non-motifs"] = result[dint]["subst. in non-motifs"]/result[dint]["frequency in non-motifs"]
                result[dint]["frequency in non-motifs"] = result[dint]["frequency in non-motifs"]/non_motif_dinucl_sum
    return(result)

def cons_by_dinucl_no_control(positions, fasta_file, correspondances_file, alignment_folder_name, dinucl):
    '''
    For each dinucleotide, calculate its conservation within a set of sites.
    '''
    #prepare output dictionary
    result = {i: {"rate": 0, "frequency": 0} for i in dinucl}
    #read in input fasta
    names, seqs = rw.read_fasta(fasta_file)
    #file that has correspondances between the gene identifiers of the focal species and of the ortholog
    correspondances = rw.read_many_fields(correspondances_file, ",")
    correspondance_dict = {}
    for i in correspondances:
        correspondance_dict[i[0]] = i[1]
    dinucl_sum = 0
    #loop over the sequences
    for pos, seq in enumerate(seqs):
        if pos % 1000 == 0:
            print(pos)
        current_name = names[pos]
        if current_name in positions:
            current_pos = positions[current_name]
            #get the identifier of the orthologous gene
            orth_idn = correspondance_dict[current_name]
            #go to the folder that contains phylip files with alignments between the CDSs of the two species (created in keep_conserved_pc), parse the relevant one
            phy_file_name = "{0}/{1}_{2}.phy".format(alignment_folder_name, current_name, orth_idn)
            aligned_sequences = [list(str(i.seq)) for i in SeqIO.parse(phy_file_name, "phylip-sequential")]
            #convert the fourfold degenerate positions obtained in the unaligned sequence to indices in the aligned sequence (i.e. account for
            #shifts due to indels)
            aligned_pos = get_aligned_positions(aligned_sequences, current_pos)
            #loop over the fourfold degenerate positions
            for pos2, site in enumerate(aligned_pos):
                change = False
                #get position in the unaligned sequence (you want to be extracting the dinucleotides
                #and you can't use the aligned sequence for that because the position might be next to a gap)
                raw_pos = current_pos[pos2]
                #get the two dinucleotides that overlap the position
                current_dinucl = [seq[(raw_pos - 1): (raw_pos + 1)], seq[raw_pos: (raw_pos + 2)]]
                #if macaque or whoever has a different base to human or whoever
                if aligned_sequences[0][site] != aligned_sequences[1][site]:
                    change = True
                #update the counters based on whether or not the dinucleotide is different in the two species
                for dint in current_dinucl:
                    dinucl_sum = dinucl_sum + 1
                    result[dint]["frequency"] = result[dint]["frequency"] + 1
                    if change:
                        result[dint]["rate"] = result[dint]["rate"] + 1
    #convert counts into frequencies
    for dint in dinucl:
        if result[dint]["frequency"] == 0:
            result[dint]["rate"] = None
        else:
            result[dint]["rate"] = result[dint]["rate"]/result[dint]["frequency"]
            result[dint]["frequency"] = result[dint]["frequency"]/dinucl_sum
    return(result)


def dS_from_hits(motifs, alignment_folder_name, input_dict_file_name, n_sim = None, simulants = None, sim_output_file_name = None, map_from_regions = False, region_name = None, method = None, no_degeneracy = False, no_overlaps = None, remove_CG = False, only_4f = False):
    '''
    Given a set of motifs, sets of simulant motifs and a set of sequences, calculate dS within the motifs and the simulants.
    '''
    #if no method is specified for calculating dS, use the Yang and Nielsen (2000) model
    if not method:
        method = "yn"
    motif_lengths = [len(i) for i in motifs]
    motifs = nc.motif_to_regex(motifs)
    #most of the input data is read from the file in input_dict_file_name (prepared in input_dict_for_dS below)
    #get the dS of sites overlapping hits to the true motifs
    if no_degeneracy:
        real_ds = calculate_dS_core_no_degeneracy(motifs, motif_lengths, input_dict_file_name, alignment_folder_name, check_density = True)
    else:
##        real_ds = calculate_dS_core(motifs, motif_lengths, input_dict_file_name, alignment_folder_name, check_density = True, map_from_regions = map_from_regions, region_name = region_name, method = method, positions_file = no_overlaps, remove_CG = remove_CG)
        real_ds = calculate_dS_core(motifs, motif_lengths, input_dict_file_name, alignment_folder_name, check_density = True, map_from_regions = map_from_regions, region_name = region_name, method = method, positions_file = no_overlaps, sequence_file = "temp_data/seq_files_for_ds/true_seq.txt", remove_CG = remove_CG, only_4f = only_4f)
    if real_ds == None or not n_sim:
        return(real_ds)
    else:
        simulants = [nc.motif_to_regex(i) for i in simulants]
        #in a parallel manner, calculate the dS at sites overlapping motifs from each of the simulant sets
        if no_degeneracy:
            kwargs_dict = dict(check_density = False)
            sim_ds = run_in_parallel(simulants, ["foo", motif_lengths, input_dict_file_name, alignment_folder_name], calculate_dS_core_no_degeneracy, kwargs_dict = kwargs_dict, onebyone = True)
        else:
##            kwargs_dict = dict(check_density = False, map_from_regions = map_from_regions, region_name = region_name, method = method, exclude_positions = no_overlaps, remove_CG = remove_CG)
            kwargs_dict = dict(check_density = False, map_from_regions = map_from_regions, region_name = region_name, method = method, exclude_positions = no_overlaps, sequence_file = "generate", remove_CG = remove_CG, only_4f = only_4f)
            sim_ds = run_in_parallel(simulants, ["foo", motif_lengths, input_dict_file_name, alignment_folder_name], calculate_dS_core, kwargs_dict = kwargs_dict, onebyone = True)
        #the reason I am using this quaint way of calculating the mean is that
        #I don't want to have to make a massive array with all the sim_ds
        #same goes for calculating the p-value
        sim_sum = 0
        sim_counter = 0
        smaller_than = 0
        with open(sim_output_file_name, "w") as output_file:
            for i in sim_ds:
                #fetch the data from the parallel subprocesses
                current_sim_ds = i.get()
                #write the simulated dS to file
                output_file.write("{0}\n".format(current_sim_ds))
                if current_sim_ds != None:
                    sim_counter = sim_counter + 1
                    sim_sum = sim_sum + current_sim_ds
                    if current_sim_ds <= real_ds:
                        smaller_than = smaller_than + 1
        if sim_sum > 0:
            mean_sim_ds = sim_sum/sim_counter
            norm_ds = (real_ds - mean_sim_ds)/mean_sim_ds
            p = (smaller_than + 1)/(sim_counter + 1)
        else:
            mean_sim_ds = None
            norm_ds = None
            p = None
        return({"dS": real_ds, "mean simulated dS": mean_sim_ds, "normalized dS": norm_ds, "effective p": p})

def effect_on_structure(pairs, pairs_coords, window_size, CDS, motifs, motifs_extended, n_sim, correspondances_file_name, alignment_folder_name, transcripts, fs, genome):
    #number of fourfold degenerate sites within window_size bp to a desert hit that differ between the two species
    site_count = 0
    #number of such sites where the orthologous base increases site accessibility
    increase_count = 0
    #same for simulants
    sim_site_counts = [0 for i in range(n_sim)]
    sim_increase_counts = [0 for i in range(n_sim)]
    motif_length = len(motifs[0])
    #position of the first motif base within the window_size long region
    first_base = int((window_size - motif_length)/2)
    last_base = first_base + motif_length - 1
    gene_name_dict = fs.get_gene_name_dict(transcripts)
    correspondances = rw.read_many_fields(correspondances_file_name, ",")
    correspondance_dict = list_to_dict(correspondances, 0, 1)
    plfold_output_folder = "temp_data/temp_folder{0}".format(random.random())
    make_dir(plfold_output_folder)
    output = {}
    real_acc = struct.accessibility_paired(pairs, window_size, motifs, output_folder = plfold_output_folder)
    output["median ratio"] = real_acc["median ratio"]
    output["p"] = real_acc["p"]
    category = 1
    for pos, pair in enumerate(pairs):
        print(pos)
        current_record = pairs_coords[pos][category]
        current_trans_id, CDS_pos, record_mapping = effect_on_structure_prepare(current_record, CDS, fs, gene_name_dict, window_size)
        current_orth_id = correspondance_dict[current_trans_id]
        phy_file_name = "{0}/{1}_{2}.phy".format(alignment_folder_name, current_trans_id, current_orth_id)
        aligned_sequences = [list(i.seq) for i in SeqIO.parse(phy_file_name, "phylip-sequential")]
        current_CDS_seq = "".join([i for i in aligned_sequences[0] if i != "-"])
        fourfold_deg = nc.get_4fold_deg(current_CDS_seq)
        CDS_pos = [i for i in CDS_pos if i in fourfold_deg]
        if CDS_pos:
            possible_pos = {j: [base_pos for base_pos, i in enumerate(aligned_sequences[0]) if (i == j) and (base_pos in fourfold_deg)] for j in nc._canon_bases_}
            true_score = real_acc["score list"][category][pos]
            CDS_pos = get_aligned_positions(aligned_sequences, CDS_pos)
            for position in CDS_pos:
                if aligned_sequences[0][position] != aligned_sequences[1][position]:
                    orth_base = aligned_sequences[1][position]
                    score, orig_pos = effect_on_structure_mutated_score(pair, category, position, record_mapping, window_size, motif_length, last_base, orth_base, output_folder = plfold_output_folder)
                    site_count = site_count + 1
                    if score > true_score:
                        increase_count = increase_count + 1
                    true_base = aligned_sequences[0][position]
                    true_coords = pairs_coords[pos][category]
                    for sim in range(n_sim):
                        sim_position = random.choice(possible_pos[true_base])
                        sim_real_score, sim_sim_score = effect_on_structure_simulation(last_base, sim_position, position, true_coords, genome, orth_base, orig_pos, window_size, motif_length, output_folder = plfold_output_folder)
                        sim_site_counts[sim] = sim_site_counts[sim] + 1
                        if sim_sim_score > sim_real_score:
                            sim_increase_counts[sim] = sim_increase_counts[sim] + 1
                            
    if site_count > 0:
        mutation_score = increase_count/site_count
        sim_mutation_scores = np.divide(np.array(sim_increase_counts), np.array(sim_site_counts))
        mutation_p = ms.calc_eff_p(mutation_score, sim_mutation_scores, greater = False)
    else:
        mutation_score = None
        mutation_p = None
    output["mutation score"] = mutation_score
    output["site count"] = site_count
    output["mutation p"] = mutation_p
    shutil.rmtree(plfold_output_folder)
    return(output)

def effect_on_structure_mutated_score(pair, category, position, record_mapping, window_size, motif_length, last_base, orth_base, output_folder = None):
    temp_seq = list(pair[category])
    orig_pos = record_mapping[position]
    temp_seq[orig_pos] = orth_base
    temp_seq = "".join(temp_seq)
    structure_dict = struct.get_bp_prob(temp_seq, int(window_size - motif_length), int(window_size - motif_length), motif_length, output_folder = output_folder)
    score = structure_dict[last_base, (motif_length - 1)]
    return(score, orig_pos)


def effect_on_structure_prepare(current_record, CDS, fs, gene_name_dict, window_size, test_mode = False):
    current_record[0] = current_record[0].lstrip("chr")
    current_CDS = CDS[current_record[3]]
    temp = "temp_data/temp_bed_file{0}.bed".format(random.random())
    rw.write_to_csv([current_record], temp, "\t")
    #the regions can also overlap with introns so you need to be sure to only get the CDS bits for the next step to work
    #remove phase information
    CDS_for_bedops = [i[0] for i in current_CDS]
    neg_strand = False
    if current_record[5] == "-":
        neg_strand = True
        CDS_for_bedops = [i for i in reversed(CDS_for_bedops)]
    record_mapping = CDS_to_bed_mapping(current_record, CDS_for_bedops)
    CDS_only_records = intersect_bed_return_bed(temp, CDS_for_bedops, use_bedops = True, chrom = current_record[0], intersect = True)
    if not CDS_only_records:
        print("No parts of the peak overlap with CDSs. This should be impossible.")
        print(CDS_for_bedops)
        print(current_record)
        raise Exception

    if len(CDS_only_records) > 1:
        print("The peak is spread out across more than one exon.")
        if not test_mode:
            print(current_record)
            print(current_CDS)
            print(CDS_only_records)
            raise Exception
    CDS_pos = []
    #there could be several if there is a really tiny intron within the region
    #but then orig_pos would be wrong so I ned to make sure there are no such cases
    for record in CDS_only_records:
        CDS_pos.extend(bed_to_CDS_indices(record, current_CDS))
    CDS_pos = sorted(CDS_pos)
    os.remove(temp)
    if len(CDS_pos) > window_size:
        print("Too many positions were retrieved in the CDS. Something's not right.")
        print(CDS_for_bedops)
        print(current_record)
        print(CDS_pos)
        raise Exception
    return(CDS_for_bedops[0][4], CDS_pos, record_mapping)

def effect_on_structure_simulation(last_base, sim_position, position, true_coords, genome, orth_base, orig_pos, window_size, motif_length, output_folder = None):
    move = sim_position - position
    sim_coords = true_coords.copy()
    sim_coords[1] = sim_coords[1] + move
    sim_coords[2] = sim_coords[2] + move
    sim_seq = get_sequence(sim_coords, genome, impose_strand = True, bed_input = True)
    local_window_size = int(window_size - motif_length)
    sim_real_score = struct.get_bp_prob(sim_seq, local_window_size, local_window_size, motif_length, output_folder = output_folder)[last_base, (motif_length - 1)]
    temp_sim_seq = list(sim_seq)
    temp_sim_seq[orig_pos] = orth_base
    temp_sim_seq = "".join(temp_sim_seq)
    sim_sim_score = struct.get_bp_prob(temp_sim_seq, local_window_size, local_window_size, motif_length, output_folder = output_folder)[last_base, (motif_length - 1)]
    return(sim_real_score, sim_sim_score)

def extend_family(blast_results, families, query):
    '''
    Given a gene identifier (query), find all genes that are connected to it
    in the BLAST results (i.e. one is a hit for the other). Add them to the current family and remove
    the relevant lines from the BLAST results.
    '''
    to_add = [i for i in blast_results if query in i]
    blast_results = [i for i in blast_results if query not in i]
    to_add = flatten(to_add)
    families[-1].extend(to_add)
    families[-1] = list(set(families[-1]))
    return(blast_results, families)

def fill_codons(motif_positions, two_fold, four_fold):
    '''
    Given a list of  motif hit positions, transform it to give a list of full codons.
    NB! This works either if you have blocks of contiguous bases (i.e. full motif hits) or fourfold degenerate positions, that is to say, there are never two adjacent positions.
    Anything else and this function will produce nonsense.
    '''
    motif_positions = [int(i) for i in motif_positions]
    positions_number = len(motif_positions)
    to_change_two_fold = []
    to_change_four_fold = []
    final_pos = []
    for i, j in enumerate(motif_positions):
        final_pos.append(j)
        if i == 0 or (j-motif_positions[i-1]) != 1:#if it's the first motif hit base in the sequence or if it's the first one in a block (the previous motif hit base wasn't contiguous)
            #if it's the second base of a codon
            if j%3 == 1:
                #add the preceding base
                del final_pos[-1]
                final_pos.append(j-1)
                final_pos.append(j)
                #if the base just added is a two-fold degenerate site in a Leucine/Arginine codon, store the position
                if j-1 in two_fold:
                    to_change_two_fold.append(j-1)
                #if it is a four-fold degenerate site in a Leu/Arg codon, store the position in another list
                elif j-1 in four_fold:
                    to_change_four_fold.append(j-1)
            #if it's the third base of a codon
            elif j%3 == 2:
                #add the two preceding bases
                del final_pos[-1]
                final_pos.extend([j-2,j-1, j])
                if j-2 in two_fold:
                    to_change_two_fold.append(j-2)
                elif j-2 in four_fold:
                    to_change_four_fold.append(j-2)
        elif (i+1 == positions_number) or (motif_positions[i+1] - j != 1):#if it's the last motif hit base of the sequence or the last one of a contiguous block
            #if it's the first base of a codon
            if j%3 == 0:
                #remove that base
                del final_pos[-1]
            #if it's the second base of a codon
            elif j%3 == 1:
                #remove both that and the preceding base
                del final_pos[-1]
                del final_pos[-1]
    return(final_pos, to_change_two_fold, to_change_four_fold)

def find_families(fasta_file_name, output_prefix):
    '''
    Given a fasta file, group the sequences into paralogous families.
    '''
    blast_results_file_name = "{0}_blast_results".format(output_prefix)
    output_prefix_short = output_prefix.split("/")
    output_prefix_short = output_prefix_short[-1]
    #run a BLAST all against all for the sequences in the fasta file
    blast_all_against_all("{0}_blast_db".format(output_prefix_short), fasta_file_name, blast_results_file_name)
    names, seqs = rw.read_fasta(fasta_file_name)

    #create an empty list for storing the indices of BLAST query - hit pairs to delete
    #open a .csv file containing the results of a BLAST and turn it into a list
    #delete all of the information except the identifiers of queries and hits
    #identify those pairs where the query and the hit come from the same sequence and delete them
    to_delete = []
    with open(blast_results_file_name) as csvfile:
        blast_results = csv.reader(csvfile, delimiter=',')
        blast_results = list(blast_results)
        print("Total number of BLAST hits.")
        print(len(blast_results))
        for i in blast_results:
            del i[2:12]
            if i[0] == i[1]:
                to_delete.append(i)
    print("Elements to delete:")
    print(len(to_delete))
    print("Unique elements to delete:")
    print(len(list(set(flatten(to_delete)))))
    for i in list(reversed(to_delete)):
        blast_results.remove(i)
            
    print("Number of results without self-matches:")
    print(len(blast_results))
    queries = [i for i,j in blast_results]
    print("Number of queries:")
    print(len(queries))
    print("Number of unique queries:")
    print(len(list(set(queries))))
    matches = [j for i,j in blast_results]
    print("Number of matches:")
    print(len(matches))
    print("Number of unique matches:")
    print(len(list(set(matches))))

    print("Genes that don't overlap between queries and matches:")
    for i in list(set(queries)):
        if i not in list(set(matches)):
            print(i)
    for i in list(set(matches)):
        if i not in list(set(queries)):
            print(i)

    #create an empty list for storing the gene families, another for storing the genes
    #that have already been analyzed within a family and a third one for storing all
    #the genes that have been analyzed across all families.
    #create a counter (fcounter) for storing the number of families that have been created
    #while there are query-hit pairs left,
    #add genes seen in the previous family to the list of all genes analyzed and then empty the
    #first list for the next family
    #pick a random query out of the remaining query-hit pairs and create a new family containing
    #just that query. This is now the current family. Increment fcounter by 1.
    #add all genes that are either hits to the current query or that the current query is a hit to
    #into the current family.
    #loop over all the genes in the current family and add everything they match or are a match
    #to into the current family
    #once you've done all the genes in a family, pick a new random query from the query-hit pairs
    #that are left and start a new family with it
    families = []
    added_something = True
    while len(blast_results) > 0:
        seen = []
        current_pair = random.choice(blast_results)
        families.append(current_pair)
        while added_something:
            length_before = len(families[-1])
            for query in families[-1]:
                if query not in seen:
                    seen.append(query)
                    [blast_results, families] = extend_family(blast_results, families, query)
            if(len(families[-1])) == length_before:
                added_something == False
                break
        
    families_file_name = "{0}_families.txt".format(output_prefix)
    families_file = open(families_file_name,"w")
    for i in range(0,len(families)):
        families_file.write("{0}\n".format(",".join(families[i])))

    #create flat version of the families list so you could count the total number of genes that have been allocated to a family
    flat_families = flatten(families)

    #these two numbers should be identical
    print("Number of genes in families:")
    print(len(flat_families))
    print("Number of unique genes in families:")
    print(len(list(set(flat_families))))

    #create a list with the sizes of all the different families
    family_sizes = [len(i) for i in families]
    print("Number of families:")
    print(len(families))
    print("Distribution of family sizes:")
    print(sorted(family_sizes))

    #close the output file
    families_file.close()

def get_aligned_nt_sequence_from_prot(nt_sequence, aligned_prot_sequence):
    '''
    input: a) a nucleotide sequence b) an aligned protein sequence
    output: the aligned nucleotide sequence corresponding to the protein sequence in b) above
    turn the nucleotide sequence into a list
    turn the protein sequence into a list
    initialize a counter that would allow you to traverse the nucleotide sequence codon by codon
    create a list for the sequence that will be output
    loop over the protein sequence
    if you come across a gap, insert three gaps into the sequence you're building
    if you come across anything else, check that the codon and the amino acid match up and insert the next codon in the input nucleotide sequence
    increment the codon counter by three
    '''
    temp_sequence = list(nt_sequence)
    temp_prot_sequence = list(aligned_prot_sequence)
    codon_counter = 0
    aligned_sequence = []
    for i in temp_prot_sequence:
        if i == "-":
            aligned_sequence.extend(["-","-","-"])
        else:
            current_codon = "".join(temp_sequence[codon_counter: codon_counter + 3])
            if (current_codon not in ["TAA","TAG","TGA"]) and (i != codon_dict[current_codon]):
                print(codon_dict)
                print(current_codon)
                print(i)
                print(nt_sequence)
                print(aligned_prot_sequence)
                print("There is a problem with the translation.")
                raise Exception
            aligned_sequence.extend(temp_sequence[codon_counter: codon_counter + 3])
            codon_counter = codon_counter + 3
    return(aligned_sequence)

##def get_aligned_positions(aligned_sequences, positions):
##    '''
##    Given two aligned DNA sequences and indices relative to an unaligned version of the first sequence,
##    convert the indices into indices relative to the aligned version of the first sequence.
##    '''
##    positions = np.array(positions)
##    if "-" in aligned_sequences[0]:#no need to make adjustments if there are no indels in the alignment
##        aligned_sequences[0] = np.array(aligned_sequences[0])
##        hyphens = np.where(aligned_sequences[0] == "-")[0] #where are the dashes in the alignment?
##        counter = 0#keeps track of by how much the current motif position needs to be shifted to the right
##        no_more_hyphens = False#goes to True when all indel positions have been processed
##        pos_index = 0#keeps track of the motif position that is currently being processed
##        repeat = False
##        while pos_index < len(positions):#until all motif positions have been processed
##            local_counter = 0#keeps track of the shift due to the current contiguous indel codon block.
##            if not no_more_hyphens:#until all indels have been taken into account
##                while hyphens[0] <= positions[pos_index]:#if any dashes precede or coincide with the first motif position (unaligned), thus causing a shift
##                    local_counter = local_counter + 1 #update this counter to account for the motif position shifting to the right
##                    contig_counter = 0#counts the number of contiguous codon indels
##                    if len(hyphens) != 1:#if the current indel codon is not the last one
##                        while hyphens[contig_counter+1] - hyphens[contig_counter] == 1:#if the following hyphen is next to the current one
##                            local_counter = local_counter + 1#update counter to account for the motif position shifting downstream
##                            contig_counter = contig_counter + 1#update counter to show you have another contiguous indel codon
##                            if contig_counter+1 == len(hyphens):#if the current indel codon is the last one, break out of the loop to avoid getting an error when it tries to check the next indel codon for contiguity.
##                                break
##                        hyphens = np.delete(hyphens, range(contig_counter))#once you've either processed all indel codons or the next indel codon is not contiguous to the current one, delete all the indel codons in the current contiguous block.
##                    hyphens = np.delete(hyphens, 0)#because contig counter gives you the number of contiguous codons - 1 so even if you've deleted 0:contig_counter, you still need to delete one more indel codon.
##                    if len(hyphens) < 1:
##                        no_more_hyphens = True
##                        break
##            counter = counter + local_counter#add the shift due to the current contiguous block of indel codons (possibly just one codon)
##            if repeat:#if this is not the only block of indel codons to precede the current motif position, only add the shift due to the current block
##                positions[pos_index] = positions[pos_index] + local_counter
##            else:#otherwise shift the motif position by the total shift that has accumulated
##                positions[pos_index] = positions[pos_index] + counter
##            if len(hyphens) > 0:#if all indel codons haven't been deleted
##                if hyphens[0] <= positions[pos_index]:#if there is another indel codon block in front of or coincident with the current motif position
##                    repeat = True #go through the loop with the next indel codon block without mocing to the next motif position
##                else:
##                    pos_index = pos_index + 1 #move to the next motif position
##                    repeat = False
##            else:#otherwise move to the next motif position and then it will shift the motif position but will skip the bit where it increases the shift (because no_more_hyphens is True)
##                pos_index = pos_index + 1
##                repeat = False
##    positions = sorted(list(set(positions)))
####    positions = [i for i in positions if i < len(aligned_sequences[0])]#to remove motif positions that overlap with the stop codon
##    excess_counter = 0
##    for pos in range(len(positions) - 1, 0, -3):
##        if positions[pos] >= len(aligned_sequences[0]) or positions[pos - 1] >= len(aligned_sequences[0]) or positions[pos - 2] >= len(aligned_sequences[0]):
##            excess_counter = excess_counter + 1
##        else:
##            break
##    for i in range(excess_counter):
##        del positions[-3:]
##    return(positions)

def get_aligned_positions(aligned_sequences, positions):
    '''
    Given two aligned DNA sequences and indices relative to an unaligned version of the first sequence,
    convert the indices into indices relative to the aligned version of the first sequence.
    '''
    if len(positions) > 0 and "-" in aligned_sequences[0]:
        max_pos = max(positions)
        corresp_dict = {}
        counter = 0
        seq = aligned_sequences[0]
        for position, base in enumerate(seq):
            if base == "-":
                pass
            else:
                corresp_dict[counter] = position
                counter = counter + 1
                if counter > max(positions):
                    break
        #if i is not in corresp_dict, that probably means that it overlaps a stop
        out_positions = [corresp_dict[i] for i in positions if i in corresp_dict]
        return(out_positions)
    else:
        return(positions)

def get_alignment(trans, alignment_folder, correspondances):
    '''
    Get the pairwise alignment corresponding to a particular transcript.
    '''
    orth_id = correspondances[trans]
    phy_file_name = "{0}/{1}_{2}.phy".format(alignment_folder, trans, orth_id)
    aligned_sequences = [str(i.seq) for i in SeqIO.parse(phy_file_name, "phylip-sequential")]
    return(aligned_sequences)

def get_LeuArg(sequence, two_or_four):
    '''
    Get the positions of first bases of Leucine and Arginine codons. Store two-fold and four-fold degenerate ones separately.
    '''
    positions = []
    if two_or_four == 2:
        codons = twofold_LeuArg
    elif two_or_four == 4:
        codons = fourfold_LeuArg
    for i in range(0, len(sequence), 3):
        if sequence[i:i+3] in codons:
            positions.append(i)
    return(positions)

def get_motif_potential(sequence, motifs, position, motif_lengths, reduction = False):
    '''
    Given a sequence (a list), a set of motifs and a position in that sequence,
    calculate what fraction of all the possible mutations at that position would increase the
    motif density of the sequence.
    '''
    all_bases = ["A", "T", "C", "G"]
    current_base = sequence[position]
    all_bases.remove(current_base)
    true_overlap = nc.get_motif_set_density(motifs, motif_lengths, sequence, concat = True)["count"]
    sequence = list(sequence)
    counter = 0
    danger_bases = []
    for base in all_bases:
        new_sequence = sequence.copy()
        new_sequence[position] = base
        new_sequence = "".join(new_sequence)
        current_overlap = nc.get_motif_set_density(motifs, motif_lengths, new_sequence, concat = True)["count"]
        if reduction:
            if current_overlap < true_overlap:
                counter = counter + 1
                danger_bases.append(base)
        else:
            if current_overlap > true_overlap:
                counter = counter + 1
                danger_bases.append(base)
    if reduction:
        return(counter)
    #because there are always three possible mutations and you want to see what fraction increase the motif density
    fraction = counter/4
    return(fraction, danger_bases)

def get_relative_SNPs(CDSs, SNP_file_name, CDS_SNP_file_name, seqs, names, genome, get_new_SNPs = False, parse_SNPs = False, remove_GT = False):
    if get_new_SNPs:
        flat_CDS = [[j[0] for j in i] for i in list(CDSs.values())]
        CDS_bed = "temp_data/temp_{0}.bed".format(random.random())
        write_features_to_bed(flat_CDS, CDS_bed, modify_chr_ids = True)
        tabix(CDS_bed, SNP_file_name, genome)
        os.remove(CDS_bed)
    if parse_SNPs:
        relative_coords_dict = {i: {"positions": [], "to_remove": [], "var_bases": [], "MAFs": [], "allele_counts": []} for i in list(CDSs.keys())}
        counter = 0
        with open(SNP_file_name) as file:
            error_counter = 0
            #loop over bed file with SNPs
            for line in file:
                if line[0] != "#":
                    if counter % 1000 == 0:
                        print(counter)
                    counter = counter + 1
                    bed_record = line.split("\t")
                    info = bed_record[6].split("$")
                    #reported ancestral base
                    ref_base = info[1]
                    #filter out anything that isn't simple SNPs
                    var_base = info[2].split(",")
                    var_base_number = len(var_base)
                    var_base = [i for i in var_base if i in nc._canon_bases_]
                    trans = bed_record[3]
                    #convert bed record to indices relative to the CDS
                    #mismatch_warning_only = True means that if the SNP coordinates don't map onto the CDS,
                    #return an error message but don't crash.
                    current_index = bed_to_CDS_indices(bed_record, CDSs[trans], mismatch_warning_only = True)[0]
                    #this is so you could remove the stuff you filter out for various reasons also from the monomorphic set
                    relative_coords_dict[trans]["to_remove"].append(str(current_index))
                    #the SNP bed files contain information on which transcript they overlap
                    if current_index != ("error"):
                        if len(ref_base) == 1:
                            #filtering out polymorphisms with more than 2 segregating alleles
                            #note that you need to check the number both before and after filtering out non-canonical bases
                            if var_base_number == 1 and len(var_base) == 1:
                            #check if the reported ancestral allele matches the base that is at that position in the CDS
                            #it's expected that sometimes it won't but it should most of the time if everything's worked out properly
                                strand = CDSs[trans][0][0][6]
                                var_base = var_base[0]
                                if strand == "-":
                                    ref_base = nc.rev_comp(ref_base)
                                    var_base = nc.rev_comp(var_base)
                                CDS_base = seqs[names.index(trans)][current_index]
                                if ref_base != CDS_base:
                                    print("PROBLEM!")
                                    print(ref_base)
                                    print(CDS_base)
                                    print(var_base)
                                    print("\n")
                                elif remove_GT and ref_base == "G" and var_base == "T":
                                    pass
                                else:
                                    relative_coords_dict[trans]["positions"].append(str(current_index))
                                    relative_coords_dict[trans]["var_bases"].append(str(var_base))
                                    MAF_MAC = info[4].split(";")
                                    MAC = MAF_MAC[0].split("=")[1]
                                    MAF = MAF_MAC[1].split("=")[1]
                                    relative_coords_dict[trans]["MAFs"].append(str(MAF))
                                    relative_coords_dict[trans]["allele_counts"].append(str(MAC))
                    else:
                        error_counter = error_counter + 1
        with open(CDS_SNP_file_name, "w") as file:
            for trans in relative_coords_dict:
                to_write = zip(relative_coords_dict[trans]["positions"], relative_coords_dict[trans]["var_bases"], relative_coords_dict[trans]["MAFs"], relative_coords_dict[trans]["allele_counts"])
                to_write = "\t".join([trans, "|".join([",".join(i) for i in to_write]), ",".join(relative_coords_dict[trans]["to_remove"])])
                file.write(to_write)
                file.write("\n")
        relative_coords_dict = {}
        print("Number of conversion errors:")
        print(error_counter)
    
    current_SNPs = rw.read_many_fields(CDS_SNP_file_name, "\t")
    current_SNPs_to_remove = list_to_dict(current_SNPs, 0, 2)
    current_SNPs_to_remove = {i: [int(j) for j in current_SNPs_to_remove[i].split(",") if j != "error"] for i in current_SNPs_to_remove if current_SNPs_to_remove[i]}
    current_SNPs = list_to_dict(current_SNPs, 0, 1)
    current_SNPs = {i: [j.split(",") for j in current_SNPs[i].split("|")] for i in current_SNPs if current_SNPs[i]}
    current_SNPs = {i: {int(j[0]): (j[1], float(j[2])) for j in current_SNPs[i]} for i in current_SNPs}
    return(current_SNPs, current_SNPs_to_remove)

def get_SNP_density(picked, motifs, seqs, names, SNP_dict, map_from_regions):
    '''
    Given sequences, motifs and SNP positions, calculate the fraction of motif bases at fourfold degenerate sites that overlap with SNPs.
    '''
    motif_lengths = [len(i) for i in motifs]
    motifs = nc.motif_to_regex(motifs)
    if not picked:
        picked = names.copy()
    SNP_dens_parallel = run_in_parallel(picked, ["foo", seqs, names, motifs, motif_lengths, SNP_dict, map_from_regions], get_SNP_density_core)
    SNP_counter = 0
    site_counter = 0
    for result in SNP_dens_parallel:
        current = result.get()
        SNP_counter = SNP_counter + current[0]
        site_counter = site_counter + current[1]
    if site_counter > 0:
        real_fraction = SNP_counter/site_counter
    else:
        real_fraction = None
    return(real_fraction)

def get_SNP_density_core(picked, seqs, names, motifs, motif_lengths, SNP_dict, map_from_regions):
    '''
    Core for get_SNP_density.
    '''
    site_counter = 0
    SNP_counter = 0
    for trans in picked:
        if map_from_regions:
            SNP_idn = map_from_regions[trans]["idn"]
        else:
            SNP_idn = trans
        if "NA" not in SNP_dict[SNP_idn]:
            seq = seqs[names.index(trans)]
            motif_pos = nc.get_motif_set_density(motifs, motif_lengths, seq, concat = True)["positions"]
            fourfold = nc.get_4fold_deg(seq)
            fourfold = [i for i in fourfold if (seq[i] != "C") and (seq[i] != "G")]
            motif_pos = [i for i in motif_pos if i in fourfold]
            site_counter = site_counter + len(motif_pos)
            if map_from_regions:
                motif_pos = region_indices_to_full_indices(motif_pos, map_from_regions[trans]["flank indices in CDS"])
            motif_pos = [i for i in motif_pos if i in SNP_dict[SNP_idn]]
            SNP_counter = SNP_counter + len(motif_pos)
    return(SNP_counter, site_counter)

def get_unaligned_pos_dict(seq):
    '''
    Given a sequence from an alignment, create a dictionary where the keys are positions in the aligned sequence
    and the values are the corresponding positions in the unaligned sequence.
    '''
    output = {}
    counter = 0
    for pos, base in enumerate(seq):
        if base == "-":
            output[pos] = None
        else:
            output[pos] = counter
            counter = counter + 1
    return(output)
    
def input_dict_for_dS(correspondances_file_name, alignment_folder_name, fasta_file_name, output_file_name, map_from_regions = False, picked = None):
    '''
    Prepare the necessary input data for dS_from_hits and write it to a file.
    '''
    if not map_from_regions:
        names, seqs = rw.read_fasta(fasta_file_name)
    else:
        names = []
        seqs = []
        for fasta in map_from_regions["fastas"]:
            curr_names, curr_seqs = rw.read_fasta(fasta)
            names.append(curr_names)
            seqs.append(curr_seqs)
        #make two dictionaries, one that maps from the fasta file names (presumably coordinates)
        #to gene identifiers and one that maps to transcript identifiers
        #the entries in the regions bed file have to have the same order as the entries in the fasta file!
        regions_bed = rw.read_many_fields(map_from_regions["regions bed file"][0], "\t")
        mapping_to_gene_ids = {}
        mapping_to_trans_ids = {}
        for pos, name in enumerate(names[0]):
            trans_id = regions_bed[pos][3]
            mapping_to_trans_ids[name] = trans_id
            mapping_to_gene_ids[name] = [key for key in map_from_regions["gene name dict"] if map_from_regions["gene name dict"][key][0] == trans_id][0]
    two_seqs = False
    if "|" in seqs[0]:
        two_seqs = True
        names = [i.replace("|", "&") for i in names]
        names = [i.replace("%", "") for i in names]
    else:
        correspondances = rw.read_many_fields(correspondances_file_name, ",")
        correspondance_dict = {}
        for i in correspondances:
            correspondance_dict[i[0]] = i[1]
    if not picked:
        if map_from_regions:
            picked = names[0]
        else:
            picked = names.copy()
    CG = ["GC", "CG"]
    motif_lengths = [2, 2]
    CG = nc.motif_to_regex(CG)
    with open(output_file_name, "w") as file:
        for idn in picked:
            current_list = [idn]
            if two_seqs:
                aligned_sequences = seqs[names.index(idn)]
                current_list.append("{0}_orth".format(idn))
            else:
                if map_from_regions:
                    orth_idn = correspondance_dict[mapping_to_gene_ids[idn]]
                    if "_" not in orth_idn:
                        orth_idn = orth_idn + "_0"
                    phy_file_name = "{0}/{1}_{2}.phy".format(alignment_folder_name, mapping_to_gene_ids[idn], orth_idn)
                else:
                    orth_idn = correspondance_dict[idn]
                    if "_" not in orth_idn:
                        orth_idn = orth_idn + "_0"
                    phy_file_name = "{0}/{1}_{2}.phy".format(alignment_folder_name, idn, orth_idn)
                current_list.append(orth_idn)
                try:
                    aligned_sequences = [str(i.seq) for i in SeqIO.parse(phy_file_name, "phylip-sequential")]
                except FileNotFoundError:
                    phy_file_name = phy_file_name[:-6] + ".phy"
                    current_list[-1] = current_list[-1][:-2]
                    aligned_sequences = [str(i.seq) for i in SeqIO.parse(phy_file_name, "phylip-sequential")]
                aligned_sequences = "|".join(aligned_sequences)
            current_list.append(aligned_sequences)
            #it's called "full_CDS" here but if the fasta file contains regions, it'll actually be a region (flank/core etc.)
            if two_seqs:
                temp = aligned_sequences.split("|")
                full_CDS = "".join([i for i in temp[0] if i in nc._canon_bases_])
                two_fold = "NA"
                four_fold = "NA"
            else:
                if not map_from_regions:
                    full_CDS = seqs[names.index(idn)]
                    #where are the positions of the first bases of twofold/fourfold Leucine/Arginine codons?
                    two_fold = ",".join([str(i) for i in get_LeuArg(list(full_CDS), 2)])
                    four_fold = ",".join([str(i) for i in get_LeuArg(list(full_CDS), 4)])
                else:
                    full_CDS = [seq[names[0].index(idn)] for seq in seqs]
                    two_fold = "|".join([",".join([str(i) for i in get_LeuArg(list(seqs), 2)]) for seq in full_CDS])
                    four_fold = "|".join([",".join([str(i) for i in get_LeuArg(list(seqs), 4)]) for seq in full_CDS])
                    full_CDS = "|".join(full_CDS)                   
            current_list.append(full_CDS)
            current_list.append(two_fold)
            current_list.append(four_fold)
            CG_pos = nc.get_motif_set_density(CG, motif_lengths, full_CDS, concat = True)["positions"]
            CG_pos = [range(i - (i%3), i - (i%3) + 3) for i in CG_pos]
            CG_pos = flatten(CG_pos)
            CG_pos = sorted(list(set(CG_pos)))
            CG_pos = ",".join([str(i) for i in CG_pos])
            current_list.append(CG_pos)
            all_4f_pos = nc.get_4fold_deg(full_CDS)
            all_4f_pos = ",".join([str(i) for i in all_4f_pos])
            current_list.append(all_4f_pos)
            if map_from_regions:
                flank_indices_in_CDS = []
                first_regions_bed = rw.read_many_fields(map_from_regions["regions bed file"][0], "\t")
                for pos, region_file in enumerate(map_from_regions["regions bed file"]):
                    if pos == 0:
                        split_idn = re.split("[:\(\)]", idn)
                        coords = split_idn[1].split("-")
                        current_coords = [[split_idn[0], coords[0], coords[1]]]
                        current_record = intersect_bed_return_bed(region_file, current_coords, bed_input = True)
                        if len(current_record) > 1:
                            print("Can't match up fasta record to bed.")
                            raise Exception
                        current_record = current_record[0]
                        current_location = first_regions_bed.index(current_record)
                    else:
                        current_record = rw.read_many_fields(region_file, "\t")[current_location]
                    current_flank_indices_in_CDS = ",".join([str(i) for i in bed_to_CDS_indices(current_record, map_from_regions["CDS"][mapping_to_trans_ids[idn]])])
                    flank_indices_in_CDS.append(current_flank_indices_in_CDS)
                current_list.append("|".join(flank_indices_in_CDS))
            current_list = "%".join(current_list)
            file.write(current_list)
            file.write("\n")

def INSIGHT(neutral_file, hit_file, freq_threshold, INSIGHT_dir, data_id):
    run_process(["cp", neutral_file, INSIGHT_dir])
    run_process(["cp", hit_file, INSIGHT_dir])
    current_dir = os.getcwd()
    os.chdir("../Software/INSIGHT")
    remove_file("{0}.ins.log".format(data_id))
    remove_file("{0}.ins.results.txt".format(data_id))
    remove_file("{0}.post.sites".format(data_id))
    run_process(["bash", "scripts/runINSIGHT-EM.sh", data_id, hit_file.split("/")[-1], neutral_file.split("/")[-1], ".", int(freq_threshold * 100)])
    os.chdir(current_dir)

def keep_conserved_pc(trans_id, orth_ids, CDS, orth_CDS, dS_threshold, alignments_folder_name):
    '''
    Given a CDS, check whether it has an ortholog in the other species to which it aligns with dS below dS_threshold and omega below 0.5 (though it doesn't have to be the same ortholog for
    the two conditions).
    '''
    
    CDS_IUPAC = Seq(CDS, IUPAC.unambiguous_dna)
    CDS_prot = CDS_IUPAC.translate()
    orth_CDS_IUPACs = [Seq(i, IUPAC.unambiguous_dna) for i in orth_CDS]
    orth_CDS_prots = [i.translate() for i in orth_CDS_IUPACs]   

    dS = []
    omegas = []
    #loop over all the orthologs in the other species
    for i in range(len(orth_CDS_prots)):
        #write the focal protein sequence and the orthologous proteins sequence to a phylip file
        prot_alignment_file_name = "temp_data/prot_alignment_file{0}.fasta".format(random.random())
        with open(prot_alignment_file_name, "w") as file:
            file.write(">{0}\n".format(trans_id))
            file.write(str(CDS_prot))
            file.write("\n>{0}\n".format(orth_ids[i]))
            file.write(str(orth_CDS_prots[i]))
        #align the two sequences using MUSCLE (protein alignment)
        input_path = "/Users/{0}/Documents/Scripts_and_data/{1}".format(user, prot_alignment_file_name)
        output_file_name = "temp_data/output_file{0}.fasta".format(random.random())
        output_path = "/Users/{0}/Documents/Scripts_and_data/{1}".format(user, output_file_name)
        muscle_object = MuscleCommandline(muscle_exe, input = input_path, out = output_path)
        stdout, stderr = muscle_object()
        remove_file(prot_alignment_file_name)
        file_from_muscle = open(output_file_name)
        muscle_string = "".join(file_from_muscle)
        file_from_muscle.close()
        remove_file(output_file_name)
        #in case the sequence is broken across multiple lines
        muscle_string = re.sub("([A-Z\-])\n([A-Z\-])","\\1\\2", muscle_string)
        aligned_prot_sequences = re.findall("^[A-Z\-]+(?=\n)", muscle_string, re.MULTILINE)
        #convert the protein alignment to a nucleotide alignment
        aligned_sequences = [[] for l in range(2)]
        aligned_sequences[0] = get_aligned_nt_sequence_from_prot(CDS, aligned_prot_sequences[0])
        aligned_sequences[1] = get_aligned_nt_sequence_from_prot(orth_CDS[i], aligned_prot_sequences[1])
        aligned_sequences_IUPAC = [Seq("".join(j),IUPAC.unambiguous_dna) for j in aligned_sequences]
        #write the nucleotide alignment to a phylip file
        alignment = MultipleSeqAlignment([SeqRecord(aligned_sequences_IUPAC[0], id = "seq1"),SeqRecord(aligned_sequences_IUPAC[1], id = "seq2")])
        phy_file_name = "{0}/{1}_{2}.phy".format(alignments_folder_name, trans_id, orth_ids[i])
        AlignIO.write(alignment, phy_file_name,"phylip-sequential")
        #run PAML codeml on the sequence, get dS and omega
        cml_input = "{0}/{1}".format(cwd, phy_file_name)
        cml_output = "{0}/phy_file{1}.out".format(cwd, random.random())
        cml_tree = "general/human_macaque.tree"
        cml = codeml.Codeml(alignment = cml_input, out_file = cml_output, tree = cml_tree)
        cml.set_options(verbose = True, seqtype = 1, runmode = 0, model = 0, NSsites = [])
        temp_dict = cml.run(command = "/Users/{0}/Documents/Software/paml4.8/bin/codeml".format(user))
        remove_file(cml_output)
        dS.append(temp_dict["NSsites"][0]["parameters"]["dS"])
        omegas.append(temp_dict["NSsites"][0]["parameters"]["omega"])
    min_dS = min(dS)
    min_omega = min(omegas)
    #if all the dSs/omegas are too high
    if min_dS >= dS_threshold or min_omega >= 0.5:
        return(False, None)
    #otherwise return the ortholog that gave the lowest dS
    else:
        best_match = orth_ids[dS.index(min_dS)]
        return(True, best_match)

def map_regions_to_CDS(fasta, bed, fs, transcripts, CDS, trans_ids = False):
    '''
    Given a fasta file, a bed file and a set of CDS coordinates, return for each record the relative coordinates to which
    the bed record maps in the corresponding full CDS.
    '''
    names, seqs = rw.read_fasta(fasta)
    regions_bed = rw.read_many_fields(bed, "\t")
    mapping_to_gene_ids = {}
    for pos, name in enumerate(names):
        mapping_to_gene_ids[name] = {}
        trans_id = regions_bed[pos][3]
        if trans_ids:
            mapping_to_gene_ids[name]["idn"] = trans_id
        else:
            gene_id = fs.convert_between_ENST_and_ENSG(trans_id, transcripts, "ENSG")
            mapping_to_gene_ids[name]["idn"] = gene_id
        current_record = regions_bed[pos]
        current_flank_indices_in_CDS = bed_to_CDS_indices(current_record, CDS[trans_id])
        mapping_to_gene_ids[name]["flank indices in CDS"] = current_flank_indices_in_CDS
    return(mapping_to_gene_ids)

def MSA_motif_cons(input_file, motifs, output_folder, lower_threshold, upper_threshold):
    motif_lengths = [len(i) for i in motifs]
    motifs = nc.motif_to_regex(motifs)
    make_dir(output_folder)
    unknown = [".", "-"]
    human = "homo_sapiens"
    output_dict = {i: [] for i in range(lower_threshold, upper_threshold + 1)}
    lines = line_count(input_file)
    print("Total: {0}.".format(lines))
    with open(input_file, "r") as file:
        for pos, line in enumerate(file):
            if pos % 1000 == 0:
                print(pos)
            #there's an empty line at the end
            if line:
                #to get rid of stuff trailing at the end
                line = line.rstrip("%\n")
                line = [i.split("|") for i in line.split("%")]
                line = [i for i in line if i[-1] == "1"]
                if len(line) >= lower_threshold:
                    dashed = {i[0]: i[2] for i in line}
                    dashless = {i: "".join([j for j in dashed[i] if j not in unknown]) for i in dashed}
                    #the fourfold deg. sites need to be calculated with regards to a trimmed exon
                    human_phase = int(line[0][-2])
                    trim_info, trimmed_human = nc.trim_sequence_report(dashless[human], human_phase)
                    fourfold_deg = nc.get_4fold_deg(trimmed_human)
                    if fourfold_deg:
                        #but then you need to know where the 4fold degs sit in the original sequence
                        fourfold_deg = [i + trim_info[0] for i in fourfold_deg]
                        human_pos = nc.get_motif_set_density(motifs, motif_lengths, dashless[human], concat = True)["positions"]
                        human_pos = [i for i in human_pos if i in fourfold_deg]
                        if human_pos:
                            human_pos = get_aligned_positions([list(dashed[human]), []], human_pos)
                            all_pos = {i: nc.get_motif_set_density(motifs, motif_lengths, dashless[i], concat = True)["positions"] for i in dashless if i != human}
                            all_pos = {i: get_aligned_positions([list(dashed[i]), []], all_pos[i]) for i in all_pos}
                            for position in human_pos:
                                site_counter = 1
                                hit_counter = 1
                                for species in all_pos:
                                    if (dashed[species][position] != "N") and (dashed[species][position] not in unknown):
                                        site_counter = site_counter + 1
                                        if position in all_pos[species]:
                                            hit_counter = hit_counter + 1
                                if site_counter >= lower_threshold:
                                    output_dict[site_counter].append(hit_counter)
    for site_number in output_dict:
        with open("{0}/{1}.csv".format(output_folder, site_number), "w") as out_file:
            for count in output_dict[site_number]:
                out_file.write(str(count))
                out_file.write("\n")

def mutation_to_motif(current_CDS, aligned_sequences, motifs, neighbours, neighbour_lengths, mutation_score, site_number):
    '''
    Given a CDS sequence, an alignment of that CDS to an orthologous CDS, a set of motifs and the set of motifs that are a single
    base substitution away from the first set of motifs, find fourfold degenerate sites within your focal species that are a single substitution
    away from one of your motifs and check whether the base that the orhologous species has at that position would generate
    a motif in human.
    '''
    bases = ["A", "T", "C", "G"]
    #you want the motif hit positions hit by hit though not motif by motif (positions of motifs that are a single base substitution away from a motif
    #in you rmain set)
    current_positions = flatten(nc.get_motif_set_density(neighbours, neighbour_lengths, current_CDS, concat = True, raw = True))
    #get the fourfold degenerate sites
    current_4f = nc.get_4fold_deg(current_CDS)
    #loop over the fourfold degenerate sites
    for syn_site in current_4f:
        #leave only motif hits that overlap with the current position
        current_hits = [i for i in current_positions if syn_site in i]
        if current_hits:
            #make a uniquified and flattened version of the motif hit positions list
            flat_positions = np.concatenate(tuple(current_positions))
            flat_positions = np.unique(flat_positions)
            aligned_sequences = [list(i) for i in aligned_sequences]
            #convert the flat positions list from indices relative to the unaligned CDS to indices relative to the aligned CDS
            aligned_positions = get_aligned_positions(aligned_sequences, flat_positions)
            #extract the motifs from the unflattened motif hit positions
            current_motifs = [[current_CDS[i] for i in j] for j in current_hits]
            #get the base that overlaps the current fourfold-degenerate site, as well as the 3 other bases
            current_base = current_CDS[syn_site]
            other_bases = bases.copy()
            other_bases.remove(current_base)
            potential_counter = 0
            danger_bases = []
            #for each of the 3 other bases, try substituting them in at the current fourfold degenerate site,
            #in all the motifs that overlap
            for base in other_bases:
                found = False
                for pos, motif in enumerate(current_motifs):
                    current_position = np.where(current_hits[pos] == syn_site)[0]
                    new_motif = motif.copy()
                    new_motif[current_position] = base
                    if new_motif in motifs:
                        found = True
                        danger_bases.append(base)
                #count how many of the 3 possible substitutions generate one of the motifs in the set
                if found:
                    potential_counter = potential_counter + 1
            motif_potential = potential_counter/4
            #if at least one possible substitution would generate a motif
            if motif_potential > 0:
                site_number = site_number + 1
                #check the base the orthologous species has at that position
                orth_base = aligned_sequences[1][aligned_positions[np.where(flat_positions == syn_site)[0]]]
                #if it's one of the bases that would generate a motif in human
                if orth_base in danger_bases:
                    #score the site between 1/4 and 3/4 depending on how many of the 3 substitutions would generate a motif
                    mutation_score = mutation_score + (1 - motif_potential)
    return(mutation_score, site_number)

def protein_alignment(fasta_id1, fasta_seq1, fasta_id2, fasta_seq2, phylip_id1, phylip_id2, phy_file_name):
    '''
    Given two CDS sequences, perform a MUSCLE protein alignment on them, reconvert the sequence back to nucleotides and write to a PHYLIP file.
    '''
    #has to be a sequence object or I can't use BioPython translate()
    fasta_seq1_IUPAC = Seq(fasta_seq1, IUPAC.unambiguous_dna)
    fasta_seq2_IUPAC = Seq(fasta_seq2, IUPAC.unambiguous_dna)
    #translate DNA to protein
    prot1 = fasta_seq1_IUPAC.translate()
    prot2 = fasta_seq2_IUPAC.translate()
    #write protein sequences to fasta
    prot_alignment_file_name = "temp_data/prot_alignment_file{0}.fasta".format(random.random())
    with open(prot_alignment_file_name, "w") as file:
        file.write(">{0}\n".format(fasta_id1))
        file.write(str(prot1))
        file.write("\n>{0}\n".format(fasta_id2))
        file.write(str(prot2))
    #align using MUSCLE
    input_path = "/Users/{0}/Documents/Scripts_and_data/{1}".format(user, prot_alignment_file_name)
    output_file_name = "temp_data/output_file{0}.fasta".format(random.random())
    output_path = "/Users/{0}/Documents/Scripts_and_data/{1}".format(user, output_file_name)
    muscle_object = MuscleCommandline(muscle_exe, input=input_path, out = output_path)
    stdout, stderr = muscle_object()
    remove_file(prot_alignment_file_name)
    file_from_muscle = open(output_file_name)
    muscle_string = "".join(file_from_muscle)
    file_from_muscle.close()
    remove_file(output_file_name)
    muscle_string = re.sub("([A-Z\-])\n([A-Z\-])","\\1\\2", muscle_string)
    aligned_prot_sequences = re.findall("^[A-Z\-]+(?=\n)",muscle_string, re.MULTILINE)
    #convert back to DNA
    aligned_sequences = [[] for l in range(2)]
    aligned_sequences[0] = get_aligned_nt_sequence_from_prot(fasta_seq1,aligned_prot_sequences[0])
    aligned_sequences[1] = get_aligned_nt_sequence_from_prot(fasta_seq2,aligned_prot_sequences[1])
    aligned_sequences_IUPAC = [Seq("".join(j),IUPAC.unambiguous_dna) for j in aligned_sequences]
    #remove stops
    aligned_sequences_IUPAC = [j[:-3] for j in aligned_sequences_IUPAC]
    #write to PHYLIP file
    alignment = MultipleSeqAlignment([SeqRecord(aligned_sequences_IUPAC[0], id = phylip_id1),SeqRecord(aligned_sequences_IUPAC[1], id = phylip_id2)])
    AlignIO.write(alignment, phy_file_name,"phylip-sequential")

def parse_degen(file_name):
    degen = rw.read_many_fields(file_name, "\t")
    degen = list_to_dict(degen, 0, 1)
    degen = {i: degen[i].split(",") for i in degen}
    for trans in degen:
        separate = [i.split(":") for i in degen[trans]]
        separate = [i for i in separate if len(i) == 2]
        degen[trans] = {int(i[0]): i[1].split("|") for i in separate}
    return(degen)

def parse_input_dict_line(line):
    '''
    Parse a line from input_dict_from_dS into a dictionary.
    '''
    output_dict = {}
    line = line.rstrip("\n")
    line = line.split("%")
    idn = line[0]
    output_dict["orth id"] = line[1]
    output_dict["aligned sequences"] = [list(i) for i in line[2].split("|")]
    if len(line) == 9:
        output_dict["flank indices in CDS"] = {}
        output_dict["full CDS"] = {}
        output_dict["two-fold"] = {}
        output_dict["four-fold"] = {}
        CG = line[6].split(",")
        all_fourfold = line[7].split(",")
        flank_indices = line[8].split("|")
        full_CDS = line[3].split("|")
        two_fold = line[4].split("|")
        four_fold = line[5].split("|")
        for pos, region in enumerate(flank_indices):
            output_dict["full CDS"][pos] = full_CDS[pos]
            output_dict["two-fold"][pos] = [int(i) for i in two_fold[pos].split(",") if len(i) > 0]
            output_dict["four-fold"][pos] = [int(i) for i in four_fold[pos].split(",") if len(i) > 0]
            output_dict["flank indices in CDS"][pos] = [int(i) for i in region.split(",")]
    else:
        CG = line[6].split(",")
        all_fourfold = line[7].split(",")
        output_dict["full CDS"] = line[3]
        if line[4] != "NA":
            output_dict["two-fold"] = [int(i) for i in line[4].split(",") if len(i) > 0]
            output_dict["four-fold"] = [int(i) for i in line[5].split(",") if len(i) > 0]
        else:
            output_dict["two-fold"] = "NA"
            output_dict["four-fold"] = "NA"
    output_dict["CG"] = [int(i) for i in CG if len(i) > 0]
    output_dict["all_fourfold"] = [int(i) for i in all_fourfold if len(i) > 0]
    return(idn, output_dict)

def reencode_disruption_pairwise(alignment, positions, degen):
    '''Take a pairwise alignment and re-encode in such a way that only motif-disrupting divergences are preserved.'''
    stops = [["T", "G", "A"], ["T", "A", "A"], ["T", "A", "G"]]
    to_remove = []
    for site in positions:
        human_base = alignment[0][site]
        alt_base = alignment[1][site]
        if alt_base != "-" and alt_base != human_base:
            if alt_base not in degen[site]:
                alignment[1][site] = human_base
                #this is because if you've just introduced a PTC, PAML will be pissed
                remainder = site%3
                if alignment[1][(site - remainder):(site - remainder + 3)] in stops:
                    to_remove.extend(range((site - remainder), (site - remainder + 3)))
    if to_remove:
        #I'm doing this rather than deleting because if you deleted, you'd also have to change the positions and that's super annoying
        for pos in to_remove:
            alignment[0][pos] = "-"
            alignment[1][pos] = "-"
    return(alignment)

def remove_degenerate_differences(aligned_sequences, aligned_positions, other_aligned_positions):
    positions_to_remove = []
    for pos, position in enumerate(other_aligned_positions):
        if position in aligned_positions:
            first_base = aligned_sequences[0][position]
            try:
                other_base = aligned_sequences[1][position]
                if first_base != other_base:
                    current_codon = nc.index_to_codon(position)
                    positions_to_remove.extend(current_codon)
            except IndexError:
                if position > (len(aligned_sequences[0]) - 3 - 1):
                    pass
                else:
                    print("Something wonky! IndexError when you try to access positions {0} in the ortholog sequence!".format(position)) 
                    raise Exception
    aligned_positions = [i for i in aligned_positions if i not in positions_to_remove]
    return(aligned_positions)

def remove_stops(temp_human, temp_macaque):
    '''
    Remove in-frame stop codons from two sequences.
    '''
    to_delete = []
    stops = [["T","A","G"],["T","A","A"],["T","G","A"]]
    for pos in range(0,len(temp_human),3):
        if temp_human[pos:pos+3] in stops:
            to_delete.extend([pos,pos+1,pos+2])
        #when there is a stop in macaque but not in human, the codon should still be removed from both sequences
        elif temp_macaque[pos:pos+3] in stops:
            to_delete.extend([pos,pos+1,pos+2])
    temp_human = [temp_human[i] for i in range(len(temp_human)) if i not in to_delete]
    temp_macaque = [temp_macaque[i] for i in range(len(temp_macaque)) if i not in to_delete]
    return(temp_human,temp_macaque)

def run_codeml(cml_input, cml_output, ctl_file = None, method = None, id1 = None, id2 = None, model = None):
    '''
    Run PAML codeml/yn00 on an alignment.
    '''
    if not ctl_file:
        ctl_file = "temp_data/temp{0}.ctl".format(random.random())
    #if you want to use the Goldman and Yang method
    if method == "gy":
        #just a tree of two species
        #in theory, I could run it pairwise without a tree but as far as I understand, this has not been implemented in BioPython
        #but this is equivalent
        cml_tree = "general/human_macaque.tree"
        cmd_name = "codeml"
        cml = codeml.Codeml(alignment = cml_input, out_file = cml_output, tree = cml_tree)
        cml.set_options(verbose = 1, seqtype = 1, runmode = 0, model = 0, NSsites = [0])
        if ctl_file:
            cml.ctl_file = ctl_file
        try:
            temp_dict = cml.run(command = cmd_name)
        except OSError:
            print(cml_input)
            print(cml_output)
            print(ctl_file)
            raise Exception
        ds = temp_dict["NSsites"][0]["parameters"]["dS"]
        dn = temp_dict["NSsites"][0]["parameters"]["dN"]
        omega = temp_dict["NSsites"][0]["parameters"]["omega"]
        tree_length = None
    elif method == "yn":
        #if you want to use the Yang and Nielsen method (faster but less accurate)
        if not id1:
            id1 = "id"
        if not id2:
            id2 = "orth_id"
        cmd_name = "yn00"
        yn = yn00.Yn00(alignment = cml_input, out_file = cml_output)
        yn.ctl_file = ctl_file
        try:
            temp_dict = yn.run(command = cmd_name)
            ds = temp_dict[id2][id1]["YN00"]["dS"]
            dn = temp_dict[id2][id1]["YN00"]["dN"]
            omega = temp_dict[id2][id1]["YN00"]["omega"]
            tree_length = None
        except IndexError:
            #this is because the BioPython script that parses the yn00 output
            #expects floats but if there aren't enough synonymous sites
            #(at least that's the way I interpret it, not quite sure though)
            #dN and dS will be nans and you'll get an IndexError
            ds = None
            omega = None
            dn = None
            tree_length = None
    elif method == "baseml":
        bml_tree = "general/human_macaque.tree"
        cmd_name = "baseml"
        bml = baseml.Baseml(alignment = cml_input, out_file = cml_output, tree = bml_tree)
        if not model:
            model = 1
        bml.set_options(verbose = 1, runmode = 0, model = model)
        bml.ctl_file = ctl_file
        temp_dict = bml.run(command = cmd_name)
        tree_length = temp_dict["tree length"]
        ds = None
        dn = None
        omega = None
    else:
        print("{0} is not a valid method.".format(method))
    remove_file(ctl_file)
    return({"dS": ds, "dN": dn, "omega": omega, "tree length": tree_length})

def tabix(bed_file, output_file, genome, vcf = None):
    '''
    Given a bed file, use tabix to get overlapping 1000Genomes SNPs.
    '''
    #divide the input bed file into smaller files
    process_number = int(os.cpu_count()/2)
    bed_file_length = line_count(bed_file)
    #if the input bed_file has fewer lines than you were planning on using cores 
    if bed_file_length <= process_number:
        process_number = 2
    lines_per_file = int(bed_file_length/process_number)
    run_process(["split", "-l", lines_per_file, bed_file, bed_file])
    bed_names = ["{0}a{1}".format(bed_file, i) for i in string.ascii_lowercase]
    if (bed_file_length%process_number) == 0:
        bed_names = bed_names[:process_number]
    else:
        bed_names = bed_names[:(process_number + 1)]
    parallel_tabix = run_in_parallel(bed_names, ["foo", genome, vcf], tabix_core, workers = process_number)
    [i.get() for i in parallel_tabix]
    output_files = ["{0}.out".format(i) for i in bed_names]
    run_process(["cat {0}??.out".format(bed_file)], file_for_output = output_file, shell = True)
    [os.remove(i) for i in bed_names]
    [os.remove(i) for i in output_files]

def tabix_core(bed_files, genome, vcf):
    '''
    The code that's parallelized in tabix above.
    '''
    for curr_bed_file in bed_files:
        curr_output_file = curr_bed_file + ".out"
        with open(curr_bed_file) as file, open(curr_output_file, "w") as file2:
            counter = 0
            for line in file:
                if counter%100 == 0:
                    print(counter)
                counter = counter + 1
                line = line.split("\t")
                chrom = line[0].lstrip("chr")
                start = int(line[1]) + 1
                end = line[2]
                trans = line[3]
                #sometimes the process fails, I think it's just a matter of connecting to the FTP server.
                #so I'm making it try and try until it gets it.
                output = ("error")
                while output == ("error"):
                    if vcf:
                        output = run_process(["tabix", vcf, "{0}:{1}-{2}".format(chrom, start, end)])
                    else:
                        if genome == "hg37":
    ##                        output = run_process(["tabix", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz", "{0}:{1}-{2}".format(chrom, start, end)])
                            output = run_process(["tabix", "general/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.gz", "{0}:{1}-{2}".format(chrom, start, end)])

                        elif genome == "hg38":
    ##                        output = run_process(["tabix", "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz".format(chrom), "{0}:{1}-{2}".format(chrom, start, end)])
                            output = run_process(["tabix", "general/1000genomes/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz".format(chrom), "{0}:{1}-{2}".format(chrom, start, end)])
    ##                        remove_file("ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz.tbi".format(chrom))
                        else:
                            print("Invalid genome!")
                            print(genome)
                            raise Exception
                #in theory, this should never happen
                if output == ("error"):
                    print("You have failed!")
                    raise Exception
                if output:
                    output = output.rstrip("\n")
                    output = output.split("\n")
                    output = [i.split("\t") for i in output]
                    output = [i for i in output if i[2][:2] == "rs"]
                    if output:
                        output = ["\t".join(["chr{0}".format(i[0]), str(int(i[1]) - 1), i[1], trans, "100", ".", "$".join([i[2], i[3], i[4], i[6], i[7]])]) for i in output]
                        output = "\n".join(output)
                        file2.write(output)
                        file2.write("\n")

def tabix_samples(genome, bed_file, output_file_name, superpop = None, subpop = None, downsample_by = None):
    '''
    Extract 1000Genomes SNPs for a subpopulation.
    '''

    sex_chromosomes = ["Y", "X"]

    if not downsample_by:
        downsample_by = 1
        
    if genome != "hg38":
        print("Subpopulation filtering has only been implemented for hg38!")
        raise Exception

    panel = rw.read_many_fields("general/1000genomes/integrated_call_samples_v3.20130502.ALL.panel", "\t")
    panel = [i for i in panel if len(i) == 4]

    if subpop:
        samples = [i[0] for i in panel if i[1] == subpop]
    elif superpop:
        samples = [i[0] for i in panel if i[2] == superpop]
    else:
        samples = [i[0] for i in panel]

    samples = np.random.choice(samples, size = int(len(samples)/downsample_by), replace = False)

    print(len(samples))

    samples = ",".join(samples)

    with open(bed_file) as file:
        sample_files = []
        counter = 0
        for line in file:
            if counter%100 == 0:
                print(counter)
            counter = counter + 1
            line = line.split("\t")
            chrom = line[0].lstrip("chr")
            if chrom in sex_chromosomes:
                print("Only autosomes are processed!")
            else:
                start = int(line[1]) + 1
                end = line[2]
                trans = line[3]
                current_vcf = "general/1000genomes/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz".format(chrom)
                temp_output_file = "temp_data/temp_vcf{0}.vcf".format(random.random())
                run_process(["tabix", "-h", current_vcf, "{0}:{1}-{2}".format(chrom, start, end)], file_for_output = temp_output_file)
##                run_process(["cp", temp_output_file, "temp_data/{0}:{1}-{2}_tabix_slice.txt".format(chrom, start, end)])
                sample_output_file = "temp_data/temp_sample_tabix{0}.txt".format(random.random())
                sample_files.append(sample_output_file)
                run_process(["perl", "../Software/vcftools/src/perl/vcf-subset", "-c", samples, temp_output_file], file_for_output = sample_output_file)
                with open(sample_output_file) as tfile:
                    string = "".join(tfile)
                    counts = re.findall("AC=\d+", string)
                    counts = [float(i.split("=")[1]) for i in counts]
                    counts = [i for i in counts if i > 216]
                    if counts:
                        print(counts)
                        print(sample_output_file)
                        print(temp_output_file)
                        raise Exception
##                        run_process(["cp", sample_output_file, "temp_data/{0}:{1}-{2}_tabix_subset.txt".format(chrom, start, end)])
                remove_file(temp_output_file)

    concat_files = ["temp_data/temp_concat_file{0}.vcf".format(random.random()), "temp_data/temp_concat_file{0}.vcf".format(random.random())]
    current_sample_files = sample_files[-10:]
    del sample_files[-10:]
    run_process(["perl", "../Software/vcftools/src/perl/vcf-concat"] + current_sample_files, file_for_output = concat_files[0])
    local_counter = 0
    while len(sample_files) > 0:
        local_counter = local_counter + 1
        current_sample_files = sample_files[-10:]
        del sample_files[-10:]
        if local_counter%2 == 0:
            current_concat_file = concat_files[0]
            previous_concat_file = concat_files[1]
        else:
            current_concat_file = concat_files[1]
            previous_concat_file = concat_files[0]
        run_process(["perl", "../Software/vcftools/src/perl/vcf-concat"] + current_sample_files + [previous_concat_file], file_for_output = current_concat_file)
    sort_file = "temp_data/temp_sort_file{0}.vcf".format(random.random())
    run_process(["vcf-sort", current_concat_file], file_for_output = sort_file)
    run_process(["bgzip", "-c", sort_file], file_for_output = output_file_name)
    run_process(["tabix", "-f", "-p", "vcf", output_file_name])
    for sample_file in sample_files:
        remove_file(sample_file)
    for concat_file in concat_files:
        remove_file(concat_file)
    remove_file(sort_file)
        
def weight_cons_by_dinucl(freqs_dict, dinucl):
    '''
    Calculate the rate of evolution in motifs vs non-motifs, controlling for any differences in dinucleotide frequencies.
    '''
    motif_count = 0
    motif_sum = 0
    non_motif_count = 0
    non_motif_sum = 0
    #loop over the dinuncleotides
    for dint in dinucl:
        #freqs_dict comes from cons_by_dinucl above
        if (freqs_dict[dint]["subst. in motifs"] != None) and (freqs_dict[dint]["subst. in non-motifs"] != None):
            #weight the rate of change by the frequency of the dinucleotide
            current_motif_rate = freqs_dict[dint]["subst. in motifs"] * freqs_dict[dint]["frequency in motifs"]
            motif_sum = motif_sum + current_motif_rate
            motif_count = motif_count + 1
            #this is not an error: you need to be taking the rate from the non-motifs but the frequency from the motifs!
            current_non_motif_rate = freqs_dict[dint]["subst. in non-motifs"] * freqs_dict[dint]["frequency in motifs"]
            non_motif_sum = non_motif_sum + current_non_motif_rate
            non_motif_count = non_motif_count + 1
    #get the average weighted rate
    result = {}
    if motif_count > 0:
        result["motif rate"] = motif_sum/motif_count
    else:
        result["motif rate"] = None
    if non_motif_count > 0:
        result["non-motif rate"] = non_motif_sum/non_motif_count
    else:
        result["non-motif rate"] = None
    if (motif_count > 0) and non_motif_count > 0:
        result["norm. rate"] = (result["motif rate"] - result["non-motif rate"])/result["non-motif rate"]
    else:
        result["norm. rate"] = None
    return(result)

def write_hits_to_phylip(fasta, positions, output_file, correspondances, alignments, degen_file, baseml = False, regions = False, fs = None, global_fasta = None):
    '''
    Given a list of positions for a series of ORFs, extract the relevant sequences from a pairwise alignment, concatenate across sequences and write to phylip format.
    '''
    if degen_file:
        degen = parse_degen(degen_file)
    correspondances = rw.read_many_fields(correspondances, ",")
    correspondances = list_to_dict(correspondances, 0, 1)
    names, seqs = rw.read_fasta(fasta)
    region_seqs = seqs.copy()
    if regions:
        names_reg = names.copy()
        names, seqs = rw.read_fasta(global_fasta)
        bed_name = fasta.replace("fasta", "bed")
        regions_bed = rw.read_many_fields(bed_name, "\t")
        transcripts = fs.get_transcripts()
        CDSs = fs.get_CDS()
        mapping_dict = map_regions_to_CDS(fasta, bed_name, fs, transcripts, CDSs, trans_ids = True)
    to_concat = [[""], [""]]

    for pos, name in enumerate(sorted(positions.keys())):
        current_pos = positions[name].copy()
        if regions:
            trans_id = regions_bed[names_reg.index(name)][3]
            flank_indices = mapping_dict[name]["flank indices in CDS"]
            name = mapping_dict[name]["idn"]
            current_pos = region_indices_to_full_indices(current_pos, flank_indices)
        seq = seqs[names.index(name)]
        if not baseml:
            twofold = get_LeuArg(list(seq), 2)
            fourfold = get_LeuArg(list(seq), 4)
            current_pos, to_change_2, to_change_4 = fill_codons(current_pos, twofold, fourfold)
        try:
            alignment = [list(str(i.seq)) for i in SeqIO.parse("{0}/{1}_{2}.phy".format(alignments, name, correspondances[name]), "phylip-sequential")]
        except FileNotFoundError:
            alignment = [list(str(i.seq)) for i in SeqIO.parse("{0}/{1}_{2}.phy".format(alignments, name, "{0}_0".format(correspondances[name])), "phylip-sequential")]            
        current_pos = get_aligned_positions(alignment, current_pos)
        if len(current_pos) % 3 != 0 and not baseml:
            print("Problem filling up codons!")
            print(name)
            print(current_pos)
            print(len(current_pos))
            raise Exception
        if not baseml:
            to_change_2 = get_aligned_positions(alignment, to_change_2)
            to_change_4 = get_aligned_positions(alignment, to_change_4)
            alignment = changeLeuArg(to_change_2, to_change_4, alignment)
        if degen_file:
            current_degen = degen[name]
            if len(current_degen) > 0:
                current_degen_keys = sorted(list(current_degen.keys()))
                current_degen_keys_aligned = get_aligned_positions(alignment, current_degen_keys)
                current_degen = {current_degen_keys_aligned[current_degen_keys.index(i)]: current_degen[i] for i in current_degen}
                alignment = reencode_disruption_pairwise(alignment, [i for i in current_pos if i in current_degen], current_degen)
        to_concat[0].append("".join([alignment[0][i] for i in current_pos]))
        to_concat[1].append("".join([alignment[1][i] for i in current_pos]))
    to_concat = ["".join(i) for i in to_concat]
    if len(to_concat[0]) > 0:
        alignment = MultipleSeqAlignment([SeqRecord(Seq(to_concat[0], IUPAC.unambiguous_dna), id = "seq1"), SeqRecord(Seq(to_concat[1], IUPAC.unambiguous_dna), id = "seq2")])
        AlignIO.write(alignment, output_file,"phylip-sequential")
    else:
        raise NoDataException
