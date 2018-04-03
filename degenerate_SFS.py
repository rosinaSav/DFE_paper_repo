'''
Author: Rosina Savisaar.
Construct a site frequency spectrum that only considers motif-disrupting SNPs.
'''

from housekeeping import list_to_dict, parse_arguments, update_counter
from mDFEest_input import get_SFS
import nucleotide_comp as nc
import random
import re
import read_and_write as rw

def check_disruption(motif_regex, current_seq, motifs, motif_lengths, fourfold_pos, full_SNPs, clean_SNPs, minor_alleles, trans, transitions_total, transitions_disr, hit_degen_file, to_remove):
    '''
    Check which SNPs are motif disrupting.
    '''
    all_sites = []
    base_dict = {}
    #write down transcript name
    hit_degen_file.write("{0}\t".format(trans))
    #for each motif
    for pos, motif_reg in enumerate(motif_regex):
        #find positions that overlap the motif
        positions = re.finditer(motif_reg, current_seq)
        #for each match
        for obj in positions:
            span = obj.span()
            #get hit positions
            current_pos = [i for i in range(span[0], span[0] + motif_lengths[pos])]
            #restrict to 4fold degenerate sites
            clean_pos = [i for i in current_pos if i in fourfold_pos]
            #if any hit positions remain
            if clean_pos:
                current_motif = motifs[pos]
                #for all hit positions
                for site in clean_pos:
                    #if the site overlaps with a SNPs
                    if (site in full_SNPs[trans]) or (trans not in to_remove) or (site not in to_remove[trans]):
                        if site not in base_dict:
                            base_dict[site] = []
                        all_sites.append(site)
                        #get the index of the hit position in the motif
                        index = current_pos.index(site)
                        #the base that is actually at that position
                        ref_allele = current_motif[index]
                        temp_motif = current_motif.copy()
                        #try substituting in each of the bases and check whether that turns the motif into a non-motif
                        for base in nc._canon_bases_:
                            if base != ref_allele:
                                transitions_total[ref_allele][base] = transitions_total[ref_allele][base] + 1
                                temp_motif[index] = base
                                if temp_motif not in motifs:
                                    transitions_disr[ref_allele][base] = transitions_disr[ref_allele][base] + 1
                                    base_dict[site].append(base)
                        if site in full_SNPs[trans]:
                            temp_motif[index] = minor_alleles[trans][site]
                            if temp_motif not in motifs:
                                #note that some sites might appear several times but because the SFS is created
                                #iterating over hit positions, this doesn't matter
                                clean_SNPs[trans][site] = full_SNPs[trans][site]
    all_sites = list(set(all_sites))
    for position in sorted(base_dict.keys()):
        hit_degen_file.write("{0}:{1},".format(position, "|".join(sorted(list(set(base_dict[position]))))))
    return(all_sites, clean_SNPs, transitions_total, transitions_disr, hit_degen_file)

def get_disrupt_bases(ref_allele, transitions):
    '''
    For each SNP, assign it as disruptive with the same probability as the frequency
    of disruptive substitutions for all possible substitutions having the same reference and substitute allele.
    Note that the latter frequencies are not based on actual substitutions/SNPs,
    they are based on all potential substitutions (so any A site might potentially substitute to T, C or G).
    '''
    disrupt_bases = []
    for base in nc._canon_bases_:
        if base != ref_allele:
            #frequency at which this SNP would be disruptive in real data
            trans_prob = transitions[ref_allele][base]
            #pick random bumber between 0 and 1. If it is below the probability,
            #assign the SNP as disruptive
            picked = random.random()
            if picked < trans_prob:
                disrupt_bases.append(base)
    return(disrupt_bases)

def get_transitions(transitions_disr, transitions_total):
    '''
    Based on a dictionary of the form [base1][base2][count of substitutions] that stores how many of each
    type of potential substitution would be motif disruptive and a similar dictionary that stores counts of all potential
    substitutions, calculate the frequency at which each type of potential subsitution is disruptive.
    '''
    transitions = {}
    for base in transitions_disr:
        transitions[base] = {}
        for base2 in transitions_disr[base]:
            #if no such substitutions were observed,
            #the probability of the substitution
            #being disruptive is set at 0
            if transitions_total[base][base2] == 0:
                transitions[base][base2] = 0
            else:
                transitions[base][base2] = transitions_disr[base][base2]/transitions_total[base][base2]
    return(transitions)

def parse_SNPs(trans_SNPs, clean_SNPs, full_SNPs, minor_alleles, trans):
    '''
    Parse an in-house SNP file into convenient dictionaries.
    '''
    full_SNPs[trans] = {}
    #for clean_SNPs you just add in the key at this stage and initialize
    #a dictionary. The contents of that dictionary will be filled in later
    #in check_disruption(...).
    clean_SNPs[trans] = {}
    minor_alleles[trans] = {}
    if trans_SNPs:
        trans_SNPs = [i.split(",") for i in trans_SNPs.split("|")]
        #this is where you get the minor allele
        minor_alleles[trans] = list_to_dict(trans_SNPs, 0, 1)
        minor_alleles[trans] = {int(i): minor_alleles[trans][i] for i in minor_alleles[trans]}
        #this is where you get the allele count
        trans_SNPs = list_to_dict(trans_SNPs, 0, 3)
        trans_SNPs = {int(i): int(trans_SNPs[i]) for i in trans_SNPs}
        full_SNPs[trans] = trans_SNPs
    return(trans_SNPs, clean_SNPs, full_SNPs, minor_alleles)       

def main():

    description = "Construct a site frequency spectrum that only considers motif-disrupting SNPs."
    args = parse_arguments(description, ["fasta", "output_file", "motif_file", "anc_file", "control_file", "SNPs_file", "N", "old_motif_format", "human", "ancestral"], ints = [6], flags = [7, 8, 9])
    fasta, output_file, motif_file, anc_file, control_file, SNPs_file, N, old_motif_format, human, ancestral = args.fasta, args.output_file, args.motif_file, args.anc_file, args.control_file, args.SNPs_file, args.N, args.old_motif_format, args.human, args.ancestral

    names, seqs = rw.read_fasta(fasta)

    #I use two different formats for storing sequence motifs,
    #got to know which on it is
    if old_motif_format:
        motifs = rw.read_names(motif_file)[1:]
        print(len(motifs))
    else:
        motifs = rw.read_motifs(motif_file)
        motifs = sorted(list(set(flatten(list(motifs.values())))))

    #get the lengths of the motifs and compile lookahead regexes
    #that recognize the whole motif but only store the position of the first bases
    #these will be needed when searchin for the motifs
    motif_lengths = [len(i) for i in motifs]
    motif_regex = nc.motif_to_regex(motifs)

    #I'm gonna treat CG and GC as two 2-bp motifs, use the same code as wehn searching for, say,
    #ESE motifs
    CG_2mers = ["CG", "GC"]
    CG_lengths = [2, 2]
    CG_regex = nc.motif_to_regex(CG_2mers)

    motifs = [list(i) for i in motifs]

    if ancestral:
        anc_pos = rw.read_pos(anc_file)

    #read in hit and control positions
    controls = rw.read_pos(control_file)
    hit_file = re.sub("controls", "hits", control_file)
    hits = rw.read_pos(hit_file)

    #read in SNP data
    SNPs = rw.read_many_fields(SNPs_file, "\t")
    #the second column in the SNPs file contains positions that need to be discarded from analysis because they contain unanalyzable SNP data
    to_remove = list_to_dict(SNPs, 0, 2)
    to_remove = {i: to_remove[i].split(",") for i in to_remove}
    to_remove = {i: [int(j) for j in to_remove[i] if j not in ["error", ""]] for i in to_remove}
    SNPs = list_to_dict(SNPs, 0, 1)

    #all the SNPs associated to a transcript
    full_SNPs = {}
    #disruptive SNPs only
    clean_SNPs = {}
    minor_alleles = {}

    #the number of hit positions where, say, a T could theoretically substitute to an A (i.e. all T positions)
    transitions_total = {i: {j: 0 for j in nc._canon_bases_} for i in nc._canon_bases_}
    #the same as above but only counting those substitutions that would turn a motif into a non-motif
    transitions_disr = {i: {j: 0 for j in nc._canon_bases_} for i in nc._canon_bases_}

    #this block of code filters the true SNPs to only leave those that are disruptive
    #and also calculates the probability of being disruptive for all potential SNPs
    with open("{0}_degen.txt".format(hit_file), "w") as hit_degen_file:
        counter = 0
        for trans in names:
            counter = update_counter(counter, 1000)
            if trans in controls:
                if trans in SNPs:
                    trans_SNPs = SNPs[trans]
                else:
                    trans_SNPs = []
                trans_SNPs, clean_SNPs, full_SNPs, minor_alleles = parse_SNPs(trans_SNPs, clean_SNPs, full_SNPs, minor_alleles, trans)
                current_seq = seqs[names.index(trans)]
                fourfold_pos = nc.get_4fold_deg(current_seq)
                #CpG filtering
                if human:
                    CG_pos = nc.get_motif_set_density(CG_regex, CG_lengths, current_seq, concat = True)["positions"]
                    fourfold_pos = [i for i in fourfold_pos if i not in CG_pos]
                if ancestral:
                    fourfold_pos = [i for i in fourfold_pos if i not in anc_pos[trans]]
                all_sites, clean_SNPs, transitions_total, transitions_disr, hit_degen_file = check_disruption(motif_regex, current_seq, motifs, motif_lengths, fourfold_pos, full_SNPs, clean_SNPs, minor_alleles, trans, transitions_total, transitions_disr, hit_degen_file, to_remove)
                hit_degen_file.write("\n")

    to_remove = {i: [j for j in to_remove[i] if j not in full_SNPs[i]] for i in to_remove if i in controls}

    hit_SFS = get_SFS(hits, clean_SNPs, to_remove, N)

    transitions = get_transitions(transitions_disr, transitions_total)
    print(transitions)

    #this block randomly assigns certain SNPs at simulant positions to be disruptive,
    #with the probability of that happening proportional to the frequency with which potential substitutions
    #of that nucleotide composition would be disruptive for true (motif) sites
    with open("{0}_degen.txt".format(control_file), "w") as control_degen_file:
        control_SNPs = {}
        counter = 0
        for trans in controls:
            control_degen_file.write("{0}\t".format(trans))
            counter = update_counter(counter, 1000)
            control_SNPs[trans] = {}
            trans_SNPs = full_SNPs[trans]
            current_seq = seqs[names.index(trans)]
            for site in controls[trans]:
                if trans not in to_remove or site not in to_remove[trans]:
                    ref_allele = current_seq[site]
                    disrupt_bases = get_disrupt_bases(ref_allele, transitions)
                    control_degen_file.write("{0}:{1},".format(site, "|".join(disrupt_bases)))
                    if site in trans_SNPs:
                        minor_allele = minor_alleles[trans][site]
                        if minor_allele in disrupt_bases:
                            control_SNPs[trans][site] = trans_SNPs[site]
            control_degen_file.write("\n")

    control_SFS = get_SFS(controls, control_SNPs, to_remove, N)

    with open(output_file, "w") as file:
        file.write("{0}\n".format(N))
        file.write(" ".join([str(i) for i in hit_SFS]))
        file.write("\n")
        file.write(" ".join([str(i) for i in control_SFS]))
        file.write("\n")    
                                                                                                                           
if __name__ == "__main__":
    main()
