'''
Author: Rosina Savisaar.
Directly compare the frequency of segregating sites/mean allele frequency between hits and controls.
'''

from housekeeping import flatten, parse_arguments, print_elements, remove_file, run_process
from INSIGHT import shuffle_dictionaries
import numpy as np
import random
import re
import read_and_write as rw
import scipy.stats

def get_chisq_site_freq(hit_data, control_data):
    '''
    Use a chi-squared test to compare the proportion of polymorphic sites
    among hits and controls.
    '''
    total_hits = len(hit_data)
    total_controls = len(control_data)
    poly_hits = len([i for i in hit_data if i[2] != "M"])
    poly_controls = len([i for i in control_data if i[2] != "M"])
    exp_poly = poly_controls/total_controls * total_hits
    mono_hits = total_hits - poly_hits
    exp_mono = total_hits - exp_poly
    chi_sq = scipy.stats.chisquare([poly_hits, mono_hits], [exp_poly, exp_mono])
    print("Hit sites: {0}, of which {1} polymorphic.".format(total_hits, poly_hits))
    print("Control sites: {0}, of which {1} polymorphic.".format(total_controls, poly_controls))
    print("Chisq: {0}.".format(chi_sq[0]))
    print("p: {0}.".format(chi_sq[1]))
    print("Difference: {0}.\n".format((poly_hits/total_hits) - (poly_controls/total_controls)))
    return((poly_hits/total_hits) - (poly_controls/total_controls))

def get_data(file):
    '''
    Read in polymorphism data from an INSIGHT input file.
    '''
    data = rw.read_many_fields(file, "\t")
    data = [i for i in data if i[0] == "site"]
    return(data)

def get_mean_freq(SFS_file):
    '''
    Use a Mann-Whitney U-test to compare MAFs in hits vs controls.
    '''
    SFS = rw.read_many_fields(SFS_file, " ")
    n = int(SFS[0][0])
    hit_freqs = flatten([[i/n for j in range(int(SFS[1][i]))] for i in range(1, len(SFS[1]))])
    control_freqs = flatten([[i/n for j in range(int(SFS[2][i]))] for i in range(1, len(SFS[2]))])
    mwu = scipy.stats.mannwhitneyu(control_freqs, hit_freqs)
    hit_median = np.median(hit_freqs)
    control_median = np.median(control_freqs)
    print("Median MAF in hits: {0}.".format(hit_median))
    print("Median MAF in controls: {0}.".format(control_median))
    print("MWU p: {0}.".format(mwu[1]))
    print("Difference: {0}.\n".format(hit_median - control_median))
    return([hit_freqs, control_freqs], hit_median - control_median)

def main():
    description = "Directly compare the frequency of segregating sites/mean allele frequency between hits and controls."
    args = parse_arguments(description, ["hit_file", "control_file", "INSIGHT_hit_file", "INSIGHT_control_file", "SFS_file", "trial_file", "trials", "shuffle"], ints = [6], flags = [7])
    hit_file, control_file, INSIGHT_hit_file, INSIGHT_control_file, SFS_file, trial_file, trials, shuffle = args.hit_file, args.control_file, args.INSIGHT_hit_file, args.INSIGHT_control_file, args.SFS_file, args.trial_file, args.trials, args.shuffle

    true_hits = rw.read_pos(hit_file)
    true_controls = rw.read_pos(control_file)

    #to store the original data in case this is a negative control and you will be shuffling
    #hits and controls
    original_INSIGHT_hit_file = INSIGHT_hit_file
    original_INSIGHT_control_file = INSIGHT_control_file

    print(hit_file)

    with open(trial_file, "w") as file:
        file.write("trial\tpoly_fraction_hits - poly_fraction_controls\tmedian_hit_MAF - median_control_MAF\n")
        for trial in range(trials):
            to_write = "{0}\t".format(trial)

            #if this is a negative control
            if shuffle:
                INSIGHT_hit_file = re.sub("_0_", "_{0}_".format(trial), original_INSIGHT_hit_file)
                INSIGHT_control_file = re.sub("_0_", "_{0}_".format(trial), original_INSIGHT_control_file)
                temp_hits_file = "temp_data/temp_hits{0}.txt".format(random.random())
                temp_controls_file = "temp_data/temp_controls{0}.txt".format(random.random())
                #shuffle hits and controls
                temp_hits, temp_controls = shuffle_dictionaries(true_hits, true_controls)
                rw.write_pos(temp_hits, temp_hits_file)
                rw.write_pos(temp_controls, temp_controls_file)
                SFS_file = "temp_data/temp_SFS_file{0}.txt".format(random.random())
                #generate an ISNIGHT input file that you could then use for the manual analysis
                run_process(["python3", "mDFEest_input.py", temp_hits_file, temp_controls_file,
                             "general/1000genomes/filtered_hg38_85_pc_multiexon_Yoruban_SNPs_relative.txt", 216,
                             SFS_file])
                remove_file(temp_hits_file)
                remove_file(temp_controls_file)
        
            hit_data = get_data(INSIGHT_hit_file)
            control_data = get_data(INSIGHT_control_file)

            poly_ratio_diff = get_chisq_site_freq(hit_data, control_data)
            to_write = to_write + "{0}\t".format(poly_ratio_diff)

            temp, median_diff = get_mean_freq(SFS_file)
            to_write = to_write + "{0}\n".format(median_diff)

            if shuffle:
                remove_file(SFS_file)
                
            file.write(to_write)
       
if __name__ == "__main__":
    main()
