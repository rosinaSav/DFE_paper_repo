'''
Author: Rosina Savisaar.
Modue that contains custom statistical methods.
'''

from housekeeping import run_process
import numpy as np
import random
import re
import scipy.misc
import sys

def calc_eff_p(real_value, sim_values, greater = True):
    '''
    Given an estimate and a series of simulated estimates, calculate and empirical effective p-value.
    If greater is True, calculate the porbbaility that a value this great or greater would have been observed by chance,
    otherwise that a value this low or lower would have been observed.
    '''
    if real_value == None or np.isnan(real_value):
        return(None)
    sim_values = [i for i in sim_values if i != None and not np.isnan(i)]
    if greater:
        more_extreme = [i for i in sim_values if i >= real_value]
    else:
        more_extreme = [i for i in sim_values if i <= real_value]
    n = len(more_extreme)
    m = len(sim_values)
    p = (n + 1)/(m + 1)
    return(p)

def calc_Zscore(real_value, sim_values):
    '''
    Convert an estimate into a Z-score based on a set of simulated values.
    '''
    if real_value == None or np.isnan(real_value):
        return(None)
    sim_values = [i for i in sim_values if i != None and not np.isnan(i)]
    average = np.mean(sim_values)
    stdev = np.std(sim_values, ddof=1)
    if stdev == 0:
        return(None)
    Zscore = (real_value - average)/stdev
    return(Zscore)

def chi_test(observed, expected):
    '''
    Given a series of observed and expected values, conduct a chi-squared test.
    '''
    observed = ",".join([str(i) for i in observed])
    expected = ",".join([str(i) for i in expected])
    string_to_R = "|".join([observed, expected])
    output = run_process(["Rscript", "R_scripts/chi_test.r", string_to_R])
    output = output.split(" ")
    output = [i for i in output if i != ""]
    chi = float((output[1].lstrip("\"")).rstrip("\""))
    p = float(((output[4].lstrip("\"")).rstrip("\n")).rstrip("\""))
    return({"chi": chi, "p": p})

def correct_multiple_testing(p_values, method):
    '''
    Given a list of p-values, correct them for multiple testing.
    '''
    p_values = [str(i) for i in p_values]
    p_values.append(method)
    string_to_R = ",".join(p_values)
    corrected_values = run_process(["Rscript", "R_scripts/holm_correct.r", string_to_R])
    corrected_values = re.findall("[\d\.]*", corrected_values, re.MULTILINE)
    corrected_values = [float(i) for i in corrected_values if "." in i]
    if (len(p_values)) - 1 != len(corrected_values):
        print("Problem correcting for multiple comparisons!")
        print(p_values)
        print(corrected_values)
        sys.exit()
    return(corrected_values)

def fishers_exact_test(observed, expected):
    '''
    Perform a Fisher's exact test on an observed and an expected proportion
    '''
    string_to_R = ",".join([str(observed[0]), str(observed[1]), str(expected[0]), str(expected[1]), "greater"])
    results = run_process(["Rscript", "R_scripts/fisher_test.r", string_to_R])
    results = results.rstrip("\"\n")
    #sometimes there's a space between the quotation marks and the newline
    results = results.rstrip("\" \n")
    results = results.split(" ")
    results = [(i.rstrip("\"")).lstrip("\"") for i in results if i != ""]
    results = results[1:]
    results = [float(i) if "Inf" not in i else i for i in results]
    return(results)

def fishers_exact_test_enrichment(element, sample, population, alt):
    '''
    Perform a Fisher's exact test to check whether a given element is enriched in a sample when compared to a population.
    '''
    N = len(population)
    n = len(sample)
    if len(sample) >= len(population):
        print("The sample has to be smaller than the population!")
        raise Exception
    K = population.count(element)
    k = sample.count(element)
    string_to_R = ",".join([str(k), str(n - k), str(K), str(N - K), alt])
    results = run_process(["Rscript", "R_scripts/fisher_test.r", string_to_R])
    results = results.rstrip("\"\n")
    #sometimes there's a space between the quotation marks and the newline
    results = results.rstrip("\" \n")
    results = results.split(" ")
    results = [(i.rstrip("\"")).lstrip("\"") for i in results if i != ""]
    results = results[1:]
    results = [float(i) if "Inf" not in i else i for i in results]
    return(results)

def normalize(real_value, simulated_values, use_median = False):
    '''
    Normalize an estimate with regards to a list of simulated estimates (normalized value = (real value - simulated mean)/simulated mean).
    '''
    if real_value == None or np.isnan(real_value):
        return(None)
    simulated_values = np.array([i for i in simulated_values if i != None and not np.isnan(i)])
    if np.sum(simulated_values) == 0:
        return(None)
    if not use_median:
        simulated_average = np.mean(simulated_values)
    else:
        simulated_average = np.median(simulated_values)
    normalized_value = np.divide((real_value - simulated_average), simulated_average)
    return(normalized_value)

def wilcoxon_signed_rank_test(vector1, vector2, alt):
    '''
    Perform a Mann-Whitney U test to compare two samples. The alternative must be one of
    "greater", "less" and "two.tailed".
    '''
    vector1 = ",".join([str(i) for i in vector1])
    vector2 = ",".join([str(i) for i in vector2])
    vectors = "|".join([vector1, vector2])
    string_to_R = "_".join([vectors, alt])
    results = run_process(["Rscript", "R_scripts/wilcoxon_signed_rank_test.r", string_to_R])
    results = results.rstrip("\n")
    results = results.split(" ")
    results = [i for i in results if i != ""]
    results = (results[1].rstrip("\"")).lstrip("\"")
    results = float(results)
    return(results)

