'''
Author: Rosina Savisaar.
Module that contains generic utility functions that make life a bit easier.
'''

import argparse
import itertools as it
import multiprocessing as mp
import os
import re
import subprocess
import sys
import time

def check_equal_lengths(several_lists):
    '''
    Check that two or more lists all have the same number of elements.
    Only works with flat lists.
    '''
    lengths = [len(i) for i in several_lists]
    if len(list(set(lengths))) != 1:
        print("Problem extracting data!")
        raise Exception

def clean_list(to_delete, list_to_clean):
    '''
    Given a list and a second list with indices, delete the elements in the first list
    specified by the indices in the second list.
    '''
    to_delete = sorted(list(set(to_delete)))
    if len(to_delete) > 0:
        for i in range(0,len(to_delete)):
            del list_to_clean[to_delete[i]-i]
    return(list_to_clean)

def fourfold_deg_list():
    '''
    Get a list of all the fourfold degenerate codons.
    '''
    codons = ["TGA","TAG","TAA","TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"]
    aas =  ["STOP","STOP","STOP","Phe", "Phe","Leu", "Leu", "Leu", "Leu", "Leu", "Leu", "Ile", "Ile", "Ile", "Met", "Val", "Val", "Val", "Val", "Ser", "Ser", "Ser", "Ser", "Pro", "Pro", "Pro", "Pro", "Thr", "Thr", "Thr", "Thr", "Ala", "Ala", "Ala", "Ala", "Tyr", "Tyr", "His", "His", "Gln", "Gln", "Asn", "Asn", "Lys", "Lys", "Asp", "Asp", "Glu", "Glu", "Cys", "Cys", "Trp", "Arg", "Arg", "Arg", "Arg", "Ser", "Ser", "Arg", "Arg", "Gly", "Gly", "Gly", "Gly"]
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

def flatten(structured_list):
    '''
    Flatten a structured list.
    '''
    flat_list = list(it.chain(*structured_list))
    return(flat_list)

def get_extension(file_name, extension_length, valid_list = None):
    '''
    Determine the extension at the end of a file name.
    '''
    extension = file_name[-extension_length:]
    if valid_list:
        if extension not in valid_list:
            print("File format must be included in {0}!".format(valid_list))
            sys.exit()
    return(extension)

def line_count(file):
    '''
    Count the number of lines in a file.
    '''
    #not using wc -l because I want the number of lines, not the number of newlines.
    output = run_process(["grep", "-c", "^", file])
    return(int(output))

def list_to_dict(input_list, index1, index2, as_list = False, uniquify = False, floatify = False, intify_key = False):
    '''
    Convert the input_list into a dictionary, with the index1th element of each sublist as the key and the index2th element as the value.
    '''
    if as_list and floatify:
        print("_as_list_ and _floatify_ can't both be True!")
        raise Exception
    output_dict = {}
    for i in input_list:
        if not as_list:
            if floatify:
                #convert the value into a float
                output_dict[i[index1]] = float(i[index2])
            else:
                output_dict[i[index1]] = i[index2]
        else:
            #if several sublists can have the same value in index1 and you want all their index2th values as a list
            if i[index1] not in output_dict:
                output_dict[i[index1]] = []
            output_dict[i[index1]].append(i[index2])
    if as_list and uniquify:
        #if the values are lists and you don't want duplicates in the value lists
        output_dict = {i: sorted(list(set(output_dict[i]))) for i in output_dict}
    if intify_key:
        output_dict = {int(i): output_dict[i] for i in output_dict}
    return(output_dict)

def make_dir(dir_name):
    '''
    Check whether a directory exists and if not, create it.
    '''
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)   

def order_sequences(ordered_identifiers,identifiers,sequences):
    '''
    Given a list of identifiers and a list of sequences where the order matches between the two lists,
    take a third list of ordered identifiers and order the sequences according to that list instead.
    '''
    ordered_sequences = [0 for i in sequences]
    for j in range(len(identifiers)):
        ordered_sequences[j] = sequences[identifiers.index(ordered_identifiers[j])]
    return(ordered_sequences)

def overlap(a, b):
    '''
    Given two lists, determine which elements appear in both.
    '''
    set_a = set(a)
    intersection = set_a.intersection(b)
    return(list(intersection))

def parse_arguments(description, arguments, floats = None, flags = None, ints = None):
    '''
    Use argparse to parse a set of input arguments from the command line.
    '''
    if not floats:
        floats = []
    if not flags:
        flags = []
    if not ints:
        ints = []
    parser = argparse.ArgumentParser(description = description)
    for pos, argument in enumerate(arguments):
        if pos in flags:
            parser.add_argument("--{0}".format(argument), action = "store_true", help = argument)
        else:
            if pos in floats:
                curr_type = float
            elif pos in ints:
                curr_type = int
            else:
                curr_type = str
            parser.add_argument(argument, type = curr_type, help = argument)
    args = parser.parse_args()
    return(args)

def pause_script():
    '''
    Pause script indefinitely.
    '''
    while True:
        time.sleep(100)
    
def print_elements(input_list):
    '''
    Take a list and print out the elements separated by carriage returns.
    '''
    for i in input_list:
        print(i)
    print("\n")

def remove_file(file_name):
    '''
    Remove a file, if it exists.
    '''
    try:
        os.remove(file_name)
    except FileNotFoundError:
        pass

def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):
    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    input_list: a list of the stuff you want to parallelize over (for example a list of gene names)
    args: a list of arguments to the function
    func: the function
    kwargs_dict: a dictionary of any keyword arguments the function might take
    workers: number of parallel processes to launch
    onebyone: if True, allocate one element from input_list to each process
    '''
    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    #in the list of arguments, I put in "foo" for the argument that corresponds to whatever is in the input_list because I couldn't be bothered to do something less stupid
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        #divide input_list into as many chunks as you're going to have processes
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = mp.Pool(workers)
    results = []
    #go over the chunks you made and laucnh a process for each
    for i in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = i
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))            
        results.append(process)
    pool.close()
    pool.join()
    return(results)

def run_process(arguments, return_string = True, input_to_pipe = None, return_error = False, file_for_input = None, file_for_output = None, univ_nl = True, shell = False):
    '''
    Run a command on the command line.
    '''
    if file_for_input:
        input_file = open(file_for_input)
        stdin_src = input_file
    else:
        stdin_src = subprocess.PIPE
    if file_for_output:
        output_file = open(file_for_output, "w")
        stdout_dest = output_file
    else:
        stdout_dest = subprocess.PIPE
    arguments = [str(i) for i in arguments]
    if shell:
        arguments = " ".join(arguments)
    process = subprocess.Popen(arguments, shell = shell, stdout = stdout_dest, stderr = subprocess.PIPE,
                               stdin = stdin_src, universal_newlines = univ_nl)
    if input_to_pipe:
        stdout, stderr = process.communicate(input_to_pipe)
    else:
        stdout, stderr = process.communicate()
    if file_for_input:
        input_file.close()
    if file_for_output:
        output_file.close()
    return_code = process.poll()
    if return_code != 0:
        print("Process failed!")
        print(" ".join(arguments))
        print(stderr)
        return("error")
    #if the process returns bytes but you want to get a string back.
    if return_string and type(stdout) == bytes:
        stdout = stdout.decode("utf-8")
    if return_error:
        return(stderr)
    else:
        return(stdout)

def update_counter(counter, step):
    '''
    Print out and update counter.
    '''
    if counter % step == 0:
        print(counter)
    counter = counter + 1
    return(counter)
