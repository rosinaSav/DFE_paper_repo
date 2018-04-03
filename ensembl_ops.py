'''
Author: Rosina Savisaar.
Functions for interacting with the Ensembl database and for processing data
retrieved from it.
'''

from bedtools_games import sort_coords
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
##import _mysql
from housekeeping import flatten, line_count, print_elements, remove_file, run_process
import random
import re

def concat_sequences(current_dict, min_species, global_name):
    '''
    Given a set of species dictionaries from filter_species, prepare the string that is to be written to the output file.
    ''' 
    if "homo_sapiens" in current_dict:
        strand = global_name[6]
        phase = global_name[-1]
        #to make it easier to loop over the sequences
        #really stupid
        seqs = {i: current_dict[i]["seqs"] for i in current_dict}
        #join and, if necessary, reorder and reverse complement the sequences
        seqs = {i: [j.upper() for j in seqs[i]] for i in seqs}
        if strand == "-":
            seqs = {i: [str(Seq(j, IUPAC.unambiguous_dna).reverse_complement()) for j in seqs[i]] for i in seqs}
            seqs = {i: list(reversed(seqs[i])) for i in seqs}
        seqs = {i: "".join(seqs[i]) for i in seqs}
        human_length = len(seqs["homo_sapiens"])
        lengths = {i: len(seqs[i]) for i in seqs}
        #only keep the sequences whose length equals the human length (in other words, they all need to have the same gabs)
        seqs = {i: seqs[i] for i in seqs if len(seqs[i]) == human_length}
        #only proceed if at least min_species species survive
        if len(seqs) >= min_species:
            #make a single string with all the different species
            outstring = "%".join(["{0}|{1}|{2}|{3}".format(i, ",".join(current_dict[i]["coords"]), seqs[i], phase) for i in sorted(list(seqs.keys()))])
            outstring = outstring + "\n"
            return(outstring)
    return(None)

def MSA_filter_by_anatomy(input_file, output_file, version):
    '''
    Given an output file from get_MSA_concat_list, filter the CDSs based on whether the exon coordinates have been conserved.
    '''
    run_process(["perl", "MSA_CDSs.pl", version, input_file, output_file])

def filter_species(current_dict, mapping_tuple):
    '''
    Given the gabs associated to a particular species and a particular genomic region, check that the sequence region forms a contiguous block in that species.
    '''
    current_coords = [[int(j) for j in i] for i in current_dict["coords"]]
    #if the start coordinate is greater than the end coordinate, then
    #it means that the sequence is on the opposite strand from human
    #so now you use what you know of what the strand in human is
    #to figure out what the strand is in this species
    if current_coords[0][1] < current_coords[0][0]:
        local_strand = mapping_tuple[1]
        current_coords = [list(reversed(i)) for i in current_coords]
    else:
        local_strand = mapping_tuple[0]
    #if there is only one gab, the coords can be used as they are
    if len(current_dict["chrom"]) == 1:
        current_dict["coords"] = [current_dict["chrom"][0], str(current_coords[0][0]), str(current_coords[0][1]), local_strand]
    #if there are several gabs but they're all on the same chromosome, you need to make sure that they are contiguous
    elif len(list(set(current_dict["chrom"]))) == 1:
        sorted_coords = sort_coords(current_coords, 0)
        differences = [sorted_coords[i + 1][0] - sorted_coords[i][1] - 1 for i in range(len(sorted_coords) - 1)]
        if sum(differences) == 0:
            #the sequences need to be sorted too
            current_dict["seqs"] = [j for (i, j) in sorted(zip([k[0] for k in current_coords], current_dict["seqs"]), key = lambda pair: pair[0])]
            current_dict["coords"] = [current_dict["chrom"][0], str(sorted_coords[0][0]), str(sorted_coords[-1][1]), local_strand]
        else:
            return(None)
    else:
        return(None)
    #if the gabs are not all on the same chromosome or are not contiguous, discard the species
    return(current_dict)

def flush_tables(host, user, password):
    mq = _mysql.connect(host = host, user = user, passwd = password)
    mq.query("FLUSH TABLES;")
    mq.close()
        
def get_MSA(coords, method, species_set, query_species, version, force_strand = True):
    '''
    Get the genome alignments that overlap a particular sequence region.
    '''
    reverse = False
    if coords[6] == "-" and force_strand:
        reverse = True
    MSA = run_process(["perl", "MSA.pl", method, species_set, version, coords[0], coords[2], coords[3], query_species])
    MSA = MSA.split("|||")
    MSA = [i.split(">") for i in MSA if i]
    MSA = [[j.split("\n") for j in i if j] for i in MSA]

    MSA_dict = {}
    for gab in MSA:
        for species in gab: 
            name = species[0]
            temp_name = name.split("/")
            true_name = temp_name[0]
            coords = "-".join(temp_name[1:])
            if true_name not in MSA_dict:
                MSA_dict[true_name] = {}
            current_seq = "".join(species[1:]).upper()
            if reverse:
                current_seq = Seq(current_seq, IUPAC.unambiguous_dna)
                current_seq = current_seq.reverse_complement()
                current_seq = str(current_seq)
            MSA_dict[true_name][coords] = current_seq       
    return(MSA_dict)

def get_MSA_concat(coords, method, species_set, query_species, version, species_names, force_strand = True, reverse = False):
    '''
    Given a coordinate list, obtain and concatenate the orthologous sequence from a GW MSA. reverse only means that it'll check strand and then
    reverse the coordinates if necessary, if the coordinates are from the plus strand it doesn't do anything.
    '''
    antisense = False
    if coords[0][6] == "-":
        antisense = True
    #so the exons would be concatenated in the right order
    if antisense and force_strand and reverse:
        coords = [i for i in reversed(coords)]
    full_seq = {i: [] for i in species_names}
    expected_length = len(coords)
    for pos, exon in enumerate(coords):
        print(exon)
        current_dict = get_MSA(exon, method, species_set, query_species, version, force_strand = force_strand)
        print("\n\n")
        print(current_dict)
        for species in current_dict:
            if species in full_seq:
                print(species)
                current_coords = list(current_dict[species].keys())
                print(current_coords)
                if len(current_coords) > 1:
                    del full_seq[species]
                else:
                    pass

def get_MSA_concat_list(input_file, output_file, min_species):
    '''
    Given the sequences and coordinates from a bunch of raw MSA objects (as returned by get_MSA_gene_list),
    filter them to only keep ones where you have a contiguous CDS of the same length as in human.
    Write an output file where each row corresponds to one CDS region and has the orthologous sequence from each of the species.
    '''
    with open(input_file) as file, open(output_file, "w") as outfile:
        current_string = []
        lines = line_count(input_file)
        strand_mapping = {"human_sense": ("1", "-1"), "human_antisense": ("-1", "1")}
        print("Total: {0}.".format(lines))
        for pos, line in enumerate(file):
            if pos % 10000 == 0:
                print(pos)
            seqs = {}
            #that means we've gotten to a new human CDS record
            if line[0] == "%":
                line = line.rstrip("\n")
                global_name = line.split("|")
                strand = global_name[6]
                #it'll always fetch the gab from the sense strand (with regards to the query species) and
                #the other species from whatever aligns to the reference strand in the query species
                #so if the gene is antisense, everything has to be flipped
                if global_name[6] == "-":
                    mapping_tuple = strand_mapping["human_antisense"]
                else:
                    mapping_tuple = strand_mapping["human_sense"]
            #if it begins with neither a percentage sign nor an asterisk,
            #that means it must be a line of sequence
            #put all those in a list for each human CDS
            elif line[0] != "*":
                current_string.append(line)
            #you've reached the end of a CDS block
            else:
                #parse the stuff you've read in into a dictionary with the species as keys
                current_dict = get_species_dict(current_string)
                #loop over the different species
                for species in list(current_dict.keys()):
                    #check that the aligning region forms a contiguous block in this species
                    current_species_dict = filter_species(current_dict[species], mapping_tuple)
                    if current_species_dict:
                        current_dict[species] = current_species_dict
                    else:
                        del current_dict[species]
                outstring = concat_sequences(current_dict, min_species, global_name)
                if outstring:
                    outfile.write(outstring)
                current_string = []

def get_MSA_gene_list(coords, coords_file, method, species_set, version, query_species, MSA_file):
    '''
    Given a dictionary of lists of lists of CDS coordinates, retrieve the Compara MSAs.
    '''
    with open(coords_file, "w") as file:
        for trans in coords:
            for exon in coords[trans]:
                phase = exon[1]
                current_coords = exon[0]
                current_coords = [str(i) for i in current_coords]
                current_coords.append(str(phase))
                current_coords = "|".join(current_coords)
                file.write(current_coords)
                file.write("\n")
    remove_file(MSA_file)
    run_process(["perl", "MSA_list.pl", method, species_set, version, coords_file, query_species, MSA_file])
    with open(MSA_file) as file:
        string = "".join(file)
    string = re.sub("([a-z])\n([a-z])", "\\1\\2", string)
    with open(MSA_file, "w") as file:
        file.write(string)

def get_pairwise_alignment(coords, coords_file, query_species, other_species, version, output_file):
    '''
    Given a list of feature coordinates and two species, get the corresponding pairwise alignments from Compara.
    '''
    #write the coordinates to file in a way that can be read by the downstream perl script
    with open(coords_file, "w") as file:
        for feature in coords:
            feature = [str(i) for i in feature]
            feature = "|".join(feature)
            file.write(feature)
            file.write("\n")
    remove_file(output_file)
    #get the alignments from the database
    run_process(["perl", "pairwise_from_ensembl.pl", coords_file, query_species, other_species, output_file, version])
    #parse them from the output file produced by the perl script
    with open(output_file) as file:
        string = "".join(file)
    string = string.split("***")
    string = [(i.rstrip("\n")).lstrip("\n") for i in string]
    string = [i.split("|||") for i in string]
    string = [[j for j in i if len(j) > 0] for i in string]
    #getting rid of cases where there's multiple GABs
    string = flatten([i for i in string if len(i) == 1])
    #write alignments to a pretty FASTA
    with open(output_file, "w") as file:
        for feature in string:
            temp = feature.split("\n")
            name = temp[0]
            name = name.split("|")
            antisense = False
            if name[6] == "-":
                antisense = True
            try:
                alignments = [temp[2], temp[3]]
                alignments = [i.split(" ") for i in alignments]
                #covert to upper case
                alignments = [([j for j in i if j][1]).upper() for i in alignments]
                #only keep alignments with no ambiguous bases in either sequence
                if "N" not in alignments[0] and "N" not in alignments[1]:
                    #reverse complement, if necessary
                    if antisense:
                        alignments = [str(Seq(i, IUPAC.unambiguous_dna).reverse_complement()) for i in alignments]
                    file.write(">{0}\n".format("|".join(name)))
                    file.write("|".join(alignments))
                    file.write("\n")
            except IndexError:
                pass
            
def get_species_dict(current_string):
    '''
    Given the lines corresponding to one CDS region from a get_MSA_gene_list output file,
    parse it into a dictionary with the species as keys.
    '''
    #to get rid of the final three pipes
    current_string = current_string[:-1]
    current_string = "".join(current_string)
    #split it up into gabs
    current_string = current_string.split("|||")
    #split by species within each gab
    current_string = [i.split(">") for i in current_string if i]
    #split into lines within each species
    current_string = [[j.split("\n") for j in i if j] for i in current_string]
    current_dict = {}
    #loop over all the current gabs
    for gab in current_string:
        #loop over the species within the gab
        for species in gab:
            name = species[0]
            #the first element in all gabs other than the first one will be empty
            #(because of how split() works)
            if name:
                #separate into the name of the species and the coordinates
                temp_name = name.split("/")
                true_name = temp_name[0]
                chrom = temp_name[-2]
                coords = temp_name[-1]
                #make a dictionary where the keys are species and then for each species,
                #you have lists to store the coordinates and the sequences where each element
                #corresponds to one gab
                if true_name not in current_dict:
                    current_dict[true_name] = {"chrom": [], "coords": [], "seqs": []}
                current_dict[true_name]["coords"].append(coords.split("-"))
                current_dict[true_name]["chrom"].append(chrom)
                current_dict[true_name]["seqs"].append("".join(species[1:]))
    return(current_dict)


def MSA_names(method, species_set, version):
    '''
    Given a Compara WGA method, a species set name and an ensembl db version, get the names of all the species in the set.
    '''
    names = run_process(["perl", "MSA_names.pl", method, species_set, version])
    names = names.rstrip(",")
    names = names.split(",")
    return(names)

def parse_MSA_output(MSA_file):
    '''
    Given an output file from get_MSA_gene_list(), parse it into a pretty dictionary.
    '''
    print("\n")
    with open(MSA_file) as file:
        raw = file.read()
    raw = raw.split("*")
    raw = [i.split("|||") for i in raw if i]
    raw = [[j.split(">") for j in i if j != "\n"] for i in raw]
    raw = [[[((k.rstrip("\n")).lstrip("\n")).split("\n") for k in j if k != "\n"] for j in i] for i in raw]
    output = {}
    sorting_help = {}
    to_remove = []
    for align in raw:
        #cause the last one is often empty
        if align:
            current_coords = align[0][0][0]
            current_coords = current_coords.split("|")
            trans = current_coords[4]
            if trans not in output:
                output[trans] = []
            current_coords[2] = int(current_coords[2])
            current_coords[3] = int(current_coords[3])
            if trans not in sorting_help:
                sorting_help[trans] = []
            sorting_help[trans].append(current_coords[2])
            if len(align) != 1:
                to_remove.append(trans)
            else:
                current_dict = {}
                current_coords[0] = current_coords[0].lstrip("%")
                current_coords = current_coords[:-1]
                current_dict["coords"] = current_coords
                align = align[0]
                current_dict["species"] = {}
                for species in align[1:]:
                    species_name = species[0].split("/")[0]
                    current_dict["species"][species_name] = {}
                    current_dict["species"][species_name]["coords"] = species[0]
                    seq = ("".join(species[1:])).upper()
                    if current_coords[6] == "-":
                        seq = str(Seq(seq, IUPAC.unambiguous_dna).reverse_complement())
                    current_dict["species"][species_name]["seq"] = seq
                output[trans].append(current_dict)
    output = {i: output[i] for i in output if i not in to_remove}
    for trans in output:
        if output[trans][0]["coords"][6] == "+":
            output[trans] = [i for (j, i) in sorted(zip(sorting_help[trans], output[trans]), key=lambda pair: pair[0])]
        elif output[trans][0]["coords"][6] == "-":
            output[trans] = [i for (j, i) in sorted(zip(sorting_help[trans], output[trans]), reverse = True, key=lambda pair: pair[0])]
        else:
            print("Invalid strand!")
            raise Exception
    return(output)
                
                
