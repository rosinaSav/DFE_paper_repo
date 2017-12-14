from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import csv
##from goatools import obo_parser
from housekeeping import flatten, print_elements, make_dir, remove_file, run_process
import itertools
import multiprocessing as mp
from operator import itemgetter
import os
import my_stats as ms
import nucleotide_comp as nc
import numpy as np
import pandas
from pyfaidx import Fasta
import random
import re
import read_and_write as rw
import statistics
import sys
import time

class Feature_Set(object):

    def __init__(self, input_file, genome):
        self.features_file_name = input_file
        self.current_dir = os.getcwd()
        self.features_file_name_full = "{0}/{1}".format(self.current_dir,input_file)
        self.gene_pattern = re.compile("ENS\w*G\d+")
        self.transcript_pattern = re.compile("ENS\w*T\d+")
        self.genome = genome
        self.gene_name_pattern = re.compile("(?<=gene_name \")[\w\d\-\.\(\)\/]*(?=\")")

    def add_dataset(self,new_dataset_name,new_features_file_name,new_input_list = None, force_new = False):
        temp_fs = Feature_Set(new_features_file_name, self.genome)
        if force_new:
            temp_fs.create_dataset(new_dataset_name,input_list = new_input_list)
        else:
            try:
                temp_fs.set_dataset(new_dataset_name,new_input_list, verbose = False) 
            except FileNotFoundError:
                temp_fs.create_dataset(new_dataset_name,input_list = new_input_list)
        self.names.extend(temp_fs.names)
        self.features.extend(temp_fs.features)
        self.dataset = "{0}_{1}".format(self.dataset,temp_fs.dataset)
        self.dataset_components.extend([temp_fs.dataset for i in temp_fs.names])
        print("The data from the dataset {0} has been added to the current dataset.".format(new_dataset_name))

    def add_families(self, families):
        try:
            current_max = len(self.families.keys())
        except AttributeError:
            current_max = 0
            self.families = {}
        for i in range(len(families)):
            self.families[i+current_max] = families[i]

    def alternative_splicing(self, gene_name_dict, CDS):
        results = {}
        for name in gene_name_dict:
            previous_coords = []
            found = False
            alt = False
            for trans in gene_name_dict[name]:
                current_CDS = CDS[trans[4]]
                current_CDS
                if current_CDS:
                    found = True
                    current_coords = [(i[0][2], i[0][3]) for i in current_CDS]
                    if previous_coords:
                        if current_coords != previous_coords:
                            alt = True
                            break
                    previous_coords = current_coords.copy()
            if not found:
                results[name] = None
            elif alt:
                results[name] = 1
            else:
                results[name] = 0
        return(results)

    def average_over_families(self, input_dict, remove_empty = False, use_nan = False):
        try:
            flat_families = flatten(self.families.values())
            output_dict = {}
            for i in input_dict:
                if i not in flat_families:
                    output_dict[i] = input_dict[i]
            for i in range(len(self.families)):
                family_name = "Family{0}".format(i)
                family_values = []
                for j in self.families[i]:
                    try:
                        family_values.append(input_dict[j])
                    except KeyError:
                        pass
                family_values = np.array(family_values)               
                family_values = [k for k in family_values if k != None]
                if family_values:
                    output_dict[family_name] = np.mean(family_values)
                else:
                    if not remove_empty:
                        if use_nan:
                            output_dict[family_name] = np.nan
                        else:
                            output_dict[family_name] = None
            return(output_dict)
        except AttributeError:
            print("No families associated to dataset!")
            raise AttributeError

    def average_over_families_2d(self, input_dict, remove_nans = False, nan_for_None = False, require_a_member = True, remove_empty = False):
        if input_dict == {}:
            return(input_dict)
        try:
            flat_families = flatten(self.families.values())
            output_dict = {}
            for i in input_dict:
                if i not in flat_families:
                    output_dict[i] = input_dict[i]
            #this is in case you need it with nans
                    #you're using i because you're taking the row length from the last item you saw
            row_length = len(input_dict[i])
            found_something = False
            for i in range(len(self.families)):
                family_name = "Family{0}".format(i)
                family_values = []
                for j in self.families[i]:
                    try:
                        family_values.append(input_dict[j])
                    except KeyError:
                        pass                
                if family_values:
                    found_something = True
                    family_values = np.array(family_values)
                    if not remove_nans:
                        output_dict[family_name] = np.mean(family_values, axis = 0)
                    else:
                        output_dict[family_name] = np.nanmean(family_values, axis = 0)
                else:
                    if not remove_empty:                        
                        if not nan_for_None:
                            output_dict[family_name] = None
                        else:
                            output_dict[family_name] = np.array([np.nan for i in range(row_length)])
            if not found_something:
                print("Nothing!")
                if require_a_member:
                    print("None of the observations belong to families!")
                    raise Exception
            return(output_dict)
        except AttributeError:
            print("No families associated to dataset!")
            raise AttributeError
       
    def bin_values(self, values_dict, bin_number, labels):
        values = [[key, value] for key, value in values_dict.items()]
        names = [i[0] for i in values]
        values = [i[1] for i in values]
        indices, bins = pandas.qcut(values, bin_number, labels = labels, retbins = True)
        output_dict = {}
        for i in labels:
            output_dict[i] = []
        for i,j in enumerate(names):
            output_dict[indices[i]].append(j)
        np.set_printoptions(suppress = True)
        return(output_dict, bins)

    def convert_between_ENST_and_ENSG(self, identifier, gene_name_dict, dest, families = False):
        #if destination is ENST, expect a gene_name_dict, otherwise expects transcripts
        if families:
            if identifier[:6] == "Family":
                return(identifier)
        if dest == "ENST":
            if gene_name_dict[identifier][0][:3] != "ENS":
                print("If dest is ENST, input gene_name_dict!")
                raise Exception
            curr_trans = gene_name_dict[identifier]
            if len(curr_trans) > 1:
                print("The identifier {0} corresponds to more than 1 transcript ID!".format(identifier))
                raise Exception
            else:
                return(curr_trans[0])
        elif dest == "ENSG":
            if gene_name_dict[identifier][0][:3] == "ENS":
                print("If dest is ENSG, input transcripts dictionary!")
                raise Exception
            return(gene_name_dict[identifier][5])
        else:
            print("Invalid destination! Valid destinations are \"ENST\" and \"ENSG\"!")
            raise Exception             

    def convert_families_to_ENST(self, families, transcripts):
        #list of ENST identifiers
        keys_list = list(transcripts.keys())
        found_counter = 0
        #loop over families
        for i in range(len(families)):
            #loop over the genes in a family
            for j in range(len(families[i])):
                counter = 0
                found = False
                #loop through all of the transcript identifiers and check whether they correspond
                #to the current gene
                while not found and counter < len(keys_list):
                    if transcripts[keys_list[counter]][5] == families[i][j]:
                        families[i][j] = keys_list[counter]
                        found = True
                        found_counter = found_counter + 1
                    counter = counter + 1
        if len(flatten(families)) != found_counter:
            print("Family identifiers could not be converted!")
            sys.exit()
        return(families)

    def convert_to_absolute_coordinates(self, feature, position, bed_input = False):
        '''
        Given a feature (1-based) and a relative position within the feature (0-based), convert that to
        absolute feature coordinates.
        '''
        if bed_input:
            strand_col = 5
            start_col = 1
            end_col = 2
            shift = 0
            for pos, line in enumerate(bed_input):
                bed_input[pos][start_col] = int(line[start_col])
                bed_input[pos][end_col] = int(line[start_col])
        else:
            strand_col = 6
            start_col = 2
            end_col = 3
            shift = 1
        feature_length = self.get_lengths({"foo": feature})["foo"]
        if feature[0][strand_col] == "+":
            pass
        elif feature[0][strand_col] == "-":
            position = feature_length - position - shift
        else:
            print("Invalid strand information!")
            raise Exception
        cum_length = 0
        #loop over the exons or whatever
        for elem in feature:
            #keep track of the cumulative length of all the exons you've seen
            new_cum_length = cum_length + (elem[end_col] - elem[start_col] + shift)
            #if the cumulative length is greater than position,
            #it must be within that exon
            if new_cum_length > position:
                return(elem[start_col] + (position - cum_length))
            cum_length = new_cum_length
        print("Position doesn't map to feature!")
        raise Exception

    def create_dataset(self, dataset_name, split_size = 100 * 1024 * 1024, input_list = None, filter_trans = True, input_type = None):
        if input_type == None:
            input_type = "transcript"
        #I stole the outline of this chunk of code from http://aamirhussain.com/2013/10/02/parsing-large-csv-files-in-python/
        if filter_trans:
            if input_list:
                self.names = input_list
            else:
                try:
                    self.names = rw.read_names("{0}_folder/{0}_trans_identifiers_wo_low_omega.txt".format(dataset_name))
                except FileNotFoundError:
                    print("No transcript identifiers associated to dataset {0}!".format(dataset_name))
        else:
            self.names = "foo"
        file_size = os.path.getsize(self.features_file_name)
        pool = mp.Pool(int(os.cpu_count()/2))
        results = []
        #you go over your file, chunk by chunk
        #at each chunk, you create a process to filter that chunk and add it to the pool
        with open(self.features_file_name) as features_file:
            #this is to skip the comment lines in the beginning
            for i in range(5):
                features_file.readline()
            cursor = features_file.tell()
            for chunk in range((file_size // split_size) + 1):
                if cursor + split_size > file_size:
                    end = file_size
                else:
                    end = cursor + split_size
                features_file.seek(end)
                #this is to ensure that each chunk ends with a full line
                features_file.readline()
                end = features_file.tell()
                process = pool.apply_async(filter_features, args=[self.features_file_name, self.names, cursor, end, self.gene_pattern, self.transcript_pattern, self.gene_name_pattern, filter_trans, input_type])
                results.append(process)
                cursor = end
        pool.close()
        pool.join()
        output_file_name = "{0}_{1}_dataset_features.bed".format(self.features_file_name[:-4],dataset_name)
        output_file = open(output_file_name, "w")
        output_writer = csv.writer(output_file, delimiter = "\t")
        output_writer.writerow(["# Lines corresponding to {0} in {1}.".format(dataset_name,self.features_file_name)])
        self.features = []
        names = []
        for i in results:
            process_file_result = i.get()
            for j in process_file_result:
                output_writer.writerow(j)
                self.features.append(j)
                names.append(j[4])
        output_file.close()
        self.dataset_file_name = output_file_name
        self.dataset_file_name_full = "{0}/{1}".format(self.current_dir,output_file_name)
        self.dataset = dataset_name
        if (not filter_trans) or (input_type == "gene"):
            names = list(set(names))
            self.names = names
            make_dir("{0}_folder".format(self.dataset))
            names_file = open("{0}_folder/{0}_trans_identifiers_wo_low_omega.txt".format(self.dataset), "w")
            for i in self.names:
                names_file.write("{0}\n".format(i))
            names_file.close()
        self.dataset_components = [self.dataset for i in self.names]
        print("Created dataset {0}.".format(dataset_name))

    def get_CDS(self):
        CDS_dict = {}
        stop_codon_dict = self.get_stop_codons()
        for i in self.names:
            CDS_dict[i] = []
        for i in self.features:
            if i[1] == "CDS":
                CDS_dict[i[4]].append(i)
        for i in CDS_dict:
            if CDS_dict[i]:
                current_stop_coords = stop_codon_dict[i]
                current_stop_coords_loc_only = [[j[2],j[3]] for j in current_stop_coords]
                CDS_dict[i] = [j for j in CDS_dict[i] if [j[2],j[3]] not in current_stop_coords_loc_only]               
                if CDS_dict[i][0][6] == "+":
                    CDS_dict[i] = sort_coords(CDS_dict[i],2)
                elif CDS_dict[i][0][6] == "-":
                    CDS_dict[i] = sort_coords(CDS_dict[i],3, reverse = True)                   
                lengths = [j[3] - j[2] + 1 for j in CDS_dict[i]]
                cumul_lengths = [sum(lengths[0:j]) for j in range(0, len(CDS_dict[i]))]
                phases = [j%3 for j in cumul_lengths]
                CDS_dict[i] = [[CDS_dict[i][j], phases[j]] for j in range(len(phases))]
        return(CDS_dict)

    def get_CDS_gene_ratio(self, CDS, transcripts, UTRs):
        ratios = {}
        CDS_lengths = self.get_lengths(CDS, CDS = True)
        transcript_lengths = self.get_lengths(transcripts, CDS = False)
        for trans in self.names:
            if len(CDS[trans]) == 0 or trans not in transcripts:
                ratios[trans] = None
            else:
                current_CDS_length = CDS_lengths[trans]
                current_transcript_length = transcript_lengths[trans]
                if trans in UTRs:
                    sum_utr = sum([i[3] - i[2] + 1 for i in UTRs[trans][5]])
                    sum_utr = sum_utr + sum([i[3] - i[2] + 1 for i in UTRs[trans][3]])
                else:
                    sum_utr = 0
                ratios[trans] = current_CDS_length/(current_transcript_length - sum_utr)
        return(ratios)

    def get_const_exons(self, gene_name_dict, CDS):
        '''
        Given a set of CDSs, only keep ones that are constitutively spliced (= present in all annotated isoforms).
        '''
        results = {}
        for name in gene_name_dict:
            all_trans = [CDS[i[4]] for i in gene_name_dict[name] if len(CDS[i[4]]) > 0]
            flat_exons = flatten(all_trans)
            flat_exon_coords = [(i[0][0], i[0][2], i[0][3]) for i in flat_exons]
            trans_number = len(all_trans)
            const_exons = [flat_exons[i] for i in range(len(flat_exons)) if flat_exon_coords.count(flat_exon_coords[i]) == trans_number]
            for trans in gene_name_dict[name]:
                current_trans = trans[4]
                results[current_trans] = [i for i in CDS[current_trans] if i in const_exons]
        return(results)

    def distance_to_nearest_junction(self, CDS, exons, stop_codons, position, sequence = None, verbose = False):
        '''
        Given a position in the CDS, return the distance to the nearest splice junction.
        '''
        #this is just to make it compatible with the format of a CDS feature
        if stop_codons:
            stop_codons = [[i, 0] for i in stop_codons]
        CDS.extend(stop_codons)
        if CDS[0][0][7] == "-":
            CDS.reverse()
        min_dist = self.distance_to_nearest_junction_core(CDS, exons, position, verbose = verbose)
        if not sequence:
            return(min_dist)
        base_pos = [pos for pos, i in enumerate(sequence) if (i == sequence[position]) and (pos != position)]
        sim_dists = []
        for sim_pos in base_pos:
            current_sim_dist = self.distance_to_nearest_junction_core(CDS, exons, sim_pos, verbose = False)
            sim_dists.append(current_sim_dist)
        p = ms.calc_eff_p(min_dist, sim_dists, greater = False)
        return(min_dist, p)
    
    def distance_to_nearest_junction_core(self, CDS, exons, position, verbose = False):
        '''
        Core for distance_to_nearest_junction.
        '''
        position = self.convert_to_absolute_coordinates([i[0] for i in CDS], position)
        for pos, exon in enumerate(exons):
            if (position >= exon[2]) and (position <= exon[3]):
                dist = [position - exon[2], exon[3] - position]
                if pos == 0:
                    return(dist[1])
                elif pos == (len(exons) - 1):
                    return(dist[0])
                return(min(dist))
        print("Can't match to junction!")
        raise Exception

    def distance_to_TSS(self, CDS, transcripts):
        distance_dict = {}
        for name in self.names:
            if transcripts[name][6] == "+":
                distance_dict[name] = CDS[name][0][0][2] - transcripts[name][2]
            elif transcripts[name][6] == "-":
                distance_dict[name] = transcripts[name][3] - CDS[name][0][0][3]
        return(distance_dict)

    def get_exon_beginnings(self, exons, CDS, file_prefix = None, write_to_fasta = False):
        #filter out genes with fewer than three exons
        ids_to_keep = []
        for i in exons:
            if len(exons[i]) >= 3:
                ids_to_keep.append(i)
        for i in list(exons.keys()):
            if i not in ids_to_keep:
                del exons[i]
            else:
                #only leave internal exons
                exons[i] = exons[i][1:-1]
                #remove exons that are too short to get 69 + 69 + 4 bp
                exons[i] = [j for j in exons[i] if (j[3] - j[2] + 1) >= 142]
        records_to_keep = []
        counter = 0
        for i in CDS.keys():
            if counter % 1000 == 0:
                print(counter)
            counter = counter + 1
            if i in exons.keys():
                exon_coords = [(j[2], j[3]) for j in exons[i]]
                #this is to get rid of exons that are only partially coding
                CDS[i] = [j for j in CDS[i] if (j[0][2], j[0][3]) in exon_coords]
                #trim remaining exons so that they would be in phase
                CDS[i] = [trim_sequence_coords(j[0], j[1]) for j in CDS[i]]
                if CDS[i]:
                    if CDS[i][0][6] == "+":
                        #also transform to BED
                        records = [[j[0], j[2] - 1, j[2] + 69 - 1, j[4], "100", j[6]] for j in CDS[i]]
                    else:
                        records = [[j[0], j[3] - 69, j[3], j[4], "100", j[6]] for j in CDS[i]]
                    records_to_keep.extend(records)
        if not file_prefix:
            return(records_to_keep)
        else:
            rw.write_to_csv(records_to_keep, "{0}_uf_filt69.bed".format(file_prefix), "\t")
            if write_to_fasta:
                fasta_from_intervals("{0}_uf_filt69.bed".format(file_prefix), "{0}_uf_filt69.fasta".format(file_prefix), self.genome, force_strand = True)
        
    def get_exon_cores_and_flanks(self, exons, CDS, file_prefix = None, write_to_fasta = False):
        #filter out genes with fewer than three exons
        ids_to_keep = []
        for i in exons:
            if len(exons[i]) >= 3:
                ids_to_keep.append(i)
        for i in list(exons.keys()):
            if i not in ids_to_keep:
                del exons[i]
            else:
                #only leave internal exons
                exons[i] = exons[i][1:-1]
                #remove exons that are too short to get 3*69 + 4 bp
                exons[i] = [j for j in exons[i] if (j[3] - j[2] + 1) >= 211]
        records_to_keep_uf = []
        records_to_keep_c = []
        records_to_keep_df = []
        counter = 0
        for i in CDS.keys():
            if counter % 1000 == 0:
                print(counter)
            counter = counter + 1
            if i in exons.keys():
                exon_coords = [(j[2], j[3]) for j in exons[i]]
                #this is to get rid of exons that are only partially coding
                CDS[i] = [j for j in CDS[i] if (j[0][2], j[0][3]) in exon_coords]
                #trim remaining exons so that they would be in phase
                CDS[i] = [trim_sequence_coords(j[0], j[1]) for j in CDS[i]]
                if CDS[i]:
                    if CDS[i][0][6] == "+":
                        #also transform to BED
                        uf_records = [[j[0], j[2] - 1, j[2] + 69 - 1, j[4], "100", j[6]] for j in CDS[i]]
                        df_records = [[j[0], j[3] - 69, j[3], j[4], "100", j[6]] for j in CDS[i]]
                        remainder = [(df_records[j][1] - uf_records[j][2] - 69)/2 for j in range(len(CDS[i]))]
                        shift = [uf_records[j][2] + remainder[j] if remainder[j]%3 == 0 else uf_records[j][2] + ((remainder[j]//3) * 3) + 3 for j in range(len(CDS[i]))]
                    else:
                        uf_records = [[j[0], j[3] - 69, j[3], j[4], "100", j[6]] for j in CDS[i]]
                        df_records = [[j[0], j[2] - 1, j[2] + 69 - 1, j[4], "100", j[6]] for j in CDS[i]]
                        remainder = [(uf_records[j][1] - df_records[j][2] - 69)/2 for j in range(len(CDS[i]))]
                        shift = [df_records[j][2] + remainder[j] if remainder[j]%3 == 0 else df_records[j][2] + ((remainder[j]//3) * 3) + 3 for j in range(len(CDS[i]))]
                    c_records = [[CDS[i][j][0], int(shift[j]), int(shift[j] + 69), CDS[i][j][4], "100", CDS[i][j][6]] for j in range(len(CDS[i]))]                    
                    records_to_keep_uf.extend(uf_records)
                    records_to_keep_c.extend(c_records)
                    records_to_keep_df.extend(df_records)
        if not file_prefix:
            return({"upstream": records_to_keep_uf, "core": records_to_keep_c, "downstream": records_to_keep_df})
        else:
            rw.write_to_csv(records_to_keep_uf, "{0}_uf.bed".format(file_prefix), "\t")
            rw.write_to_csv(records_to_keep_df, "{0}_df.bed".format(file_prefix), "\t")
            rw.write_to_csv(records_to_keep_c, "{0}_c.bed".format(file_prefix), "\t")
            if write_to_fasta:
                fasta_from_intervals("{0}_uf.bed".format(file_prefix), "{0}_uf.fasta".format(file_prefix), self.genome, force_strand = True)
                fasta_from_intervals("{0}_df.bed".format(file_prefix), "{0}_df.fasta".format(file_prefix), self.genome, force_strand = True)
                fasta_from_intervals("{0}_c.bed".format(file_prefix), "{0}_c.fasta".format(file_prefix), self.genome, force_strand = True)

    def get_exons_coding(self, exons, CDS):
        #filter out genes with fewer than three exons
        ids_to_keep = []
        for i in exons:
            if len(exons[i]) >= 3:
                ids_to_keep.append(i)
        for i in list(exons.keys()):
            if i not in ids_to_keep:
                del exons[i]
            else:
                #only leave internal exons
                exons[i] = exons[i][1:-1]
        counter = 0
        to_keep = {}
        for i in CDS.keys():
            if counter % 1000 == 0:
                print(counter)
            counter = counter + 1
            if i in exons.keys():
                exon_coords = [(j[2], j[3]) for j in exons[i]]
                #this is to get rid of exons that are only partially coding
                CDS[i] = [j for j in CDS[i] if (j[0][2], j[0][3]) in exon_coords]
                #trim remaining exons so that they would be in phase
                CDS[i] = [trim_sequence_coords(j[0], j[1]) for j in CDS[i]]
                to_keep[i] = CDS[i]
        return(to_keep)

    def get_exon_numbers(self, exons):
        exon_number_dict = {}
        for i in exons:
            if len(exons[i]) > 0:
                exon_number_dict[i] = len(exons[i])
            else:
                exon_number_dict[i] = None
        return(exon_number_dict)

    def get_exon_ranks(self, exons, reverse = False, limit = None):
        ranks_dict = {}
        if limit:
            max_rank = limit
        else:
            exon_numbers = [len(i) for i in exons]
            max_rank = max(exon_numbers) - 1
        if reverse:
            [exons[i].reverse() for i in exons]
        for rank in range(max_rank):
            ranks_dict[rank] = [exons[i][rank] for i in exons if len(exons[i]) > rank]
        return(ranks_dict)

    def get_exons(self):
        exons_dict = {}
        for i in self.names:
            exons_dict[i] = []
        for i in self.features:
            if i[1] == "exon":
                exons_dict[i[4]].append(i)
        for i in self.names:
            if exons_dict[i]:
                if exons_dict[i][0][6] == "+":
                    exons_dict[i] = sort_coords(exons_dict[i], 2)
                elif exons_dict[i][0][6] == "-":
                    exons_dict[i] = sort_coords(exons_dict[i], 3, reverse = True)                    
        return(exons_dict)

    def get_exon_sizes_by_rank(self, exons, rank):
        exon_sizes = {}
        for i in self.names:
            if len(exons[i]) <= rank:
                exon_sizes[i] = None
            else:
                if exons[i][0][6] == "+":
                    current_coords = sort_coords(exons[i], 2, reverse = False)
                elif exons[i][0][6] == "-":
                    current_coords = sort_coords(exons[i], 3, reverse = True)
                else:
                    print("Invalid strand information.")
                exon_sizes[i] = current_coords[rank][3] - current_coords[rank][2] + 1
        return(exon_sizes)

    def get_family_from_id(self,name):
        for family, ids in self.families.items():
            if name in ids:
                return(family)
        return(None)

    def get_flanking_intron_sizes(self, exons):
        intron_sizes = {}
        for i in self.names:
            if len(exons[i]) < 2:
                intron_sizes[i] = None
            else:
                intron_sizes[i] = {}
                if exons[i][0][6] == "+":
                    current_intron_sizes = [exons[i][j+1][2] - exons[i][j][3] - 1 for j in range(len(exons[i]) - 1)]
                elif exons[i][0][6] == "-":
                    current_intron_sizes = [exons[i][j][2] - exons[i][j+1][3] - 1 for j in range(len(exons[i]) - 1)]
                if len([k for k in current_intron_sizes if k < 1]) != 0:
                    print("Returned negative intron size!")
                    sys.exit()
                intron_sizes[i]["upstream"] = [None]
                intron_sizes[i]["upstream"].extend(current_intron_sizes)
                intron_sizes[i]["downstream"] = current_intron_sizes
                intron_sizes[i]["downstream"].append(None)
        return(intron_sizes)

    def get_gene_name_dict(self, transcripts):
        gene_name_dict = {}
        transcripts = list(transcripts.values())
        for i in transcripts:
            if i[5] in gene_name_dict.keys():
                gene_name_dict[i[5]].append(i[4])
            else:
                gene_name_dict[i[5]] = [i[4]]
        return(gene_name_dict)

    def get_intron_densities(self, exons, CDS):
        intron_densities = {}
        for i in self.names:
            if not exons[i] or not CDS[i]:
                intron_densities[i] = None
            else:
                current_intron_number = len(exons[i]) - 1
                current_CDS_lengths = [j[0][3] - j[0][2] + 1 for j in CDS[i]]
                current_CDS_length = sum(current_CDS_lengths)
                intron_densities[i] = current_intron_number/current_CDS_length
        return(intron_densities)

    def get_introns(self):
        introns_dict = {}
        for i in self.names:
            introns_dict[i] = []
        for i in self.features:
            if i[1] == "intron":
                introns_dict[i[4]].append(i)
        for i in self.names:
            if introns_dict[i]:
                if introns_dict[i][0][6] == "+":
                    introns_dict[i] = sort_coords(introns_dict[i], 2)
                elif introns_dict[i][0][6] == "-":
                    introns_dict[i] = sort_coords(introns_dict[i], 3, reverse = True)                    
        return(introns_dict)

    def get_introns_from_coords(self, exons, upstream = None, downstream = None):
        if upstream and downstream:
            print("You cannot restrict both the upstream and downstream distance to the exon!")
            raise Exception
        introns = {}
        for trans in exons:
            upstream_local = upstream
            downstream_local = downstream
            introns[trans] = []
            if exons[trans][0][6] == "-":
                upstream_local = downstream
                downstream_local = upstream
            if len(exons[trans]) > 1:
                if exons[trans][0][2] > exons[trans][1][2]:
                    exons[trans] = [i for i in reversed(exons[trans])]
            starts = [i[3] + 1 for i in exons[trans][:-1]]
            ends = [i[2] - 1 for i in exons[trans][1:]]
            if upstream_local:
                starts = [i - upstream_local + 1 for i in ends]
            if downstream_local:
                ends = [i + downstream_local - 1 for i in starts]
            coords = zip(starts, ends)
            template = exons[trans][0]
            for coord in coords:
                current_intron = template.copy()
                current_intron[1] = "intron"
                current_intron[2] = coord[0]
                current_intron[3] = coord[1]
                introns[trans].append(current_intron)
        return(introns)

    def get_lengths(self, feature, CDS = False):
        lengths_dict = {}
        if not CDS:
            #this is in case it's, for instance, a transcript feature where the values are lists
            #rather than lists of lists
            if type(feature[list(feature.keys())[0]][0]) != list:
                feature = {i:[feature[i]] for i in feature}
            for i in feature:
                length = 0
                for j in feature[i]:
                    length = length + int(j[3]) - int(j[2]) + 1
                lengths_dict[i] = length
        else:
            for i in feature:
                length = 0
                for j in feature[i]:
                    length = length + int(j[0][3]) - int(j[0][2]) + 1
                #this is to include the stop as well
                length = length + 3
                lengths_dict[i] = length
        return(lengths_dict)

    def get_mean_intron_sizes(self, transcripts, exons):
        mean_intron_sizes = {}
        for i in self.names:
            transcript_length = transcripts[i][3] - transcripts[i][2] + 1
            exonic_length = sum([(j[3] - j[2] + 1) for j in exons[i]])
            intronic_length = transcript_length - exonic_length
            if len(exons[i]) > 1:
                mean_intron_sizes[i] = intronic_length/(len(exons[i]) - 1)
                if mean_intron_sizes[i] < 0:
                    print("Returned negative mean intron size!")
                    print(exons[i])
                    print(transcripts[i])
                    sys.exit()
            elif len(exons[i]) == 1:
                mean_intron_sizes[i] = 0
            else:
                mean_intron_sizes[i] = None
        return(mean_intron_sizes)

    def get_singletons(self):
        flat_families = flatten(list(self.families.values()))
        singletons = [i for i in self.names if i not in flat_families]
        return(singletons)

    def get_start_codons(self):
        start_codon_dict = {}
        for i in self.names:
            start_codon_dict[i] = []
        for i in range((len(self.features)-1),-1,-1):
            if self.features[i][1] == "start_codon":
                start_codon_dict[self.features[i][4]].append(self.features[i])
        return(start_codon_dict)

    def get_stop_codons(self):
        stop_codon_dict = {}
        for i in self.names:
            stop_codon_dict[i] = []
        for i in range((len(self.features)-1),-1,-1):
            if self.features[i][1] == "stop_codon":
                stop_codon_dict[self.features[i][4]].append(self.features[i])
        return(stop_codon_dict)

    def get_transcripts(self, obligatory_coords = True):
        transcript_dict = {}
        for i in range((len(self.features)-1),-1,-1):
            if self.features[i][1] == "transcript":
                transcript_dict[self.features[i][4]] = self.features[i]
        if obligatory_coords:
            if len(transcript_dict.keys()) != len(self.names):
                print("Not all transcript identifiers have associated transcript coordinates!")
                sys.exit()
        return(transcript_dict)

    def get_UTRs(self, CDS):
        UTRs_dict = {}
        for i in self.names:
            UTRs_dict[i] = {5: [], 3: []}
        for i in self.features:
            if i[1] == "UTR":
                if i[6] == "+":
                    if CDS[i[4]][0][0][2] > i[3]:
                        UTRs_dict[i[4]][5].append(i)
                    elif CDS[i[4]][-1][0][3] < i[2]:
                        UTRs_dict[i[4]][3].append(i)
                    else:
                        print("Nonsensical UTR coordinates!")
                        sys.exit()
                if i[6] == "-":
                    if CDS[i[4]][0][0][3] < i[2]:
                        UTRs_dict[i[4]][5].append(i)
                    elif CDS[i[4]][-1][0][2] > i[3]:
                        UTRs_dict[i[4]][3].append(i)
                    else:
                        print("Nonsensical UTR coordinates!")
                        sys.exit()                        
        return(UTRs_dict)

    def get_UTRs_new(self, exons, CDS):
        UTRs = {}
        for trans in exons:
            if trans in CDS:
                UTRs[trans] = {}
                antisense = False
                if exons[trans][0][6] == "-":
                    antisense = True
                if len(exons[trans]) > 1:
                    if exons[trans][0][2] > exons[trans][1][2]:
                        exons[trans] = list(reversed(exons[trans]))
                        CDS[trans] = list(reversed(CDS[trans]))
                exon_starts = [i[2] for i in exons[trans]]
                exon_ends = [i[3] for i in exons[trans]]
                CDS_start = CDS[trans][0][0][2]
                CDS_end = CDS[trans][-1][0][3]
                upstream = []
                downstream = []
                rank = 0
                while CDS_start > exon_starts[rank]:
                    template = exons[trans][rank].copy()
                    template[1] = "UTR"
                    if exon_ends[rank] < CDS_start:
                        upstream.append(template)
                    else:
                        template[3] = CDS_start - 1
                        upstream.append(template)
                    rank = rank + 1
                    if rank >= len(exon_starts):
                        break
                rank = -1
                while CDS_end < exon_ends[rank]:
                    template = exons[trans][rank].copy()
                    template[1] = "UTR"
                    if exon_starts[rank] > CDS_end:
                        downstream.append(template)
                    else:
                        template[2] = CDS_end + 1
                        downstream.append(template)
                    rank = rank - 1
                    if -rank > len(exon_starts):
                        break
                if antisense:
                    UTRs[trans][5] = downstream
                    UTRs[trans][3] = list(reversed(upstream))
                else:
                    UTRs[trans][5] = upstream
                    UTRs[trans][3] = list(reversed(downstream))
        return(UTRs)
       
    def load_sequences(self, feature, file = None):
        if file:
            names, sequences = rw.read_fasta(file)
        else:
            names, sequences = rw.read_fasta("{0}_{1}_{2}.fasta".format(self.features_file_name[:-4], self.dataset, feature))
        feature_dict = {}
        name_pattern = re.compile("E[A-Z0-9]*(?=_)")
        for i in self.names:
            feature_dict[i] = []
        for i,j in enumerate(names):
            current_name = re.findall(name_pattern, j)
            current_name = current_name[0]
            feature_dict[current_name].append(sequences[i])
        return(feature_dict)

    def map_from_regions_to_exons(self, bed_record, exons):
        '''
        Take a bed-record and a list of 1-based exon coordinates from the same
        transcript and figure out which exon the bed record region belongs to.
        '''
        if exons[0][0] != bed_record[0]:
            print("Chromosome IDs don't match up!")
            print(exons[0][0])
            print(bed_record[0])
            raise Exception
        elif exons[0][6] != bed_record[5]:
            print("Strands don't match up!")
            print(exons[0][6])
            print(bed_record[5])
            raise Exception
        #convert to base 1
        bed_start = int(bed_record[1]) + 1
        bed_end = int(bed_record[2])
        matches = []
        for pos, exon in enumerate(exons):
            if bed_end > exon[3]:
                pass
            elif bed_start < exon[2]:
                pass
            else:
                matches.append(pos)
        if len(matches) == 0:
            return(None)
        elif len(matches) > 1:
            print("The bed record matches to more than one exon!")
            print(bed_record)
            print(exons)
            raise Exception
        else:
            return(matches[0])      

    def pick_random_members(self, seed = None):
        if seed:
            random.seed(a = seed)
        picked = []
        flat_families = flatten(list(self.families.values()))
        for i in self.names:
            if i not in flat_families:
                picked.append(i)
        for i in self.families:
            picked.append(random.choice(self.families[i]))
        expected_size = len(self.names) - len(flat_families) + len(self.families)
        if len(picked) != expected_size:
            print("Problem picking random members!")
            print("Expected {0}, got {1}!".format(expected_size, len(picked)))
            sys.exit()
        random.seed(a = None)
        return(picked)

    def pick_shortest_pc_trans(self, gene_name_dict, CDS):
        picked_transcripts = []
        for i in gene_name_dict:
            current_CDS = {}
            for j in gene_name_dict[i]:
                current_CDS[j] = CDS[j]
            current_pc_trans = [j for j in current_CDS if current_CDS[j]]
            if not current_pc_trans:
                pass
            elif len(current_pc_trans) == 1:
                picked_transcripts.append(current_pc_trans[0])
            else:
                current_pc_numbers = np.array([len(current_CDS[j]) for j in current_pc_trans])
                shortest_trans = current_pc_trans[np.nonzero(current_pc_numbers == np.min(current_pc_numbers))[0][0]]
                picked_transcripts.append(shortest_trans)
        return(picked_transcripts)

    def set_dataset(self, dataset_name, input_list = None, verbose = True):
        dataset_file_name = "{0}_{1}_dataset_features.bed".format(self.features_file_name[:-4],dataset_name)
        try:
            self.features = rw.read_many_fields(dataset_file_name, "\t")
            self.features = [[i[0], i[1], int(i[2]), int(i[3]), i[4], i[5], i[6], i[7]] for i in self.features[1:]]
            self.dataset = dataset_name
            self.dataset_file_name = "{0}_{1}_dataset_features.bed".format(self.features_file_name[:-4],dataset_name)
            self.dataset_file_name_full = "{0}/{1}".format(self.current_dir,self.dataset_file_name)
            if input_list:
                self.names = input_list
            else:
                try:
                    self.names = rw.read_names("{0}_folder/{0}_trans_identifiers_wo_low_omega.txt".format(dataset_name))
                except FileNotFoundError:
                    names = [i[4] for i in self.features]
                    self.names = list(set(names))
            self.dataset_components = [self.dataset for i in self.names]
            if verbose:
                print("Dataset has been set to {0}.".format(dataset_name))
        except FileNotFoundError:
            if verbose:
                print("Dataset {0} does not exist!".format(dataset_name))
            else:
                raise FileNotFoundError

    def trans_per_gene(self, gene_name_dict):
        picked = []
        for gene in gene_name_dict:
            to_add = random.choice(gene_name_dict[gene])
            picked.append(to_add[4])
        return(picked)

    def write_CDS(self, coords, file = None, input_list = None):
        if not file:
            file = "{0}_{1}_CDS.fasta".format(self.features_file_name[:-4], self.dataset)
        output_file = open(file, "w")
        if not input_list:
            input_list = list(coords.keys())
        for i in input_list:
            for j in range(len(coords[i])):
                self.write_sequence(coords[i][j][0], "{0}_CDS_{1}".format(i, j), output_file, impose_strand = True)
        output_file.close()

    def write_full_CDS(self, coords, file = None, input_list = None, check_ORF = True, bare_name = False, gene_name = False, PTC_check = False):
        '''Take a set of coordinates containing the CDS of a given transcript, get the corresponding
        sequence (concatenated) and write it to file. If input_list, only do it for those CDS specified in the list.'''
        stop_codon_dict = self.get_stop_codons()
        if not file:
            file = "{0}_{1}_full_CDS.fasta".format(self.features_file_name[:-4], self.dataset)
        output_file = open(file, "w")
        if not input_list:
            input_list = list(coords.keys())
        to_write = []
        for i in input_list:
            #Ensembl considers that the stop codon is not part of the CDS, except if a stop codon is
            #split in two between two exons, in which case the 3' chunk will be labelled both as CDS and as stop.
            #These have to be filtered out.
            current_coords = [coords[i][j][0] for j in range(len(coords[i]))]
            current_stop_coords = stop_codon_dict[i]
            current_stop_coords_loc_only = [[j[2],j[3]] for j in current_stop_coords]
            current_coords = [j for j in current_coords if [j[2],j[3]] not in current_stop_coords_loc_only]
            if current_coords[0][6] == "-": 
                current_stop_coords = sort_coords(current_stop_coords, 2, reverse = True)
            else:
                current_stop_coords = sort_coords(current_stop_coords, 2, reverse = False)
            try:
                #you need a list in case the stop codon is split in half and there are therefore two stop codon features for the transcript.
                append_list = [get_sequence(j, self.genome, impose_strand = True) for j in current_stop_coords]
            except KeyError:
                append_list = []
                print(current_coords)
            append_string = "".join(append_list)
            if gene_name:
                current_name = coords[i][0][0][5]
            else:
                current_name = i
            if bare_name:
                self.write_sequence(current_coords, current_name, output_file, append_string = append_string, impose_strand = True, check_ORF = check_ORF, PTC_check = PTC_check)
            else:
                self.write_sequence(current_coords, "{0}_full_CDS".format(current_name), output_file, append_string = append_string, impose_strand = True, check_ORF = check_ORF, PTC_check = PTC_check)
        output_file.close()       
        
    def write_sequence(self, coords, name, file, impose_strand = True, check_ORF = False, append_string = None, PTC_check = False):
        '''Get the sequence at coords and write it to pre-existing file file. If impose_strand, return reverse complement if
        the feature is on the - strand. If check_ORF, verify that the returned sequence is a multiple of 3 long,
        begins with a start codon and ends with a stop codon. If append_string, append a particular string to the sequence.
        coords can be either a flat list giving the coordinates of one feature or a list of lists giving the coordinates of many features.
        In the latter case, the sequence will be concatenated.'''
        #to avoid circular import
        from conservation import check_ORF_integrity
        sequence = get_sequence(coords, self.genome, impose_strand)
        if sequence == None:
            pass
        else:
            if append_string:
                sequence = "".join([sequence, append_string])
            if check_ORF:
                ORF_check = check_ORF_integrity(sequence, PTC_check = PTC_check)
                if not ORF_check[0]:
                    pass
                else:
                    if file:
                        file.write(">{0}\n".format(name))
                        file.write("{0}\n".format(sequence))
                    else:
                        return(sequence)
            else:
                if file:
                    file.write(">{0}\n".format(name))
                    file.write("{0}\n".format(sequence))
                else:
                    return(sequence)

def bed_from_coords(coords, bed):
    with open(bed, "w") as file:
        for coord in coords:
            if coord:
                current_list = [coord[0], coord[2] - 1, coord[3], coord[4], "100", coord[6]]
                current_list = [str(i) for i in current_list]
                file.write("\t".join(current_list))
                file.write("\n")

def filter_features(file_name,input_list,start,stop,gene_pattern,transcript_pattern,gene_name_pattern, filter_trans, input_type):
    '''Read the file file_name between bytes start and stop. Keep only those features that correspond to
    transcript IDs in input_list.'''
    with open(file_name) as input_file:
        input_file.seek(start)
        lines = input_file.readlines(stop - start - 1)
        to_keep = []
        counter = 0
        if input_type == "gene":
            temp_transcript_pattern = transcript_pattern
            temp_gene_pattern = gene_pattern
            transcript_pattern = temp_gene_pattern
            gene_pattern = temp_transcript_pattern
        for i in lines:
            counter = counter + 1
            i = i.split("\t")
            text_column = i[8]
            match_obj = re.search(transcript_pattern, text_column)
            try:
                current_id = match_obj.group(0)
                if filter_trans:
                    go = current_id in input_list
                else:
                    go = True
                if go:
                    match_obj = re.search(gene_pattern, text_column)
                    try:
                        current_gene_ID = match_obj.group(0)
                        match_obj = re.search(gene_name_pattern, text_column)
                        try:
                            current_gene_name = match_obj.group(0)
                            if "\"" in current_gene_name:
                                print("Error extracting gene names!")
                                print(current_gene_name)
                                sys.exit()
                        except AttributeError:
                            current_gene_name = "no_gene_name"
                            print(text_column)
                        if input_type == "gene":
                            to_keep.append([i[0],i[2],int(i[3]),int(i[4]),current_gene_ID,current_id,i[6],current_gene_name])
                        else:
                            to_keep.append([i[0],i[2],int(i[3]),int(i[4]),current_id,current_gene_ID,i[6],current_gene_name])
                    except AttributeError:
                        pass
            except AttributeError:
                pass
    return(to_keep)

class Feature_Subset(Feature_Set):

    def __init__(self, parent_set, subset_list):
        '''Take the Feature_Set parent_set and create a similar object that only has the features corresponding to the identifiers in subset_list.'''
        Feature_Set.__init__(self, parent_set.features_file_name, parent_set.genome)
        self.names = subset_list
        self.dataset = parent_set.dataset
        self.features = [i for i in parent_set.features if i[4] in subset_list]
        self.genome = parent_set.genome
        try:
            self.families = parent_set.families
        except AttributeError:
            pass

    def create_new_features_file(self, output_file_name):
        regex = "\|".join(self.names)
        regex = regex + "\|#"
        run_process(["grep", regex], file_for_input = self.features_file_name, file_for_output = output_file_name)
        self.features_file_name = output_file_name

    def pick_random_members(self):
        '''Return a list that has all the singleton identifiers and a random member from each family.
        Because this is a Feature_Subset method, it checks that the family contains members of the subset first.'''
        picked = []
        flat_families = flatten(list(self.families.values()))
        for i in self.names:
            if i not in flat_families:
                picked.append(i)
        for i in self.families:
            remaining_members = [j for j in self.families[i] if j in self.names]
            if remaining_members:
                picked.append(random.choice(remaining_members))
        return(picked)

def add_trans_IDs_to_bed(bed_file_name, transcripts, output_file_name, two_strands = False):
    file_prefix = bed_file_name[:-4]
    features_bed_file_name = "temp_bed_file{0}.bed".format(random.random())
    unsorted_bed_file_name = "{0}_unsorted.bed".format(features_bed_file_name[:-4])
    write_features_to_bed(transcripts, unsorted_bed_file_name)
    run_process(["sort-bed", unsorted_bed_file_name], file_for_output = features_bed_file_name)
    if two_strands:
        separate_strands(features_bed_file_name, features_bed_file_name[:-4], 6)
        pos_file_name = "{0}_pos.bed".format(features_bed_file_name[:-4])
        neg_file_name = "{0}_neg.bed".format(features_bed_file_name[:-4])
    with open(bed_file_name) as input_file, open(output_file_name, "w") as output_file:
        bed_reader = csv.reader(input_file, delimiter = "\t")
        bed_writer = csv.writer(output_file, delimiter = "\t")
        for i in bed_reader:
            temp_bed_file_name = "temp_bed_file{0}.bed".format(random.random())
            rw.write_to_csv([i], temp_bed_file_name, "\t")
            if two_strands:
                if i[5] == "+":
                    bedtools_output = run_bedops(pos_file_name, temp_bed_file_name, chrom = i[0])
                else:
                    bedtools_output = run_bedops(neg_file_name, temp_bed_file_name, chrom = i[0])
            else:
                bedtools_output = run_bedops(features_bed_file_name, temp_bed_file_name, chrom = i[0])
            if bedtools_output:
                bedtools_output = bedtools_output.splitlines()
                bedtools_output = [j.split("\t") for j in bedtools_output]
                trans_IDs = [j[3] for j in bedtools_output]
                trans_IDs = list(set(trans_IDs))
                if len(trans_IDs) == 1:
                    trans_ID = trans_IDs[0]
                else:
                    trans_ID = random.choice(trans_IDs)
            else:
                trans_ID = "None"
            row_to_write = i.copy()
            try:
                row_to_write[3] = trans_ID
            except IndexError:
                row_to_write.append(trans_ID)
            bed_writer.writerow(row_to_write)
            remove_file(temp_bed_file_name)
    remove_file(features_bed_file_name)
    remove_file(unsorted_bed_file_name)
    if two_strands:
        remove_file(pos_file_name)
        remove_file(neg_file_name)

def bed_lengths(input_data):
    #test whether you have a file name or a list of bed records
    file = False
    if type(input_data) == str:
        input_data = open(input_data)
        file = True
    total_length = 0
    for line in input_data:
        if file:
            line = line.split("\t")
        total_length = total_length + (int(line[2]) - int(line[1]))
    if file:
        input_data.close()
    return(total_length)

def bed_to_CDS_indices(bed_record, CDS, check_density = False, sequence = None, motifs = None, motif_lengths = None, mismatch_warning_only = False):
    #remove phase information
    CDS = [i[0] for i in CDS]
    #convert to base 1
    hit_start = int(bed_record[1]) + 1
    hit_end = int(bed_record[2])
    #figure out which CDS region the feature maps to
    CDS_match = [i for i in range(len(CDS)) if CDS[i][2] <= hit_start and CDS[i][3] >= hit_end]
    if len(CDS_match) != 1:
        print("Bed file record and CDS coordinates don't match up!")
        print(bed_record)
        print(CDS)
        if not mismatch_warning_only:
            raise Exception
        return(["error"])
    CDS_match = CDS_match[0]
    cum_length = int(np.sum(np.array([i[3] - i[2] + 1 for i in CDS[:CDS_match]])))
    if CDS[0][6] == "+":
        relative_hit_start = cum_length + (hit_start - CDS[CDS_match][2])
    elif CDS[0][6] == "-":
        relative_hit_start = cum_length + (CDS[CDS_match][3] - hit_end)
    indices = list(range(relative_hit_start, relative_hit_start + hit_end - hit_start + 1))
    if check_density:
        test_sequence = "".join([sequence[i] for i in indices])
        test_dens = nc.get_motif_set_density(motifs, motif_lengths, test_sequence, simulants = None, concat = False, n_sim = False)
        print(test_dens[0])
        print(bed_record[-1])
        if test_dens[0] != 1.0 or len(test_dens[1]) != len(indices):
            print("Problem with motif hit positions!")
            print(bed_record)
            print(motifs)
            print(test_sequence)
            print(sequence)
            print(test_dens)
            raise Exception
    return(indices)

def bed_to_feature(bed_record, feature_type, gene_id = None, gene_description = None):
    '''
    Convert a bed record to a home-spun feature.
    '''
    if not gene_id:
        gene_id = "foo"
    if not gene_description:
        gene_description = "foo"
    feature = [bed_record[0], feature_type, int(bed_record[1]) + 1, int(bed_record[2]), bed_record[3], gene_id, bed_record[5], gene_description]
    return(feature)

def CDS_to_bed_mapping(bed_record, CDS):
    '''
    Given a bed record and (flat) CDS features, generate a hash with the mappings from CDS positions to relative positions in the bed record.
    '''
    if CDS[0][0] != bed_record[0]:
        print("Chromosomes don't match!")
        print(bed_record)
        print(CDS)
        raise Exception
    if CDS[0][6] != bed_record[5]:
        print("Strands don't match!")
        print(bed_record)
        print(CDS)
        raise Exception
    bed_range = list(range(int(bed_record[1]), int(bed_record[2])))
    exon_ranges = [list(range(i[2] - 1, i[3])) for i in CDS]
    exon_ranges = sorted(flatten(exon_ranges))
    if bed_record[5] == "-":
        exon_ranges = list(reversed(exon_ranges))
        bed_range = list(reversed(bed_range))
    output_dict = {}
    for pos, site in enumerate(exon_ranges):
        if site not in bed_range:
            output_dict[pos] = None
        else:
            output_dict[pos] = bed_range.index(site)
    return(output_dict)

def convert_coords(input_bed, output_bed, assembly1, assembly2):
    '''
    Convert the coordinates in a bed file to those from another ssembly using CrossMap.
    '''
    cwd = os.getcwd()
    user = re.search("(?<=/Users/)\w*(?=/)", cwd)
    user = user.group(0)
    os.chdir("/Users/{0}/Documents/Software/CrossMap-0.2.2/bin".format(user))
    #because the second assembly starts with a capital letter in the chain file names
    assembly2 = list(assembly2)
    assembly2[0] = assembly2[0].upper()
    assembly2 = "".join(assembly2)
    run_process(["python", "CrossMap.py", "bed", "../data/{0}To{1}.over.chain.gz".format(assembly1, assembly2), "../../../Scripts_and_data/{0}".format(input_bed), "../../../Scripts_and_data/{0}".format(output_bed)])
    os.chdir(cwd)

def convert_pos_to_bed(positions, bed_data):
    '''
    Given a feature from a bed file and contiguous relative 0-based coordinates within the feature, convert coordinates to bed format.
    '''
    result = bed_data.copy()
    if bed_data[5] == "+":
        reverse = False
    elif bed_data[5] == "-":
        reverse = True
    else:
        print("Invalid strand information!")
        print(bed_data)
        raise Exception
    if reverse:
        end = int(bed_data[2])
        start_pos = end - 1 - positions[-1]
        end_pos = end - positions[0]
    else:
        start = int(bed_data[1])
        start_pos = start + positions[0]
        end_pos = start + positions[-1] + 1
    result[1] = start_pos
    result[2] = end_pos
    return(result)

def convert_to_relative_coordinates(start_coord, end_coord, bed_data):
    '''Takes the start and end coordinates of a feature (in 1-based chromosomal coordinates)
    and bed data as a list of lists (in 0-based chromosomal coordinates) and returns a
    list of ranges that indicate where the intervals in the bed file map in the feature
    (in 0-based Python coordinates, with the feature starting at 0).'''
    ranges = []
    length = end_coord - start_coord + 1
    for i in bed_data:
        current_start = int(i[1]) - start_coord + 1
        if current_start < 0:
            current_start = 0
        current_end = int(i[2]) - start_coord + 1
        if current_end > length:
            current_end = length
        ranges.append(range(current_start, current_end))
    return(ranges)

def coords_from_bed(bed, feature_type):
    '''
    Takes a bed file and reads in the coordinates as features.
    '''
    coords = []
    with open(bed) as file:
        for line in file:
            if line:
                line = line.split("\t")
                chrom = line[0].lstrip("chr")
                coords.append([chrom, feature_type, int(line[1]) + 1, int(line[2]), line[3], ".", line[5].rstrip("\n"), "."])
    return(coords)

def fasta_from_intervals(bed_file, fasta_file, genome, force_strand = True, names = False):
    '''
    Takes a bed file and creates a fasta file with the corresponding sequences.
    If names == False, the fasta record names will be generated from the sequence coordinates.
    If names == True, the fasta name will correspond to whatever is in the 'name' field of the bed file
    '''
    genome_fasta = "Genomes/{0}/{0}.fa".format(genome)
    bedtools_args = ["bedtools", "getfasta", "-s", "-fi", genome_fasta, "-bed", bed_file, "-fo", fasta_file]
    if not force_strand:
        del bedtools_args[2]
    if names:
        bedtools_args.append("-name")
    run_process(bedtools_args)
    names, seqs = rw.read_fasta(fasta_file)
    seqs = [i.upper() for i in seqs]
    rw.write_to_fasta(names, seqs, fasta_file)

def get_GO_terms(gene_names, obo_file_name, go_file_name, pool = False):
    '''Take a list of HGNC gene names and return a set of GO annotations for each.'''
    annotations_dict = {i: [] for i in gene_names}
    obo = obo_parser.GODag(obo_file_name)

    with open(go_file_name) as go_file:
        for line in go_file:
            line = line.split("\t")
            #comment lines
            if line[0][0] != "!":
                gene_name = line[2]
                if gene_name in gene_names:
                    GO = line[4]
                    annotation = obo[GO].name
                    annotations_dict[gene_name].append(annotation)

    #uniquify the annotations for each gene
    annotations_dict = {i: list(set(annotations_dict[i])) for i in annotations_dict}

    if not pool:
        return(annotations_dict)
    else:
        annotations_list = list(annotations_dict.values())
        annotations_list = flatten(annotations_list)
        return(annotations_list)
                
def get_introns_GTF(input_file_name, output_file_name, metadata = True):
    '''Take a GTF features file and create a GTF file that has the coordinates of the introns.'''
    temp_exons_file_name = "temp_data/temp_exons_file{0}.gtf".format(random.random())
    temp_transcripts_file_name = "temp_data/temp_transcripts_file{0}.gtf".format(random.random())
    temp_exons_file_name_full = "{0}/{1}".format(os.getcwd(),temp_exons_file_name)
    temp_transcripts_file_name_full = "{0}/{1}".format(os.getcwd(),temp_transcripts_file_name)
    temp_transcripts_file = open(temp_transcripts_file_name, "w")
    temp_exons_file = open(temp_exons_file_name, "w")
    temp_transcripts_writer = csv.writer(temp_transcripts_file, delimiter = "\t")
    temp_exons_writer = csv.writer(temp_exons_file, delimiter = "\t")
    input_file = open(input_file_name)
    input_reader = csv.reader(input_file, delimiter = "\t")
    meta_data = []
    for i in range(5):
        meta_data.append(next(input_reader)[0])
    for i in input_reader:
        try:
            if i[2] == "exon":
                temp_exons_writer.writerow(i)
            elif i[2] == "transcript":
                temp_transcripts_writer.writerow(i)
        except IndexError:
            print(i)
            pass
    input_file.close()
    temp_exons_file.close()
    temp_transcripts_file.close()
    #the "-s" flag means that the exon and the transcript have to be on the same strand.
    bedtools_output = run_process(["/usr/local/bin/subtractBed", "-a", temp_transcripts_file_name_full, "-b", temp_exons_file_name_full, "-s"])
    #to put back the metadata that bedtools gets rid of
    meta_data = "\n".join(meta_data)
    if metadata:
        bedtools_output = meta_data + "\n" + bedtools_output
    sed_output = run_process(["sed", "s/transcript/intron/"], input_to_pipe = bedtools_output)
    sed_output2 = run_process(["sed", "s/\"(\")/\1/"], input_to_pipe = bedtools_output)
    rw.write_all(sed_output, output_file_name)
    output_data = rw.read_many_fields(output_file_name, delimiter = "\t")
    for i in range(5, len(output_data)):
        output_data[i][3] = int(output_data[i][3]) + 1
    rw.write_to_csv(output_data, output_file_name, delimiter = "\t")
    os.remove(temp_exons_file_name)
    os.remove(temp_transcripts_file_name)

def get_sequence(coords, genome, impose_strand = False, bed_input = False):
    '''Get the sequence corresponding to the 1-based coords of a feature from a genome sequence.
    If a list of features is given, the resulting sequences will be concatenated in the order in which they appear
    in the list (reversed order if impose_strand is True and the feature is on the - strand).
    OPTIONS
    impose_strand: if True, the reverse complement of the genome sequence will be returned for features on the - strand.
    False by default.'''
    file_name = "Genomes/{0}/{0}.fa".format(genome)
    genome_seq = Fasta("Genomes/{0}/{0}.fa".format(genome))
    if bed_input:
        strand_col = 5
    else:
        strand_col = 6
    #check whether only a single feature or a list of features was supplied. Convert to list if it's the former.
    if type(coords[0]) == str:
        coords = [coords]
    complete_sequence = []
    if impose_strand and coords[0][strand_col] == "-":
        coords = list(reversed(coords))
    for i in coords:
        chr_name = i[0]
        #pyfaidx works with 0-based coordinates
        try:
            if bed_input:
                sequence = genome_seq[chr_name.lstrip("chr")][(int(i[1])):(int(i[2]))]
            else:
                sequence = genome_seq[chr_name][(i[2]-1):(i[3])]
        except KeyError:
            print("{0} not in genome fasta!".format(chr_name))
            return(None)
        complete_sequence.append(str(sequence).upper())
    sequence = "".join(complete_sequence)
    genome_seq.close()
    if impose_strand:
        if coords[0][strand_col] == "-":
            sequence = Seq(sequence, IUPAC.unambiguous_dna)
            sequence = sequence.reverse_complement()
            sequence = str(sequence)
    return(sequence)

def intersect_bed(input_list, bed_file, overlap = False, write_both = False, sort = False,
                  force_strand = False, modify_chr_ids = False, no_name_check = False, no_dups = True, chrom = None, use_bedops = False, bed_input = False):
    '''Use bedtools/bedops to intersect coordinates from a 1-based features list with a 0-based bedfile.
    Return those features in the input list that overlap with coordinates in the bed file.
    OPTIONS
    overlap: minimum overlap required as a fraction of the features in the input list
    write_both: if True, return also the lines in the bed file that overlap the features in the input list.
    False by default.
    no_name_check: if set to False, chromosome names aren't checked for consistency between files. Only set
    to True when you're sure you know what you're doing.
    no_dups: keeps those entries in the input list that have any overlaps to elements in the bed file
    without duplicating input_list features in the case of multiple overlaps.
    see intersect_bed_return_bed for other options.
    bed_input: input_list is a list of bed-records (base 0 ) rather than features (base 1)'''
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    temp_file_name_full = "{0}/{1}".format(os.getcwd(),temp_file_name)
    bed_file_name_full = "{0}/{1}".format(os.getcwd(),bed_file)
    if not bed_input:
        if modify_chr_ids:
            input_list = [[modify_chr(i[0]), i[2] - 1, i[3], i[4], ".", i[6]] for i in input_list]
        else:
            input_list = [[i[0], i[2] - 1, i[3], i[4], ".", i[6]] for i in input_list]
    rw.write_to_csv(input_list, temp_file_name, delimiter = "\t")
    if use_bedops:
        bedtools_output = run_bedops(temp_file_name_full, bed_file_name_full, force_strand, write_both, chrom, overlap, sort)
    else:
        bedtools_output = run_bedtools(temp_file_name_full, bed_file_name_full, force_strand, write_both, chrom, overlap, sort, no_name_check, no_dups)
    os.remove(temp_file_name)
    bedtools_output = bedtools_output.splitlines()
    bedtools_output = [i.split("\t") for i in bedtools_output]
    if write_both:
        bedtools_output = [[i[0], int(i[1]) + 1, int(i[2]), i[3], i[5], i[6:]] for i in bedtools_output]
    else:
        bedtools_output = [[i[0], int(i[1]) + 1, int(i[2]), i[3], i[5]] for i in bedtools_output]
    return(bedtools_output)

def intersect_bed_return_bed(bed_file, input_list, use_bedops = False, overlap = False, write_both = False, sort = False, output_file = None,
                             force_strand = False, modify_chr_ids = False, no_name_check = False, no_dups = True, chrom = None, bed_input = False, intersect = False):
    '''Use bedtools/bedops to intersect coordinates from a 1-based features list with a 0-based bedfile.
    Return those lines in the bed file that overlap with features in the features list.
    OPTIONS
    overlap: minimum overlap required as a fraction of the clusters in the bed file.
    write_both: if True, return also the features that overlap the clusters in the bed file (only
    valid when using bedtools).
    False by default.
    force_strand: check that the feature and the bed interval are on the same strand (only valid with bedtools)
    no_name_check: set to False so that chromosome names wouldn't be checked for consistency (only valid with bedtools)
    chrom: search only on a specific chromosome (only valid with bedops).
    see intersect_bed for other options'''
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    temp_file_name_full = "{0}/{1}".format(os.getcwd(), temp_file_name)
    bed_file_name_full = "{0}/{1}".format(os.getcwd(), bed_file)
    if not bed_input:
        if modify_chr_ids:
            input_list = [[modify_chr(i[0]), i[2] - 1, i[3], i[4], ".", i[6]] for i in input_list]
        else:
            input_list = [[i[0], i[2] - 1, i[3], i[4], ".", i[6]] for i in input_list]
    rw.write_to_csv(input_list, temp_file_name, delimiter = "\t")
    if output_file:
        temp_file_name2 = "temp_data/IBRB_temp_bed_file{0}.bed".format(random.random())
        temp_file_name_full2 = "{0}/{1}".format(os.getcwd(), temp_file_name2)
    else:
        temp_file_name_full2 = None
    if use_bedops:
        bedtools_output = run_bedops(bed_file_name_full, temp_file_name_full, force_strand, write_both, chrom, overlap, sort, output_file = temp_file_name_full2, intersect = intersect)
    else:
        bedtools_output = run_bedtools(bed_file_name_full, temp_file_name_full, force_strand, write_both, chrom, overlap, sort, no_name_check, no_dups, output_file = temp_file_name_full2, intersect = intersect)
    os.remove(temp_file_name)
    if output_file:
        run_process(["mv", temp_file_name2, output_file])
        pass
    else:
        bedtools_output = bedtools_output.splitlines()
        bedtools_output = [i.split("\t") for i in bedtools_output]
    return(bedtools_output)

def map_relative_to_feature(dictionary, CDSs):
    mapped_dictionary = {i: {} for i in dictionary}
    for trans in dictionary:
        if dictionary[trans]:
            dictionary[trans] = sorted(dictionary[trans])
            cumul_length = 0
            if CDSs[trans][0][0][6] == "+":
                mult_factor = 1
                start_coord = 2
            else:
                mult_factor = -1
                start_coord = 3
            CDS_pos = 0
            for position in dictionary[trans]:
                found = False
                while not found:
                    try:
                        current_CDS = CDSs[trans][CDS_pos][0]
                        current_length = current_CDS[3] - current_CDS[2] + 1
                        if (current_length + cumul_length) > position:
                            mapped_dictionary[trans][position] = (position * mult_factor) - (cumul_length * mult_factor) + current_CDS[start_coord]
                            found = True
                        else:
                            cumul_length = cumul_length + current_length
                            CDS_pos = CDS_pos + 1
                    except IndexError:
                        #this is when the position maps to a stop
                        diff = position - cumul_length
                        if diff < 3:
                            found = True
                        else:
                            print("CDS and position coordinates don't match up!")
                            print(CDSs[trans])
                            print(cumul_length)
                            print(CDS_pos)
                            raise IndexError
    return(mapped_dictionary)

def matches_to_bed(bed_record, matches):
    first_cluster_start = [matches[0]]
    last_cluster_end = [matches[-1]]
    cluster_starts = [j for i, j in enumerate(matches[1:]) if j - matches[i] != 1]
    first_cluster_start.extend(cluster_starts)
    cluster_starts = first_cluster_start
    cluster_ends = [j for i, j in enumerate(matches[:-1]) if matches[i+1] - j != 1]
    cluster_ends.extend(last_cluster_end)
    strand = bed_record[5]
    if strand == "+":
        feature_start = int(bed_record[1])
        cluster_ranges = [(cluster_starts[i] + feature_start, cluster_ends[i] + feature_start + 1) for i in range(len(cluster_starts))]
    else:
        feature_end = int(bed_record[2])
        cluster_ranges = [(feature_end - cluster_ends[i] - 1, feature_end - cluster_starts[i]) for i in range(len(cluster_starts))]
    chrom = bed_record[0]
    name = bed_record[3]
    new_records = [[chrom, int(i[0]), int(i[1]), name, 100, strand] for i in cluster_ranges]
    return(new_records)

def merge_intervals(input_file, output_file, force_strand = True):
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    bedtools_args = ["mergeBed", "-s", "-c", "4,5,6", "-o", "distinct,distinct,distinct", "-i", input_file]
    if not force_strand:
        del bedtools_args[1]
    run_process(bedtools_args, file_for_output = temp_file_name)
    bed_data = rw.read_many_fields(temp_file_name, "\t")
    bed_data = [[i[0], i[1], i[2], i[4], i[5], i[3]] for i in bed_data]
    rw.write_to_csv(bed_data, output_file, "\t")
    os.remove(temp_file_name)

def modify_chr(chr_id):
    '''Take a bare chromosome name and prefix "chr" to it, except if it's one of
    the weird ones.'''
    if chr_id[0] not in ["C", "c", "J", "G"]:
        chr_id = "chr{0}".format(chr_id)
    return(chr_id)

def region_indices_to_full_indices(region_indices, location_indices):
    start_pos = location_indices[0]
    full_indices = [i + start_pos for i in region_indices]
    return(full_indices)

def run_bedops(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, output_file = None, intersect = False):
    sorting_temp_file_name = "{0}/temp_data/temp_sort_{1}.bed".format(os.getcwd(), random.random())
    sorting_temp_file_name2 = "{0}/temp_data/temp_sort2_{1}.bed".format(os.getcwd(), random.random())
    if intersect:
        command = "--intersect"
    else:
        command = "--element-of"
    if sort:
        run_process(["sort-bed", A_file], file_for_output = sorting_temp_file_name)
        run_process(["sort-bed", B_file], file_for_output = sorting_temp_file_name2)
        bedops_args = ["bedops", "--chrom", "foo", command, "1", sorting_temp_file_name, sorting_temp_file_name2]
    else:
        bedops_args = ["bedops", "--chrom", "foo", command, "1", A_file, B_file]
    if overlap:
        bedops_args[4] = overlap
    if chrom:
        bedops_args[2] = chrom
    else:
        del bedops_args[1:3]
    if intersect:
        del bedops_args[4]
    if force_strand:
        print("Bedops can't search by strand! Either use bedtools or divide input data by strand!")
        raise Exception
    if write_both:
        print("Bedops can't write both features!")
        raise Exception
    bedops_output = run_process(bedops_args, file_for_output = output_file)
    remove_file(sorting_temp_file_name)
    remove_file(sorting_temp_file_name2)
    return(bedops_output)

def run_bedtools(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, no_name_check = False, no_dups = True, hit_number = False, output_file = None, intersect = False):
    if write_both:
        write_option = "-wo"
    elif hit_number:
        write_option = "-c"
    else:
        write_option = "-wa"
    bedtools_args = ["intersectBed", "-a", A_file,"-b", B_file, write_option]
    if overlap:
        bedtools_args.extend(["-f", str(overlap)])
    if force_strand:
        bedtools_args.append("-s")
    if no_name_check:
        bedtools_args.append("-nonamecheck")
    if no_dups:
        bedtools_args.append("-u")
    if chrom:
        print("Bedtools cannot be restricted to a single chromosome. Use bedops!")
        raise Exception
    if sort:
        print("Bedtools is not set up to sort your bed file!")
        raise Exception
    if intersect:
        print("Bedtools has not been set up to take intersections. Either implement it or use bedops!")
    bedtools_output = run_process(bedtools_args, file_for_output = output_file)
    return(bedtools_output)

def separate_strands(input_file, output_file_prefix, strand_column):
    pos_file_name = "{0}_pos.bed".format(output_file_prefix)
    run_process(["awk", "${0} == \"+\"".format(strand_column)], file_for_input = input_file, file_for_output = pos_file_name)
    neg_file_name = "{0}_neg.bed".format(output_file_prefix)
    run_process(["awk", "${0} == \"-\"".format(strand_column)], file_for_input = input_file, file_for_output = neg_file_name)

def setup(features_file, genome, dataset, families_file = None):
    '''
    Set up a Feature_Set object and attribute families, if required.
    '''
    fs = Feature_Set(features_file, genome)
    fs.set_dataset(dataset)
    if families_file:
        families = rw.read_families(families_file)
        fs.add_families(families)
    return(fs)

def sort_coords(coords_to_sort, index, reverse = False):
    '''Sort a list of features based on the coordinate in column 'index'.
    OPTIONS
    reverse: if True, sort the list in descending order. False by default.'''
    if reverse:
        sorted_coords = sorted(coords_to_sort, key = lambda x:x[index], reverse = True)
    else:
        sorted_coords = sorted(coords_to_sort, key = lambda x:x[index])
    return(sorted_coords)

def sort_bed(input_file_name, output_file_name):
    temp_file_name = "temp_data/temp_sorted_bed.bed"
    run_process(["sort-bed", input_file_name], file_for_output = temp_file_name)
    run_process(["mv", temp_file_name, output_file_name])

def trim_sequence_coords(coords, phase):
    if coords[6] == "+":
        strand = "pos"
    elif coords[6] == "-":
        strand = "neg"
    if strand == "pos":
        if phase == 0:
            pass
        elif phase == 1:
            coords[2] = coords[2] + 2
        elif phase == 2:
            coords[2] = coords[2] + 1
        else:
            print("Invalid phase information!")
            sys.exit()
        seq_length = coords[3] - coords[2] + 1
        if seq_length%3 == 0:
            pass
        elif seq_length%3 == 1:
            coords[3] = coords[3] - 1
        else:
            coords[3] = coords[3] - 2
    elif strand == "neg":
        if phase == 0:
            pass
        elif phase == 1:
            coords[3] = coords[3] - 2
        elif phase == 2:
            coords[3] = coords[3] - 1
        else:
            print("Invalid phase information!")
            sys.exit()
        seq_length = coords[3] - coords[2] + 1
        if seq_length%3 == 0:
            pass
        elif seq_length%3 == 1:
            coords[2] = coords[2] + 1
        else:
            coords[2] = coords[2] + 2
    return(coords)

def uniquify_lines(input_file, output_file, sort = False):
    '''
    Uniquify the lines in a file.
    '''
    with open(input_file) as file:
        input_data = file.readlines()
    input_data = list(set(input_data))
    temp_file_name = output_file
    if sort:
        temp_file_name = "temp_data/temp_uniquify.bed"
    with open(temp_file_name, "w") as file:
        for line in input_data:
            file.write(line)
    if sort:
        sort_bed(temp_file_name, output_file)
        os.remove(temp_file_name)
        
def write_features_to_bed(features, bed_file_name, bare = False, modify_chr_ids = False):
    '''
    Take either a feature dictionary or a list of features
    and make a bed file from the coordinates.
    '''
    if type(features) == dict:
        features = list(features.values())
    if type(features[0][0]) == list:
        features = flatten(features)
    if not bare:
        if modify_chr_ids:
            features = [[modify_chr(i[0]), i[2] - 1, i[3], i[4], ".", i[6]] for i in features]
        else:
            features = [[i[0], i[2] - 1, i[3], i[4], ".", i[6]] for i in features]
    else:
        if modify_chr_ids:
            features = [[modify_chr(i[0]), i[2] - 1, i[3], i[4]] for i in features]
        else:
            features = [[i[0], i[2] - 1, i[3], i[4]] for i in features]
    rw.write_to_csv(features, bed_file_name, delimiter = "\t")





        
        

        
