'''
Author: Rosina Savisaar.
Take a positions file with hits within exonic subregions and convert them to full CDS indices.
'''

from bedtools_games import setup
import conservation
from housekeeping import parse_arguments
import re
import read_and_write as rw

def convert_region_to_CDS_func(pos_dict, bed, CDSs):
    '''
    Given a bed file containing the coordinates of core/flank regions and a dictionary of hit positions
    (in the subregion),
    as well as a list of list of correpsonding full CDS coordinates,
    determine where the cores/flanks map in the full CDSs.
    '''
    #assumes there is a fasta file that has the same name as the bed file,
    #except for the extension
    fasta = re.sub("bed", "fasta", bed)
    #this gives you where the flanks/cores map in the full CDS
    #with the region (core/flank) identifiers as keys
    mapping_dict = conservation.map_regions_to_CDS(fasta, bed, None, None, CDSs, trans_ids = True)
    out_pos = {}
    for region in pos_dict:
        #in the output dictionary you want transcript identifiers as keys
        idn = mapping_dict[region]["idn"]
        if idn not in out_pos:
            out_pos[idn] = []
        #you determine where the subregion starts in the full sequence
        start = mapping_dict[region]["flank indices in CDS"][0]
        #and then based on that, convert the positions in the subregion to positions in the
        #full CDS
        out_pos[idn].extend([i + start for i in pos_dict[region]])
    out_pos = {i: sorted(out_pos[i]) for i in out_pos}
    return(out_pos)

def main():
    description = "Take a positions file with hits within exonic subregions and convert them to full CDS indices."
    args = parse_arguments(description, ["positions_file", "bed_file", "genome", "features_file", "dataset", "output_file"])
    positions_file, bed_file, genome, features_file, dataset, output_file = args.positions_file, args.bed_file, args.genome, args.features_file, args.dataset, args.output_file

    #set up data
    fs = setup(features_file, genome, dataset)
    CDSs = fs.get_CDS()

    pos_dict = rw.read_pos(positions_file)

    #do actual work
    converted_pos = convert_region_to_CDS_func(pos_dict, bed_file, CDSs)

    #write output to file
    rw.write_pos(converted_pos, output_file)

if __name__ == "__main__":
    main()
