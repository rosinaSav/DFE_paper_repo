from bedtools_games import setup
import conservation
from housekeeping import parse_arguments
import re
import read_and_write as rw

def convert_region_to_CDS_func(pos_dict, bed, CDSs):
    fasta = re.sub("bed", "fasta", bed)
    mapping_dict = conservation.map_regions_to_CDS(fasta, bed, None, None, CDSs, trans_ids = True)
    out_pos = {}
    for region in pos_dict:
        idn = mapping_dict[region]["idn"]
        if idn not in out_pos:
            out_pos[idn] = []
        start = mapping_dict[region]["flank indices in CDS"][0]
        out_pos[idn].extend([i + start for i in pos_dict[region]])
    out_pos = {i: sorted(out_pos[i]) for i in out_pos}
    return(out_pos)

def main():
    description = "Take a positions file with hits within exonic subregions and convert them to full CDS indices."
    args = parse_arguments(description, ["positions_file", "bed_file", "genome", "features_file", "dataset", "output_file"])
    positions_file, bed_file, genome, features_file, dataset, output_file = args.positions_file, args.bed_file, args.genome, args.features_file, args.dataset, args.output_file

    fs = setup(features_file, genome, dataset)
    CDSs = fs.get_CDS()

    pos_dict = rw.read_pos(positions_file)

    converted_pos = convert_region_to_CDS_func(pos_dict, bed_file, CDSs)

    rw.write_pos(converted_pos, output_file)

if __name__ == "__main__":
    main()
