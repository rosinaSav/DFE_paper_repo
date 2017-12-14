from housekeeping import parse_arguments
from INSIGHT import reduce_dict
import numpy as np
import read_and_write as rw

def main():
    description = "Take a hits file and a control file and shuffle which elements are in which."
    args = parse_arguments(description, ["input_hits", "input_controls", "output_hits", "output_controls", "hit_reduce", "control_reduce"], floats = [4, 5])
    input_hits, input_controls, output_hits, output_controls, hit_reduce, control_reduce = args.input_hits, args.input_controls, args.output_hits, args.output_controls, args.hit_reduce, args.control_reduce

    hits = rw.read_pos(input_hits)
    controls = rw.read_pos(input_controls)

    if hit_reduce > 0:
        hits = reduce_dict(hits, hit_reduce)
        controls = reduce_dict(controls, control_reduce)
        rw.write_pos(hits, output_hits)
        rw.write_pos(controls, output_controls)
    else:
        with open(output_hits, "w") as hits_o, open(output_controls, "w") as controls_o:
            for gene in hits:
                hit_length = len(hits[gene])
                combined = hits[gene] + controls[gene]
                current_hits_o = sorted(np.random.choice(combined, size = hit_length, replace = False))
                current_controls_o = sorted([i for i in combined if i not in current_hits_o])
                hits_o.write("{0}\t{1}\n".format(gene, ",".join([str(i) for i in current_hits_o])))
                controls_o.write("{0}\t{1}\n".format(gene, ",".join([str(i) for i in current_controls_o])))
    
if __name__ == "__main__":
    main()
