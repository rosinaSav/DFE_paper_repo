'''
Author: Rosina Savisaar.
Run mDFEest with shuffled input to check the false positive rate.
'''

from CpG_sim_Keightley import mDFEest
from housekeeping import parse_arguments, remove_file, run_process
import random

def main():
    description = "Run mDFEest with shuffled input to check the false positive rate."
    args = parse_arguments(description, ["hits_file", "controls_file", "output_file", "n_sim", "SNP_file", "SNP_number", "hit_reduce", "control_reduce", "const_pop"], ints = [3, 5], floats = [6, 7], flags = [8])
    hits_file, controls_file, output_file, n_sim, SNP_file, SNP_number, hit_reduce, control_reduce, const_pop = args.hits_file, args.controls_file, args.output_file, args.n_sim, args.SNP_file, args.SNP_number, args.hit_reduce, args.control_reduce, args.const_pop

    with open(output_file, "w") as file:
        for sim in range(n_sim):
            print(sim)

            temp_hits_file = "temp_data/hits_file{0}.txt".format(random.random())
            temp_controls_file = "temp_data/controls_file{0}.txt".format(random.random())
            temp_input_file = "temp_data/input_file{0}.txt".format(random.random())

            #shuffle hits and controls for negative control
            run_process(["python3", "shuffle_hits_and_controls.py", hits_file, controls_file, temp_hits_file, temp_controls_file, hit_reduce, control_reduce])

            #generate multiDFEest input file
            run_process(["python3", "mDFEest_input.py", temp_hits_file, temp_controls_file, SNP_file, SNP_number, temp_input_file])

            output = mDFEest("beta", temp_input_file, pop_change = True)

            print(output)
            print(output["Nes_0.0_0.1"])
            print(output["Nes_0.1_1.0"])

            file.write("{0}\t{1}\t{2}".format(sim, output["Nes_0.0_0.1"], output["Nes_0.1_1.0"]))

            #if you also want to run with fixed population size
            if const_pop:
                output = mDFEest("beta", temp_input_file, pop_change = False)

                file.write("{0}\t{1}\t{2}".format(sim, output["Nes_0.0_0.1"], output["Nes_0.1_1.0"]))

            file.write("\n")

            remove_file(temp_hits_file)
            remove_file(temp_controls_file)
            remove_file(temp_input_file)

if __name__ == "__main__":
    main()
