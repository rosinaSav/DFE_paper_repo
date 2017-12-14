from housekeeping import parse_arguments, remove_file, run_process
import read_and_write as rw

def mDFEest(model, input_file, n_spikes = None, repetitions = None, fold_SFS = True, pop_change = False, seed = None):
    flags = []

    if fold_SFS:
        fold_SFS = 1
    else:
        fold_SFS = 0
    #this looks weird but is normal: this value will be the value of conpop in the multiDFE call, meaning it'll be 1 with constant population size
    if pop_change:
        pop_change = 0
    else:
        pop_change = 1

    if model == "lognormal":
        model_code = 4
        par_number = 2
    elif model == "gamma":
        model_code = 2
        par_number = 2
    elif model == "beta":
        model_code = 3
        par_number = 2
    elif model == "spikes":
        model_code = 0
        if not n_spikes:
            print("To be able to use a spikes model, you need to specify the number of spikes.")
            raise Exception
        par_number = (2 * n_spikes) - 1
        flags = ["-ranrep", repetitions, "-nspikes", n_spikes]
    elif model == "steps":
        model_code = 1
        if not n_spikes:
            print("To be able to use a steps model, you need to specify the number of steps.")
            raise Exception
        par_number = (2 * n_spikes) - 1
        flags = ["-ranrep", repetitions, "-nspikes", n_spikes]
    elif model == "six_spikes":
        model_code = 5
        par_number = 5
        flags = ["-ranrep", repetitions]
    else:
        print("{0} is not a valid model name!".format(model))
        raise Exception

    input_file_short = input_file.split("/")
    input_file_short = input_file_short[-1]

    if not os.path.exists("../multidfe/{0}".format(input_file_short)):
        run_process(["cp", input_file, "../multidfe"])
    MDE_output = "{0}.MAXL.out".format(input_file_short)
    current_dir = os.getcwd()
    os.chdir("../multidfe")
    arguments = ["./MultiDFE", "-N1", 100, "-conpop", pop_change, "-sfsfold", fold_SFS, "-selmode", model_code, "-file", input_file_short]
    if seed:
        seed_string = "GSL_RNG_SEED={0}".format(seed)
        arguments = [seed_string] + arguments
    arguments.extend(flags)
    print(" ".join([str(i) for i in arguments]))
    run_process(arguments)
    output = rw.read_many_fields(MDE_output, "\t")[0]
    output = [i.split(":") for i in output if ":" in i]
    output = {i[0]: float(i[1]) for i in output}
    ll = output["L"]
    print("\n")
    print(par_number)
    print(ll)
    AIC = (2 * par_number) - (2 * ll)
    output["AIC"] = AIC
    if n_spikes:
        output["model"] = "{0}_{1}".format(model, n_spikes)
    else:
        output["model"] = model
    remove_file(MDE_output)
    os.chdir(current_dir)
    return(output)

def write_mDFEest_output(output, file, change_mode):
    file.write("{0}\t".format(output["model"]))
    file.write("{0}\t".format(str(change_mode)))
    file.write("{0}\t".format(output["AIC"]))
    file.write("{0}\t".format(output["Nes_0.0_0.1"]))
    file.write("{0}\t".format(output["Nes_0.1_1.0"]))
    file.write("{0}\t".format(output["Nes_1.0_10.0"]))
    file.write("{0}\t".format(output["Nes_10.0_100.0"]))
    file.write(str(output))
    file.write("\n") 

def main():
    description = "Run mDFEest."
    args = parse_arguments(description, ["hit_file", "control_file", "SNP_file", "SNP_number", "input_file", "output_file", "seed", "fixed_model", "new_input", "shuffle", "fix_pop_change"], ints = [3], flags = [8, 9, 10])
    hit_file, control_file, SNP_file, SNP_number, input_file, output_file, seed, fixed_model, new_input, shuffle, fix_pop_change = args.hit_file, args.control_file, args.SNP_file, args.SNP_number, args.input_file, args.output_file, args.seed, args.fixed_model, args.new_input, args.shuffle, args.fix_pop_change

    if new_input:
        remove_file("../multidfe/{0}".format(input_file.split("/")[-1]))
        arguments = ["python3", "mDFEest_input.py", hit_file, control_file, SNP_file, SNP_number, input_file]
        if shuffle:
            arguments.append("--shuffle")
        run_process(arguments)
    
    if seed == "None":
        seed = None
    else:
        seed = float(seed)

    if fix_pop_change:
        pop_change = [True]
    else:
        pop_change = [False, True]

    if fixed_model == "None":
        allowed = ["lognormal", "gamma", "beta", "spikes", "steps", "fixed six spikes"]
        spike_range = [2, 6]
    else:
        allowed = [fixed_model]
        spike_range = [2, 3]

    with open(output_file, "w") as file:
        file.write("model\tpop_change\tAIC\tNes_0.0_0.1\tNes_0.1_1.0\tNes_1.0_10.0\tNes_10.0_100.0\traw\n")
        for change_mode in pop_change:
    
            print("\nPopulation expansion: {0}.".format(str(change_mode)))

            if "lognormal" in allowed:
                print("lognormal model:")
                output = mDFEest("lognormal", input_file, pop_change = change_mode, seed = seed)
                print(output)
                write_mDFEest_output(output, file, change_mode)

            if "gamma" in allowed:
                print("gamma model:")
                output = mDFEest("gamma", input_file, pop_change = change_mode, seed = seed)
                print(output)
                write_mDFEest_output(output, file, change_mode)

            if "beta" in allowed:
                print("beta model:")
                output = mDFEest("beta", input_file, pop_change = change_mode, seed = seed)
                print(output)
                write_mDFEest_output(output, file, change_mode)

            for spike_number in range(spike_range[0], spike_range[1]):

                if "spikes" in allowed:
                    print("{0}-spikes model:".format(spike_number))
                    output = mDFEest("spikes", input_file, n_spikes = spike_number, seed = seed, repetitions = 10, pop_change = change_mode)
                    print(output)
                    write_mDFEest_output(output, file, change_mode)

                if "steps" in allowed:
                    print("{0}-steps model:".format(spike_number))
                    output = mDFEest("steps", input_file, n_spikes = spike_number, seed = seed, repetitions = 10, pop_change = change_mode)
                    print(output)
                    write_mDFEest_output(output, file, change_mode)

            if "fixed six spikes" in allowed:
                print("fixed six spikes model:")
                output = mDFEest("six_spikes", input_file, pop_change = change_mode, seed = seed)
                print(output)
                write_mDFEest_output(output, file, change_mode)

if __name__ == "__main__":
    main()

    

