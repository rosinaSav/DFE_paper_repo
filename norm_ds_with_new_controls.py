from bedtools_games import Feature_Set
import conservation
from housekeeping import flatten, parse_arguments, remove_file, run_in_parallel, run_process, update_counter
from INSIGHT import parse_basinhoppin_pos, shuffle_dictionaries
import my_stats as ms
import nucleotide_comp as nc
import numpy as np
import random
import read_and_write as rw

def CpG_frequency(fasta, hits, controls):
    '''
    Compare the CpG frequency at hit vs control sites.
    '''
    names, seqs = rw.read_fasta(fasta)
    seqs = {names[i]: seqs[i] for i in range(len(names))}
    hit_site_counter = 0
    hit_CpG_counter = 0
    control_site_counter = 0
    control_CpG_counter = 0
    for name in hits:
        seq = seqs[name]
        current_true_dints = [seq[i - 1: i + 1] for i in hits[name] if i != 0] + [seq[i: i + 2] for i in hits[name] if i != (len(seq) - 1)]
        current_control_dints = [seq[i - 1: i + 1] for i in controls[name] if i != 0] + [seq[i: i + 2] for i in controls[name] if i != (len(seq) - 1)]
        hit_site_counter = hit_site_counter + len(current_true_dints)
        control_site_counter = control_site_counter + len(current_control_dints)
        hit_CpG_counter = hit_CpG_counter + len([i for i in current_true_dints if i == "CG" or i == "GC"])
        control_CpG_counter = control_CpG_counter + len([i for i in current_control_dints if i == "CG" or i == "GC"])
    hit_freq = hit_CpG_counter/hit_site_counter
    control_freq = control_CpG_counter/control_site_counter
    print("Hit CpG frequency: {0}.".format(hit_freq))
    print("Control CpG frequency: {0}.".format(control_freq))
    return(hit_freq, control_freq)

def get_control_sites(fasta, genome, feature_set, families_file, dataset, temp_motifs_file, hit_file, control_file, error_file, anc_CG_file, high_CG_file, flags):
    '''
    Given motifs and sequences, pick control sites using the optimization method.
    '''
    arguments = ["python3", "pick_control_sites.py", fasta, genome, feature_set, families_file, dataset,
                 temp_motifs_file, 10, hit_file, 500, 10, control_file, error_file,
                 "None", anc_CG_file, high_CG_file, "--old_motif_format"]

    arguments = arguments + flags
    run_process(arguments)

def get_density(fasta, motifs, fs):
    '''
    Determine motif density.
    '''
    densities_output_file_name = "temp_data/temp{0}.txt".format(random.random())
    sim_densities_output_file_name = "temp_data/temp{0}.txt".format(random.random())
    positions_output_file_name = "temp_data/temp{0}.txt".format(random.random())
    sim_positions_output_folder_name = "temp_data/temp{0}".format(random.random()) 
    density_output = nc.get_sequence_set_density(fasta, None, motifs, None, None, densities_output_file_name, sim_densities_output_file_name, positions_output_file_name, sim_positions_output_folder_name, feature_set = fs, concat = False, positions = False)
    #median density across CDSs
    print("Median density: {0}.".format(density_output["median density"]))
    density_output = nc.get_sequence_set_density(fasta, None, motifs, None, None, densities_output_file_name, sim_densities_output_file_name, positions_output_file_name, sim_positions_output_folder_name, feature_set = fs, concat = True, positions = False)
    #overall density across all CDSs
    print("Overall density: {0}.".format(density_output["density"]))

def get_new_method_results(hit_file, control_file, hit_phylip, control_phylip, correspondances, alignments, fasta, baseml = False, return_CpG = False, global_fasta = None, return_overall = False, motifs = None, fs = None, regions = False):
    '''
    Calculate normalized dS.
    '''
    if "_degen.txt" in hit_file:
        degen_hits_file = hit_file
        degen_controls_file = control_file
        hit_file = hit_file[:-10]
        control_file = control_file[:-10]
    else:
        degen_hits_file = None
        degen_controls_file = None
    hits = parse_basinhoppin_pos(hit_file)
    controls = parse_basinhoppin_pos(control_file)

    try:
        #write control and hit sequences to PHYLIP files
        conservation.write_hits_to_phylip(fasta, hits, hit_phylip, correspondances, alignments, degen_hits_file, baseml = baseml, fs = fs, regions = regions, global_fasta = global_fasta)
        conservation.write_hits_to_phylip(fasta, controls, control_phylip, correspondances, alignments, degen_controls_file, baseml = baseml, fs = fs, regions = regions, global_fasta = global_fasta)

        #if you're doing nucleotide-based rather than codon-based
        if baseml:
            method = "baseml"
            statistic = "tree length"
        else:
            method = "gy"
            statistic = "dS"

        #if you want to return the density * normalized dS statistic, you need the density
        if return_overall:
            density = nc.get_sequence_set_density(fasta, None, motifs, None, False, "temp_data/temp_dens1.txt", "temp_data/temp_dens2.txt", "temp_data/temp_pos.txt", None, feature_set = fs, concat = True, positions = False)["density"]
            print("Density: {0}.".format(density))

        #get dS estimates from PAML
        hit_ds = conservation.run_codeml(hit_phylip, "temp_data/temp_{0}.phy".format(random.random()), method = method)[statistic]
        control_ds = conservation.run_codeml(control_phylip, "temp_data/temp_{0}.phy".format(random.random()), method = method)[statistic]

        remove_file(control_phylip)

        hit_freq, control_freq = CpG_frequency(fasta, hits, controls)
        print("Hit dS: {0}.".format(hit_ds))
        print("Control dS: {0}.".format(control_ds))
        norm_ds = (hit_ds - control_ds)/control_ds
        print("Normalized dS: {0}.\n".format(norm_ds))

        if return_overall:
            overall = norm_ds * density
            print("Overall decrease: {0}.\n".format(overall))
            return(norm_ds, density, overall)
              
        if return_CpG:
            return(norm_ds, hit_freq, control_freq)
        return((hit_ds - control_ds)/control_ds)
    except conservation.NoDataException:
        print("No input sequence available.")
        if return_CpG:
            return(None, None, None)
        return(None)

def get_sim_p(norm_ds, hit_file, control_file, correspondances, alignments, fasta, n_sim, sim_ds_file = None, baseml = False, reverse_site_numbers = False):
    '''
    Get an empirical p-value for the normalized dS estimate.
    '''
    if "_degen.txt" in hit_file:
        degen_hits_file = hit_file
        degen_controls_file = control_file
        hit_file = hit_file[:-10]
        control_file = control_file[:-10]
    else:
        degen_hits_file = None
        degen_controls_file = None
    hits = parse_basinhoppin_pos(hit_file)
    controls = parse_basinhoppin_pos(control_file)

    if baseml:
        method = "baseml"
        statistic = "tree length"
    else:
        method = "gy"
        statistic = "dS"

    sim_norm_ds = []

    simulations = [i for i in range(n_sim)]
    #parallelize
    result = run_in_parallel(simulations, ["foo", hits, controls, fasta, correspondances, alignments, method, statistic, reverse_site_numbers, degen_hits_file, degen_controls_file], get_sim_p_core)
    for run in result:
        sim_norm_ds.extend(run.get())

    if sim_ds_file:
        rw.write_names([str(i) for i in sim_norm_ds], sim_ds_file)

    p = ms.calc_eff_p(norm_ds, sim_norm_ds, greater = False)
    Z = ms.calc_Zscore(norm_ds, sim_norm_ds)
    sd = np.std(sim_norm_ds)
    print("p: {0}".format(p))
    CI_low = norm_ds - sd
    CI_high = norm_ds + sd
    print("CI: {0} - {1}".format(CI_low, CI_high))
    return(p, CI_low, CI_high, sd, Z)

def get_sim_p_core(simulations, hits, controls, fasta, correspondances, alignments, method, statistic, reverse_site_numbers, degen_hits_file, degen_controls_file):
    '''
    Core function for get_sim_p.
    '''
    sim_norm_ds = []
    counter = 0
    for sim in simulations:
        counter = update_counter(counter, 10)

        if not reverse_site_numbers:
            temp_hits, temp_controls = shuffle_dictionaries(hits, controls)
        else:
            temp_controls, temp_hits = shuffle_dictionaries(hits, controls)

        hit_phylip = "temp_data/temp{0}.phy".format(random.random())
        control_phylip = "temp_data/temp{0}.phy".format(random.random())

        conservation.write_hits_to_phylip(fasta, temp_hits, hit_phylip, correspondances, alignments, degen_hits_file)
        conservation.write_hits_to_phylip(fasta, temp_controls, control_phylip, correspondances, alignments, degen_controls_file)

        hit_ds = conservation.run_codeml(hit_phylip, "temp_data/temp_{0}.phy".format(random.random()), method = method)[statistic]
        control_ds = conservation.run_codeml(control_phylip, "temp_data/temp_{0}.phy".format(random.random()), method = method)[statistic]
        sim_norm_ds.append((hit_ds - control_ds)/control_ds)

        remove_file(hit_phylip)
        remove_file(control_phylip)
    return(sim_norm_ds)    

def main():
    description = "Use different sets of hit and control sites to calculate the normalized dS of a dataset."
    args = parse_arguments(description, ["dataset", "feature_set", "genome", "families_file", "fasta", "hit_file_prefix", "motifs_file", "correspondances", "alignments", "suffix", "trials", "trial_file", "old_trial_file", "region_fasta", "old_motif_format", "nonsense", "no_families", "newest_only", "top_set_only", "calc_p", "reverse_site_numbers", "matched", "degen", "regions"], ints = [10], flags = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23])
    dataset, feature_set, genome, families_file, fasta, hit_file_prefix, motifs_file, correspondances, alignments, suffix, trials, trial_file, old_trial_file, region_fasta, old_motif_format, nonsense, no_families, newest_only, top_set_only, calc_p, reverse_site_numbers, matched, degen, regions = args.dataset, args.feature_set, args.genome, args.families_file, args.fasta, args.hit_file_prefix, args.motifs_file, args.correspondances, args.alignments, args.suffix, args.trials, args.trial_file, args.old_trial_file, args.region_fasta, args.old_motif_format, args.nonsense, args.no_families, args.newest_only, args.top_set_only, args.calc_p, args.reverse_site_numbers, args.matched, args.degen, args.regions

    n_sim = 1000

    print(suffix)

    fs = Feature_Set(feature_set, genome)
    fs.set_dataset(dataset)
    if no_families:
        picked = fs.names
    else:
        families = rw.read_families(families_file)
        fs.add_families(families)
        picked = fs.pick_random_members()       

    hit_phylip = "temp_data/temp_{0}.phy".format(random.random())
    control_phylip = "temp_data/temp_control_{0}.phy".format(random.random())

    if not nonsense:   
        if old_motif_format:
            motifs = rw.read_names(motifs_file)[1:]
        else:
            motifs = rw.read_motifs(motifs_file)
            if top_set_only:
                summary_data = rw.read_many_fields("RBP/RBP_hg38_introncontaining_new.txt", "\t")
                summary_dict = list_to_dict(summary_data, 0, 4, floatify = True)
                motifs = {RBP: motifs[RBP] for RBP in motifs if (summary_dict[RBP] < 0.1)}
            motifs = list(set(flatten(motifs.values())))

    if reverse_site_numbers:
        site_number_suffix = "_reversed_site_numbers_"
    else:
        site_number_suffix = ""

    if matched:
        matched_suff = "_matched"
    else:
        matched_suff = ""

    if degen:
        degen_suff = "_degen.txt"
    else:
        degen_suff = ""

    with open(trial_file, "w") as trial_out:

        trial_out.write("trial\tA\tT\tC\tG\told\told_no_hum_CG\tnew_no_human_CG\tnew_no_hum_no_anc_CG\tnew_w_CG\tnew_no_anc_CG\tnew_no_anc_CG_macaque\tnewer_no_human_CG\tnewer_no_hum_no_anc_CG\tnewer_w_CG\tnewer_no_anc_CG\n")
        if old_trial_file != "None":
            old_trials = rw.read_many_fields(old_trial_file, "\t")
            old_trials = old_trials[1:]
            old_trials = [i[1:5] for i in old_trials]
            seed_kmers = 1
        else:
            seed_kmers = None
        
        for trial in range(trials):

            print(trial)

            trial_output = [trial]

            if nonsense:
                if old_trial_file != "None":
                    scaled_comp = [float(i) for i in old_trials[trial]]
                else:
                    comp = [random.random() for i in range(4)]
                    scaled_comp = [i/np.sum(comp) for i in comp]
                comp_dict = {i: scaled_comp[pos] for pos, i in enumerate(nc._canon_bases_)}
                motifs, obtained_dict = nc.kmers_from_nc(6, 50, comp_dict = comp_dict, return_freqs = True, seed = seed_kmers)
                motifs = ["motifs"] + motifs
                trial_output = trial_output + [obtained_dict[i] for i in nc._canon_bases_]
                temp_motifs_file = "temp_data/temp_motifs.txt"
                rw.write_names(motifs, temp_motifs_file)

            print("===NEW METHOD WITH NO ANCESTRAL CpG (MACAQUE, BIG TREE, CONTEXT), REPLACEMENT CONTROL===")
            hit_file = "{0}_hits_no_anc_CG_only_macaque_big_context{1}_replace.txt{2}".format(hit_file_prefix, matched_suff, degen_suff)
            control_file = "{0}_controls_no_anc_CG_only_macaque_big_context{1}_replace.txt{2}".format(hit_file_prefix, matched_suff, degen_suff)
            if nonsense:
                hit_file = "temp_data/temp_hits{0}.txt".format(random.random())
                control_file = "temp_data/temp_controls{0}.txt".format(random.random())
                error_file = "temp_data/temp_error{0}.txt".format(random.random())
                get_control_sites(fasta, genome, feature_set, families_file, dataset, temp_motifs_file, hit_file, control_file, error_file, "DFE/for_everybody/filtered_hg38_85_pc_multiexon_anc_CG_big_context_threshold05.txt", ["--leave_CG", "--context", "--remove_ancestral_CpG", "--macaque_anc", "--big_tree", "--replacement_control"])
            get_density(fasta, motifs, fs)
            norm_ds = get_new_method_results(hit_file, control_file, hit_phylip, control_phylip, correspondances, alignments, fasta, regions = regions, global_fasta = region_fasta, fs = fs)
            trial_output.append(norm_ds)
            if calc_p:
                p, low_CI, high_CI, sd, Z = get_sim_p(norm_ds, hit_file, control_file, correspondances, alignments, fasta, n_sim, reverse_site_numbers = reverse_site_numbers, sim_ds_file = "{0}{1}_sim_norm_ds_no_anc_CG_only_macaque_big_context{2}_replace.txt{3}".format(hit_file_prefix, site_number_suffix, matched_suff, degen_suff))
                

            trial_output = "\t".join([str(i) for i in trial_output])
            trial_out.write(trial_output)
            trial_out.write("\n")

            remove_file(hit_phylip)

if __name__ == "__main__":
    main()
