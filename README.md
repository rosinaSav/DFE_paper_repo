This is the source code for Savisaar and Hurst paper on the distribution of fitness effects at ESE sites.

bedtools_games.py, conservation.py, ensembl_ops.py, housekeeping.py, nucleotide_comp.py, 
read_and_write.py and my_stats.py are custom Python modules. Most of the code in these modules is not relevant to the paper.
However, most of the scripts used to perform the analyses described in the paper use one or more of the functions contained in this module.
This is why these modules are included in the repository.

INSIGHT.py: for performing INSIGHT analysis

Keightley.py: for performing multiDFE analysis

Keightley_manual.py: for comparing the frequency of segregating sites and minor allele frequencies between ESE and control sites

Keightley_stability.py: to determine the false positive rate of multiDFE

convert_region_to_CDS.py: to convert motif hit/control positions specified relative to an exonic subregion into coordinates relative to a full CDS

degenerate_SFS.py: to modify an SFS so as to remove the effects of motif degeneracy. Also creates hit/control position files that specify the motif degeneracy of the possible substitutions at each site.

get_core_and_flank_files.py: to extract exonic subregion sequences and coordinates.

mDFEest_input.py: to prepare input file for multiDFE.

norm_ds_with_new_controls.py: to calculate normalized dS and motif density.

pick_control_sites.py: to pick control sites to go with a set of motif hits sites

shuffle_hits_and_controls: given a set of motif hit and control sites, generate new files where the hits and controls have been shuffled.

