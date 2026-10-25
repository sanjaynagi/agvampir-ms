# Neighbour-joining tree of Ag-vampIR samples with Ag1000G reference samples

Uses only the amplicon SNP sites (biallelic, at most 5% missing, segregating in the QC-passed samples; 2,901 sites, 2,758 used after missingness filters in the combined set).
Reference samples are from Ag1000G (Ag3) sample sets with unrestricted use only, subsampled by taxon and country (229 samples).
Ag-vampIR samples are subsampled to at most 25 per cohort (250). Distance is the mean absolute difference in ALT-allele dosage over jointly called sites; the tree is built with `anjl` (canonical NJ) and drawn with an equal-angle layout.

1. `1_extract_amplicon_genotypes.py` (run on the server with the amplicon VCF) writes `amp_dosage.npz`
2. `2_fetch_ag1000g_genotypes.py` pulls the reference genotypes with `malariagen_data` (needs `/tmp/ag3_meta_unrestricted.csv`, the sample metadata for the unrestricted sets)
3. `3_distance_and_nj.py` computes distances and the NJ layout
4. `4_draw_tree.py` draws `figures/figure2C_ag1000g_nj.png`

Paths are as run for the paper and will need editing.
