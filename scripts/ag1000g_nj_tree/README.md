# Neighbour-joining tree of Ag-vampIR samples with Ag1000G reference samples

Sites: biallelic amplicon SNPs, at most 5% missing and segregating in the QC-passed amplicon samples (2,901), then filtered to at most 2% missing in both sets after masking amplicon genotypes with depth below 10 (from allelic depths) and Ag1000G calls carrying a third base (1,765 sites).
Samples: up to 25 per cohort (203 kept after a 10% per-sample missingness filter) and 229 Ag1000G reference samples from sample sets with unrestricted use only, subsampled by taxon and country.
Distance: mean absolute difference in ALT-allele dosage over jointly called sites. Tree: canonical neighbour-joining (`anjl`), equal-angle layout.

1. `1_extract_amplicon_genotypes.py` (run on the server with the amplicon VCF) writes `amp_dosage_dp.npz`
2. `2_fetch_ag1000g_genotypes.py` pulls reference genotypes with `malariagen_data` (needs the metadata of unrestricted sets as `/tmp/ag3_meta_unrestricted.csv`)
3. `3_distance_and_nj.py MIN_DP SITE_MISSING` filters, computes distances and the NJ layout, and prints nearest-reference diagnostics
4. `4_draw_tree.py` draws `figures/figure2C_ag1000g_nj.png`
5. `5_fetch_matched_obuasi_wgs.py` fetches the matched whole-genome data for 40 Obuasi individuals (amplicon vs WGS concordance: 98.8% of genotypes agree, 99.2% at depth 10 or more; each individual's own WGS is its nearest neighbour for 30 of 30)

Paths are as run for the paper and will need editing.
