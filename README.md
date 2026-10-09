# agvampir-ms

Code, configuration, supplementary material and the results book for the manuscript
**"Targeted genomic surveillance of insecticide resistance in African malaria vectors"** (Nagi et al.).

**Results book:** https://sanjaynagi.github.io/agvampir-ms/

The analysis uses [AmpSeeker](https://github.com/sanjaynagi/AmpSeeker) (commit `864ba69`) with the
Ag-vampIR panel on 933 samples from ten cohorts (two field populations and eight laboratory colonies;
784 samples pass QC). The earlier four-cohort book is archived at
https://sanjaynagi.github.io/agvampir002-results/.

## Contents

| Path | What it holds |
|---|---|
| `config/` | AmpSeeker configuration for the combined run (`config.yaml`), sample metadata, target BED file and colour scheme |
| `notebooks/` | The AmpSeeker analysis notebooks as executed for this paper (QC, population structure, allele frequencies, species ID, kdr analysis), with outputs cleared; the outputs are in the results book |
| `analysis-notebooks/` | Additional analysis notebooks (taxon classifier, phenotype plots, VCF-to-array loading) |
| `scripts/` | `34mb_gene_annotation.py` (Supplementary Data 4) and `mapping_rates_by_taxon.py` (mapping rates in Supplementary Table 2) |
| `figures/` | Main Figures 2 and 3A–B |
| `supplement/` | Supplementary Tables 1–6, Data 1–4, Figures 1–6 and Texts 1–3 |

## Reproducing the analysis

1. Clone AmpSeeker and check out `864ba69`.
2. Place the FASTQs in `resources/reads/{sample_id}_{1,2}.fastq.gz` and copy `config/` into the AmpSeeker `config/` directory.
3. Run the workflow as described in the AmpSeeker documentation. The results book is built by the workflow.
4. Mapping rates: run `samtools flagstat` on each BAM, then `python scripts/mapping_rates_by_taxon.py flagstat_dir "supplement/Supplementary Data 2.tsv" out.tsv`.
5. Gene annotation: `python scripts/34mb_gene_annotation.py AgamP4.12.gff3 "supplement/Supplementary Data 3.xlsx" out.tsv`.

## Data

Raw sequencing reads are not stored here (SRA BioProject PRJNA1207724; see the manuscript's Data availability statement).
Per-sample taxon assignments and QC status are in `supplement/Supplementary Data 2.tsv`.

## Not yet included

The genotype–phenotype association and bioassay prediction notebooks for the Siaya colony.
