"""Mean and SD of AgamP4 mapping rate per taxon (Reviewer 2), for a supplementary table.

Needs per-sample `samtools flagstat` output, which is not in results-combined (BAMs are not stored). Generate it with, for example:
    for b in results/alignments/*.bam; do samtools flagstat "$b" > flagstat/$(basename "$b" .bam).flagstat; done
Usage: python mapping_rates_by_taxon.py flagstat_dir "Supplementary Data 2.tsv" out.tsv
"""
import re, sys
from pathlib import Path
import pandas as pd

fdir, meta, out = sys.argv[1:4]
rows = []
for f in Path(fdir).glob("*.flagstat"):
    m = re.search(r" mapped \((\d+\.\d+)%", f.read_text())  # first match is the all-reads line
    rows.append({"sample_id": f.stem, "mapped_pct": float(m.group(1))})
d = pd.DataFrame(rows).merge(pd.read_csv(meta, sep="\t"), on="sample_id")
d = d[d.qc_pass]
res = (d.groupby(d.taxon.fillna("unassigned")).mapped_pct.agg(n="size", mean_mapped_pct="mean", sd_mapped_pct="std").round(2))
res.to_csv(out, sep="\t"); print(res)
