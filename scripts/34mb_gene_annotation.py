"""Gene annotation of the 34 Mb sweep locus on 2L (AgamP4.12), for Supplementary Data 4.

Region: span of the 14 FDR-significant SNPs at the locus in the Siaya association study
(Supplementary Data 3), 2L:34,100,727-34,118,260, plus 100 kb flanks. Genes are taken from
the AgamP4.12 BASEFEATURES GFF3. SNP consequences come from the ANN field in Supplementary Data 3.
Usage: python 34mb_gene_annotation.py GFF3 "Supplementary Data 3.xlsx" out.tsv
"""
import re, sys
import pandas as pd

gff, sd2, out = sys.argv[1:4]
CHROM, FLANK = "2L", 100_000

def go_names(ids):
    """GO term names from the QuickGO API (empty labels if offline)."""
    import json, urllib.request
    if not ids: return {}
    req = urllib.request.Request("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + ",".join(sorted(ids)),
                                 headers={"Accept": "application/json"})
    try:
        return {r["id"]: r["name"] for r in json.load(urllib.request.urlopen(req, timeout=30))["results"]}
    except Exception as e:
        print("GO label lookup failed:", e); return {}

d = pd.read_excel(sd2)
d["chrom"] = d.snp.str.extract(r"snp_(\w+):")[0]
d["pos"] = d.snp.str.extract(r":(\d+)_")[0].astype(int)
loc = d[(d.chrom == CHROM) & (d.pos.between(34_100_000, 34_119_000))].copy()
sig = loc[loc.fdr_sig]
lo, hi = int(sig.pos.min()), int(sig.pos.max())
wlo, whi = lo - FLANK, hi + FLANK

def attrs(s):
    return dict(kv.split("=", 1) for kv in s.strip().split(";") if "=" in kv)

genes, mrna = {}, {}
with open(gff) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[0] != CHROM: continue
        s, e = int(p[3]), int(p[4])
        if e < wlo or s > whi: continue
        a = attrs(p[8])
        if p[2] in ("gene", "ncRNA_gene"):
            genes[a["ID"]] = dict(gene_id=a["ID"], start=s, end=e, strand=p[6], biotype=a.get("biotype", ""),
                                  description=re.sub(r"\s*\[Source:.*\]", "", a.get("description", "")))
        elif p[2] == "mRNA":
            mrna.setdefault(a["Parent"], set()).update(a.get("Ontology_term", "").split(",") if a.get("Ontology_term") else [])

rows = []
for g in sorted(genes.values(), key=lambda x: x["start"]):
    inside = loc[(loc.pos >= g["start"]) & (loc.pos <= g["end"])]
    sig_inside = inside[inside.fdr_sig]
    dist = 0 if ((sig.pos >= g["start"]) & (sig.pos <= g["end"])).any() else int(min(abs(sig.pos - g["start"]).min(), abs(sig.pos - g["end"]).min()))
    terms = sorted(mrna.get(g["gene_id"], set()))
    ann = []
    for _, r in inside.iterrows():
        a = str(r.ANN).split(",")[0].split("|")
        ann.append(f"{r.pos}:{a[1]}" + (f" ({a[10]})" if len(a) > 10 and a[10] else ""))
    rows.append({**g, "length_bp": g["end"] - g["start"] + 1,
                 "overlaps_significant_snp": len(sig_inside) > 0,
                 "n_significant_snps_in_gene": len(sig_inside),
                 "distance_to_nearest_significant_snp_bp": dist,
                 "go_terms": ";".join(terms),
                 "genotyped_snp_consequences": "; ".join(ann)})
GO = go_names({t for r in rows for t in re.findall(r"GO:\d+", r["go_terms"])})
for r in rows:
    r["go_terms"] = "; ".join(f"{t} {GO.get(t, '')}".strip() for t in re.findall(r"GO:\d+", r["go_terms"]))
res = pd.DataFrame(rows)
res.insert(0, "contig", CHROM)
res.to_csv(out, sep="\t", index=False)
print(f"Associated-SNP span {CHROM}:{lo:,}-{hi:,}; window {wlo:,}-{whi:,}; {len(res)} genes; "
      f"{int(res.overlaps_significant_snp.sum())} overlap significant SNPs")
