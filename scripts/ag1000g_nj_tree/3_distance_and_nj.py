import matplotlib as mpl; mpl.use("Agg")
import numpy as np, pandas as pd, json, pickle, anjl, sys
HERE="/private/tmp/ag3tree"; MIN_DP=int(sys.argv[1]) if len(sys.argv)>1 else 10; SITE_MISS=float(sys.argv[2]) if len(sys.argv)>2 else 0.02; SAMPLE_MISS=0.10
rng=np.random.default_rng(42)
b=np.load(f"{HERE}/amp_dosage_dp.npz"); da=b["dos"].astype(np.int16); dp=b["dp"]; samples=b["samples"]
da[(dp<MIN_DP)]=-1                                     # low-depth amplicon genotypes -> missing
r=np.load(f"{HERE}/ref_dosage.npz",allow_pickle=True); dr=r["dos"].astype(np.int16); order=r["order"].astype(str); keep=r["keep"]
meta=pd.read_csv("/Users/sanjinho/Projects/vectors/AmpSeeker-AgVampIR/results-combined/config/metadata.tsv",sep="\t").set_index("sample_id").loc[samples]
amp_sub=[]
for loc,g in meta.groupby("location"):
    ids=g.index.to_numpy(); amp_sub+=list(rng.choice(ids,size=min(25,len(ids)),replace=False))
ai=[int(np.where(samples==s)[0][0]) for s in amp_sub]
A=da[:,ai]; R=dr
cand=keep&((R<0).mean(axis=1)<=0.20)&((A<0).mean(axis=1)<=0.20)
# drop poorly genotyped samples on candidate sites, then apply the strict site filter on the remaining samples
sa=(A[cand]<0).mean(axis=0)<=SAMPLE_MISS; sr=(R[cand]<0).mean(axis=0)<=SAMPLE_MISS
print(f"MIN_DP={MIN_DP} site missingness<={SITE_MISS}; samples kept: amplicon {sa.sum()}/{len(sa)}, reference {sr.sum()}/{len(sr)}")
A2=A[:,sa]; R2=R[:,sr]; amp2=[s for s,k in zip(amp_sub,sa) if k]; ref2=[s for s,k in zip(order,sr) if k]
sites=keep&((R2<0).mean(axis=1)<=SITE_MISS)&((A2<0).mean(axis=1)<=SITE_MISS)
print("sites used:",sites.sum(),"of",len(sites))
D=np.hstack([A2[sites],R2[sites]]).T; lab=amp2+ref2; kind=["amp"]*len(amp2)+["ref"]*len(ref2)
valid=D>=0; n=len(D); dist=np.zeros((n,n),dtype=np.float32)
for i in range(n):
    both=valid[i]&valid; dist[i]=(np.abs(D[i]-D)*both).sum(1)/np.maximum(both.sum(1),1)
dist=(dist+dist.T)/2; np.fill_diagonal(dist,0)
# diagnostics: nearest reference sample of the expected taxon vs reference-reference nearest neighbour
refm=pd.read_csv("/tmp/ag3_meta_unrestricted.csv").set_index("sample_id")
rt=np.array([refm.loc[s,"taxon"] for s in ref2]); na=len(amp2)
dar=dist[:na,na:]; drr=dist[na:,na:].copy(); np.fill_diagonal(drr,np.inf)
print("median nearest-neighbour distance among reference samples:",round(float(np.median(drr.min(1))),4))
print("median nearest reference distance for each cohort (any reference | closest taxon):")
for loc in meta.loc[amp2].location.unique():
    idx=np.where(meta.loc[amp2].location.values==loc)[0]
    nn=dar[idx].min(1); tx=pd.Series(rt[dar[idx].argmin(1)]).value_counts().to_dict()
    print(f"  {loc:11s} n={len(idx):2d} nearest-ref {np.median(nn):.4f}  nearest taxa {tx}")
het_a=((A2[sites]==1).sum(0)/np.maximum((A2[sites]>=0).sum(0),1)).mean(); het_r=((R2[sites]==1).sum(0)/np.maximum((R2[sites]>=0).sum(0),1)).mean()
print(f"mean heterozygous fraction: amplicon {het_a:.3f}, reference {het_r:.3f}")
Z=anjl.canonical_nj(dist); out=anjl.layout_equal_angle(Z)
np.save(f"{HERE}/dist_strict.npy",dist); json.dump({"lab":lab,"kind":kind},open(f"{HERE}/labels.json","w")); pickle.dump(out,open(f"{HERE}/layout.pkl","wb"))
