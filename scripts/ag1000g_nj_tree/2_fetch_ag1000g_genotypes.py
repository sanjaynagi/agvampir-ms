"""NJ tree of Ag-vampIR amplicon samples together with Ag1000G (Ag3) reference samples, using amplicon SNP sites only.
Reference samples are drawn only from sample sets with unrestricted use."""
import matplotlib as mpl; mpl.use("Agg")
import json, sys, os
import numpy as np, pandas as pd
import malariagen_data
HERE=os.path.dirname(os.path.abspath(__file__)); rng=np.random.default_rng(42)

# ---- amplicon genotypes
z=np.load(f"{HERE}/amp_dosage.npz"); dos_amp=z["dos"].astype(np.int8); amp_samples=z["samples"]; chrom=z["chrom"]; pos=z["pos"]; ref_a=z["ref"]; alt_a=z["alt"]
meta=pd.read_csv("/Users/sanjinho/Projects/vectors/AmpSeeker-AgVampIR/results-combined/config/metadata.tsv",sep="\t").set_index("sample_id").loc[amp_samples]
amp_sub=[]
for loc,g in meta.groupby("location"):
    ids=g.index.to_numpy(); amp_sub+=list(rng.choice(ids,size=min(25,len(ids)),replace=False))
amp_idx=[int(np.where(amp_samples==s)[0][0]) for s in amp_sub]
print("amplicon samples used:",len(amp_sub),"SNP sites:",len(pos))

# ---- reference samples (unrestricted Ag3 sample sets only)
ag3=malariagen_data.Ag3(results_cache="/tmp/ag3_cache",pre=True)
m=pd.read_csv("/tmp/ag3_meta_unrestricted.csv")
caps={"gambiae":(4,60),"coluzzii":(4,50),"arabiensis":(4,50),"pwani":(0,30),"bissau":(0,15),"gcx4":(0,8),"melas":(0,8),"fontenillei":(0,3),"quadriannulatus":(0,3),"merus":(0,3)}
ref_ids=[]
for tx,(per_c,cap) in caps.items():
    d=m[m.taxon==tx]
    if per_c:
        d=d.sample(frac=1,random_state=1).groupby("country").head(per_c)
    ref_ids+=list(d.sample(n=min(cap,len(d)),random_state=2).sample_id)
ref=m[m.sample_id.isin(ref_ids)].copy(); print("reference samples:",len(ref)); print(ref.taxon.value_counts().to_dict())
sets=sorted(ref.sample_set.unique())

# ---- windows around amplicon sites
df=pd.DataFrame({"c":chrom,"p":pos,"i":np.arange(len(pos))}).sort_values(["c","p"])
df["w"]=((df.groupby("c").p.diff().fillna(1e9))>1500).cumsum()
dos_ref=np.full((len(pos),len(ref_ids)),-1,dtype=np.int8); keep_site=np.zeros(len(pos),bool); order=None
other_n=np.zeros(len(pos),int); called_n=np.zeros(len(pos),int)
for w,g in df.groupby("w"):
    c=g.c.iloc[0]; s,e=int(g.p.min()),int(g.p.max())
    ds=ag3.snp_calls(region=f"{c}:{s}-{e}",sample_sets=sets,sample_query="sample_id in @ids",sample_query_options={"local_dict":{"ids":ref_ids}},site_mask=None)
    sid=ds["sample_id"].values.astype(str)
    if order is None: order=sid
    vp=ds["variant_position"].values; va=ds["variant_allele"].values.astype(str); gt=ds["call_genotype"].values  # variants, samples, ploidy
    col=np.array([int(np.where(order==x)[0][0]) for x in sid])
    for i,p in zip(g.i,g.p):
        k=np.searchsorted(vp,p)
        if k>=len(vp) or vp[k]!=p: continue
        al=va[k]
        
        if (al==ref_a[i]).sum()==0: continue
        a=gt[k]  # samples, ploidy
        miss=(a<0).any(axis=1)
        idx=np.clip(a,0,len(al)-1); base=al[idx]                      # samples x ploidy, allele strings
        other=((base!=ref_a[i])&(base!=alt_a[i])&(a>=0)).any(axis=1)  # carries a third base
        cnt=(base==alt_a[i]).sum(axis=1).astype(np.int8)
        cnt[miss]=-1; cnt[other]=-1                                    # third-base calls treated as missing
        other_n[i]=int(other.sum()); called_n[i]=int((~miss).sum())
        dos_ref[i,col]=cnt; keep_site[i]=True
    print("window",w,c,s,e,"done",flush=True)
np.savez_compressed(f"{HERE}/ref_dosage.npz",dos=dos_ref,order=order,keep=keep_site,other_n=other_n,called_n=called_n)
