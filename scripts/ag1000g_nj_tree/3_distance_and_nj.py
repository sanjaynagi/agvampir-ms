import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd, json, anjl
HERE="/private/tmp/ag3tree"; rng=np.random.default_rng(42)
z=np.load(f"{HERE}/amp_dosage.npz"); da=z["dos"]; samples=z["samples"]
r=np.load(f"{HERE}/ref_dosage.npz",allow_pickle=True); dr=r["dos"]; order=r["order"].astype(str); keep=r["keep"]
meta=pd.read_csv("/Users/sanjinho/Projects/vectors/AmpSeeker-AgVampIR/results-combined/config/metadata.tsv",sep="\t").set_index("sample_id").loc[samples]
amp_sub=[]
for loc,g in meta.groupby("location"):
    ids=g.index.to_numpy(); amp_sub+=list(rng.choice(ids,size=min(25,len(ids)),replace=False))
ai=[int(np.where(samples==s)[0][0]) for s in amp_sub]
refm=pd.read_csv("/tmp/ag3_meta_unrestricted.csv").set_index("sample_id").loc[order]
# the ref window loop saved columns in `order`
sites=keep & ((dr<0).mean(axis=1)<=0.10) & ((da[:,ai]<0).mean(axis=1)<=0.10)
print("sites used:",sites.sum(),"of",len(keep))
D=np.hstack([da[sites][:,ai],dr[sites]]).astype(np.int16).T   # samples x sites
lab=list(amp_sub)+list(order); kind=["amp"]*len(amp_sub)+["ref"]*len(order)
valid=(D>=0); n=len(D); dist=np.zeros((n,n),dtype=np.float32)
for i in range(n):
    both=valid[i]&valid
    diff=np.abs(D[i]-D)*both
    dist[i]=diff.sum(axis=1)/np.maximum(both.sum(axis=1),1)
dist=(dist+dist.T)/2; np.fill_diagonal(dist,0)
Z=anjl.canonical_nj(dist)
out=anjl.layout_equal_angle(Z)
print(type(out),[type(o) for o in out]); [print(o.columns.tolist(), len(o)) for o in out]
np.save(f"{HERE}/dist.npy",dist); json.dump({"lab":lab,"kind":kind},open(f"{HERE}/labels.json","w"))
import pickle; pickle.dump(out,open(f"{HERE}/layout.pkl","wb"))
