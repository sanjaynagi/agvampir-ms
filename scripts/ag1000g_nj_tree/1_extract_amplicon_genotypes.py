"""Same site selection as extract_amp_sites.py, but also keep per-genotype depth (DP) for a depth filter."""
import numpy as np, pandas as pd, allel, os
RES="/mnt/user_shares/snagi/lstm_scratch/network_scratch/ampseq-agvampir-combined/results"
OUT=os.path.expanduser("~/ag3-tree"); meta=pd.read_csv(f"{RES}/config/metadata.qcpass.tsv",sep="\t")
v=allel.read_vcf(f"{RES}/vcfs/amplicons/ag-vampir-combined.annot.vcf",fields=["samples","variants/CHROM","variants/POS","variants/REF","variants/ALT","calldata/GT","calldata/AD"])
samples=v["samples"]; keep=np.isin(samples,meta.sample_id)
gt=allel.GenotypeArray(v["calldata/GT"]).compress(keep,axis=1); dp=np.clip(v["calldata/AD"],0,None).sum(axis=2)[:,keep]; samples=samples[keep]
chrom=v["variants/CHROM"]; pos=v["variants/POS"]; ref=v["variants/REF"]; alt=v["variants/ALT"]
n_alt=(alt!="").sum(axis=1)
snp=(np.char.str_len(ref.astype(str))==1)&(np.char.str_len(alt[:,0].astype(str))==1)&(n_alt==1)&np.isin(chrom,["2L","2R","3L","3R","X"])
miss=gt.is_missing().mean(axis=1); seg=gt.count_alleles(max_allele=1).is_segregating()
m=snp&(miss<=0.05)&seg
g=gt.compress(m,axis=0); d=dp[m]
dos=np.where(g.is_missing(),-1,g.to_n_alt()).astype(np.int8)
np.savez_compressed(f"{OUT}/amp_dosage_dp.npz",dos=dos,dp=np.clip(d,-1,32000).astype(np.int16),samples=samples.astype(str),chrom=chrom[m].astype(str),pos=pos[m],ref=ref[m].astype(str),alt=alt[m,0].astype(str))
print("saved",dos.shape)
