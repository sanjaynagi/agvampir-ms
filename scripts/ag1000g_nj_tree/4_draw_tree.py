import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd, json, pickle
from matplotlib.lines import Line2D
out=pickle.load(open('layout.pkl','rb')); L=json.load(open('labels.json')); lab=np.array(L['lab']); kind=np.array(L['kind'])
nodes,leaves,lines=out
assert leaves.id.max()==len(lab)-1 and leaves.id.min()==0
leaves=leaves.sort_values('id').reset_index(drop=True)
col=json.load(open('/Users/sanjinho/Projects/vectors/AmpSeeker-AgVampIR/results-combined/config/metadata_colours.json'))['location']
meta=pd.read_csv('/Users/sanjinho/Projects/vectors/AmpSeeker-AgVampIR/results-combined/config/metadata.tsv',sep='\t').set_index('sample_id')
ref=pd.read_csv('/tmp/ag3_meta_unrestricted.csv').set_index('sample_id')
tcol={'gambiae':'#DDCC77','coluzzii':'#CC6677','arabiensis':'#88CCEE','pwani':'#117733','bissau':'#AA4499','gcx4':'#999933','melas':'#332288','fontenillei':'#888888','merus':'#888888','quadriannulatus':'#888888'}
tname={'gambiae':'An. gambiae','coluzzii':'An. coluzzii','arabiensis':'An. arabiensis','pwani':'Pwani (Gcx3)','bissau':'Bissau','gcx4':'Gcx4','melas':'An. melas','fontenillei':'other (Ag1000G)'}
fig,ax=plt.subplots(figsize=(9,9))
ax.plot(lines.x,lines.y,color='#b0b0b0',lw=0.5,zorder=1)
isref=kind[leaves.id.values]=='ref'
for tx in tcol:
    sel=isref&np.array([ref.loc[s,'taxon']==tx if k=='ref' else False for s,k in zip(lab[leaves.id.values],kind[leaves.id.values])])
    ax.scatter(leaves.x[sel],leaves.y[sel],s=22,c=tcol[tx],marker='s',linewidths=0,alpha=0.9,zorder=2)
am=~isref
c=[col[meta.loc[s,'location']] for s in lab[leaves.id.values][am]]
ax.scatter(leaves.x[am],leaves.y[am],s=34,c=c,edgecolors='black',linewidths=0.4,zorder=3)
ax.set_aspect('equal'); ax.axis('off')
h1=[Line2D([],[],marker='s',ls='',color=tcol[t],label=tname[t],markersize=7) for t in ['gambiae','coluzzii','arabiensis','pwani','bissau','gcx4','melas','fontenillei']]
h2=[Line2D([],[],marker='o',ls='',markerfacecolor=col[l],markeredgecolor='black',label=l.replace('_',' '),markersize=7) for l in col]
l1=ax.legend(handles=h1,title='Ag1000G reference',loc='upper left',frameon=False,fontsize=8,title_fontsize=9)
ax.add_artist(l1); ax.legend(handles=h2,title='Ag-vampIR cohort',loc='lower left',frameon=False,fontsize=8,title_fontsize=9)
dst='/Users/sanjinho/Projects/vectors/AmpSeeker-AgVampIR/results-combined/'
plt.savefig(dst+'figure2C_ag1000g_nj.png',dpi=300,bbox_inches='tight'); plt.savefig(dst+'figure2C_ag1000g_nj.svg',bbox_inches='tight')
print('saved',(~isref).sum(),'amplicon and',isref.sum(),'reference tips')
