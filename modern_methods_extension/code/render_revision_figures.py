"""Render revision figures from recorded outputs, preserving failures and denominators."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from matplotlib.patches import FancyBboxPatch
ext=Path(__file__).resolve().parent.parent; package=ext.parents[1]
out=ext/'revision_figures'; out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
manifest=[]
def save(fig,name,sources):
    for suffix in ('png','svg'):
        path=out/f'{name}.{suffix}'; fig.savefig(path,dpi=600,facecolor='white')
        manifest.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
          sources=[str(s.relative_to(package)) for s in sources],figure_size_inches=fig.get_size_inches().tolist(),dpi=600))
    plt.close(fig)

# Fixed fixtures are not independent stochastic observations.
fig,ax=plt.subplots(figsize=(6.69,3.7),layout='constrained'); ax.axis('off')
labels=['Core contract tests','Distinct fault fixtures','Repeated determinism checks','Cross-decomposition comparisons',
        'Historical consumer experiment','CAD/HF technical attempts','Independent biological validation']
values=['93 / 93','25 unique fixtures','25 fixtures x 100 repeats','80 rank comparisons','100 data pairs x 5 conditions',
        '192; model failures retained','Not performed']
for i,(label,value) in enumerate(zip(labels,values)):
    yy=1-(i+.6)/len(labels); ax.text(.02,yy,label,va='center'); ax.text(.55,yy,value,va='center',weight='bold',fontsize=8.5)
    ax.axhline(yy-.05,color='#dddddd',linewidth=.6)
save(fig,'Fig3',[package/'manuscript/main.tex'])

counts_path=package/'supplementary/code_package/ESM_1_Code_and_Reproducibility_Package/external_validation/audit_and_results/formal_model_terminal_counts.tsv'
counts=pd.read_csv(counts_path,sep='\t')
fig,axes=plt.subplots(1,2,figsize=(6.69,3.6),layout='constrained',gridspec_kw={'width_ratios':[1,1.7]})
axes[0].axis('off'); axes[0].set_title('a',loc='left',weight='bold')
axes[0].text(.02,.86,'24 regions\n4 GTEx tissues\n2 GWAS outcomes',va='top',linespacing=1.6)
axes[0].text(.02,.39,'192 / 192 READY\n192 exact reruns\n15 comparable outputs',va='top',linespacing=1.6)
codes=['OK','OK_NO_COMPARABLE_CREDIBLE_SETS','E_SUSIE_NONCONVERGENCE']
names=['Comparable CS pairs','No comparable CS pairs','Nonconvergent']
colours=['#0072B2','#8CBED6','#D55E00']; hatches=['','///','xx']; bottom=np.zeros(2)
for code,label,colour,hatch in zip(codes,names,colours,hatches):
    values=np.array([counts.loc[(counts.stratum==s)&(counts.stable_code==code),'count'].sum() for s in ['CAD','HF']])
    axes[1].bar(['CAD','HF'],values,bottom=bottom,label=label,color=colour,hatch=hatch,edgecolor='white')
    for j,n in enumerate(values):
        if n: axes[1].text(j,bottom[j]+n/2,str(n),ha='center',va='center',fontsize=9)
    bottom+=values
axes[1].set(ylim=(0,103),ylabel='Technical attempts'); axes[1].set_title('b',loc='left',weight='bold')
axes[1].legend(loc='upper center',bbox_to_anchor=(.45,-.1),frameon=False,fontsize=7.5)
save(fig,'Fig7',[counts_path])

source=ext/'ld_guardnet/ood_predictions.tsv'; data=pd.read_csv(source,sep='\t')
classes=['valid','order_mismatch','sign_corruption','ld_mismatch']
y=data.true_class.to_numpy(); p=data[['prob_'+c for c in classes]].to_numpy()
fig=plt.figure(figsize=(6.69,5.5),layout='constrained')
panel=fig.subplot_mosaic([['a','b'],['c','c']])
axes=[panel[key] for key in ('a','b','c')]
axes[0].axis('off'); axes[0].set_title('a',loc='left',weight='bold')
for i,text in enumerate(['Five node features','Signed LD channels','Two 48-unit blocks','Mean + max pooling','Four-state softmax']):
    yy=.93-i*.19
    axes[0].text(.5,yy,text,ha='center',va='center',bbox=dict(boxstyle='round',fc='#EAF2F8',ec='#555555'),fontsize=8)
    if i<4: axes[0].annotate('',xy=(.5,yy-.13),xytext=(.5,yy-.055),arrowprops=dict(arrowstyle='->'))
cm=confusion_matrix(y,p.argmax(1),normalize='true')
axes[1].imshow(cm,vmin=0,vmax=1,cmap='Blues'); axes[1].set_title('b',loc='left',weight='bold')
labels=['Valid','Value order','Sign error','LD mismatch']
axes[1].set_xticks(range(4),labels,rotation=38,ha='right',fontsize=8); axes[1].set_yticks(range(4),labels,fontsize=8)
axes[1].set(xlabel='Predicted',ylabel='Generated state')
for i in range(4):
    for j in range(4): axes[1].text(j,i,f'{cm[i,j]:.2f}',ha='center',va='center',color='white' if cm[i,j]>.55 else 'black',fontsize=8)
axes[2].set_title('c',loc='left',weight='bold'); axes[2].plot([0,1],[0,1],'--',color='gray',linewidth=.8)
bin_rows=[]; groups=list(data.groupby('replicate').indices.values()); rng=np.random.default_rng(411)
boot_indices=[np.concatenate([groups[k] for k in rng.integers(0,len(groups),len(groups))]) for _ in range(500)]
for prefix,label,marker,colour in [('prob_','LD-GuardNet','o','#0072B2'),('baseline_prob_','Logistic','s','#D55E00')]:
    scores=1-data[prefix+'valid'].to_numpy(); truth=(y!=0).astype(float)
    edges=np.unique(np.quantile(scores,np.linspace(0,1,11))); bins=np.digitize(scores,edges[1:-1],right=True)
    xs=[]; means=[]; lows=[]; highs=[]
    for k in sorted(set(bins)):
        mask=bins==k; estimates=[]
        for ix in boot_indices:
            subset=ix[bins[ix]==k]
            if len(subset): estimates.append(truth[subset].mean())
        lo,hi=np.quantile(estimates,[.025,.975]); xs.append(scores[mask].mean()); means.append(truth[mask].mean()); lows.append(lo); highs.append(hi)
        bin_rows.append(dict(model=label,bin=k,n=int(mask.sum()),mean_prediction=xs[-1],observed=means[-1],low=lo,high=hi))
    axes[2].errorbar(xs,means,yerr=[np.maximum(0,np.array(means)-lows),np.maximum(0,np.array(highs)-means)],
                     marker=marker,color=colour,label=label,linewidth=.8,markersize=3,capsize=1)
axes[2].set(xlim=(0,1),ylim=(0,1),xlabel='Predicted anomaly score',ylabel='Observed anomaly fraction'); axes[2].legend(frameon=False,fontsize=7)
pd.DataFrame(bin_rows).to_csv(out/'Fig9_calibration_bins.tsv',sep='\t',index=False)
save(fig,'Fig9',[source])

summary_path=ext/'revision_validation/all_model_challenge_summary.tsv'; s=pd.read_csv(summary_path,sep='\t')
rates_path=ext/'revision_validation/all_model_region_fault_metrics.tsv'; rates=pd.read_csv(rates_path,sep='\t')
models=['signed_graph','unsigned_graph','pooled_mlp','logistic','susie_kriging']
model_labels=['Signed LD graph','Unsigned LD graph','Pooled MLP','Logistic','SuSiE kriging']
faults=['valid','sign_1pct','sign_5pct','sign_10pct','sign_20pct','block_permutation','random_permutation','ld_scramble']
fault_labels=['Valid\n(FPR)','Sign\n1%','Sign\n5%','Sign\n10%','Sign\n20%','Block\nvalues','Random\nvalues','LD\ncontent']
fig,axes=plt.subplots(2,1,figsize=(6.69,5.4),layout='constrained',gridspec_kw={'height_ratios':[1,1.45]})
for i,model in enumerate(models):
    row=s.set_index('model').loc[model]
    axes[0].errorbar(row.anomaly_auc,4-i,xerr=[[row.anomaly_auc-row.region_bootstrap_low],[row.region_bootstrap_high-row.anomaly_auc]],fmt='o',color='#0072B2',capsize=3)
axes[0].set_yticks(range(5),model_labels[::-1]); axes[0].set(xlim=(.5,1),xlabel='Anomaly AUROC (region-bootstrap 95% interval)'); axes[0].set_title('a',loc='left',weight='bold')
matrix=np.array([[rates.loc[(rates.model==model)&(rates.fault==fault),'flag_rate'].mean() for fault in faults] for model in models])
axes[1].imshow(matrix,vmin=0,vmax=1,cmap='Blues',aspect='auto'); axes[1].set_title('b',loc='left',weight='bold')
axes[1].set_xticks(range(8),fault_labels,fontsize=7.5); axes[1].set_yticks(range(5),model_labels,fontsize=8)
for i in range(5):
    for j in range(8): axes[1].text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if matrix[i,j]>.6 else 'black')
axes[1].set_xlabel('Flagged fraction; 440 examples/state; thresholds from validation-valid inputs')
save(fig,'Fig10',[summary_path,rates_path])
(out/'figure_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Rendered Fig3, Fig7, Fig9, Fig10 with source tables and manifests.')
