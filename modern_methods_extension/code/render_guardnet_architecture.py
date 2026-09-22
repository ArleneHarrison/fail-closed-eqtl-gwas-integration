"""Code-grounded, editable vector diagram of the implemented LD-GuardNet."""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import torch

HERE=Path(__file__).resolve().parent
SOURCE=HERE/"run_ld_guardnet_benchmark.py"
spec=importlib.util.spec_from_file_location("guardnet_architecture_source",SOURCE)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
OUT=HERE.parent/"architecture_v4";OUT.mkdir(exist_ok=True)
V=HERE.parents[2]

# Verify tensor dimensions and parameter count against the actual network.
model=m.SignedLDGuardNet(5,48,.2).eval()
shapes={}
for name,layer in (("block1",model.layer1),("block2",model.layer2),("classifier",model.classifier)):
    layer.register_forward_hook(lambda module,inputs,output,n=name: shapes.update({n:{"input":list(inputs[0].shape),"output":list(output.shape)}}))
eye=torch.eye(100).unsqueeze(0)
with torch.no_grad():
    logits=model(torch.zeros(1,100,5),torch.tensor([0]),eye,eye)
assert shapes["block1"]=={"input":[1,100,15],"output":[1,100,48]}
assert shapes["block2"]=={"input":[1,100,144],"output":[1,100,48]}
assert shapes["classifier"]=={"input":[1,96],"output":[1,4]}
assert sum(p.numel() for p in model.parameters())==12772
probe=np.array([[1.,.3,-.2],[.3,1.,.1],[-.2,.1,1.]])
saved_panels=m.PANELS;m.PANELS=("audit",)
ap,an=m.signed_operators({"audit":probe});m.PANELS=saved_panels
for sign,actual in ((1,ap[0].numpy()),(-1,an[0].numpy())):
    raw=probe.copy();np.fill_diagonal(raw,0)
    w=np.maximum(sign*raw,0)+np.eye(3)
    inv=1/np.sqrt(w.sum(1))
    np.testing.assert_allclose(actual,inv[:,None]*w*inv[None,:],atol=1e-7)

INK="#203448";BLUE="#176A96";ORANGE="#A85418";TEAL="#14665C";GREY="#63717E"
with plt.rc_context({"font.family":"DejaVu Sans","font.size":7.2,
                     "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
                     "mathtext.fontset":"dejavusans","axes.unicode_minus":False}):
    fig=plt.figure(figsize=(170/25.4,176/25.4),facecolor="white")
    ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,170);ax.set_ylim(176,0);ax.axis("off")
    def txt(x,y,s,size=7.2,weight="normal",color=INK,ha="center",va="center"):
        return ax.text(x,y,s,fontsize=size,fontweight=weight,color=color,ha=ha,va=va,linespacing=1.35)
    def box(x,y,w,h,fill="#F4F7FA",edge="#BCC7D1",lw=.65):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.25,rounding_size=1.2",
                                   facecolor=fill,edgecolor=edge,linewidth=lw))
    def arrow(x1,y1,x2,y2,color=GREY,style="-",lw=.8):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=7,
                                    linewidth=lw,color=color,linestyle=style,shrinkA=0,shrinkB=0))
    def panel(letter,title,y):
        txt(2,y,letter,10.5,"bold",ha="left")
        txt(9,y,title,8.4,"bold",ha="left")

    panel("a","Aligned regional inputs and two signed LD operators",4)
    box(3,11,39,17)
    txt(22.5,16,"Two-trait z scores",7.2,"bold")
    txt(22.5,22.5,r"$z^{(1)}, z^{(2)}$  |  $p=100$ variants",7.0)
    arrow(42.5,19.5,50.5,19.5)
    box(51,11,116,17,"#EEF5FA","#9DBBD0")
    txt(109,16,"Within-region, per-trait standardization",7.2,"bold")
    txt(109,23,r"$X_j=[\tilde z_{1j},\tilde z_{2j},|\tilde z_{1j}|,|\tilde z_{2j}|,"
        r"\mathrm{sign}(\tilde z_{1j}\tilde z_{2j})]$     $X\in\mathbb{R}^{p\times5}$",7.4)
    box(3,33,39,22)
    txt(22.5,38,"Aligned LD matrix R",7.4,"bold")
    txt(22.5,45,"Same variant order\nSet diagonal to zero",7.0)
    arrow(42.5,44,49,44)
    # Shared split into separate positive and negative weight channels.
    ax.plot([48,48,111,111],[44,30.5,30.5,44],color=GREY,lw=.65)
    arrow(48,44,51,44,color=BLUE)
    arrow(111,44,115,44,color=ORANGE,style="--")
    box(51,33,54,22,"#EEF5FA","#8BB4CE")
    txt(78,37.5,"Positive LD channel",7.4,"bold",BLUE)
    txt(78,44,r"$W_+=\max(R_{\rm off},0)$",7.5)
    txt(78,50.5,r"$A_+=D_+^{-1/2}(W_++I)D_+^{-1/2}$",7.1)
    box(115,33,52,22,"#FCF3EB","#D0A47D")
    txt(141,37.5,"Negative LD channel",7.4,"bold",ORANGE)
    txt(141,44,r"$W_-=\max(-R_{\rm off},0)$",7.5)
    txt(141,50.5,r"$A_-=D_-^{-1/2}(W_-+I)D_-^{-1/2}$",7.1)
    txt(85,58,"Both channels include self-loops; D is the corresponding row-sum degree matrix.",6.8,color=GREY)

    panel("b","Two message-passing blocks (same operators, different learned weights)",65)
    def block(x,y,w,label,din):
        box(x,y,w,41,"#FFFFFF","#9CAEBB",.8)
        txt(x+w/2,y+4,label,7.5,"bold")
        slot=(w-8)/3
        labels=[r"$H_k$",r"$A_+H_k$",r"$A_-H_k$"]
        colors=[GREY,BLUE,ORANGE]
        for j in range(3):
            sx=x+2+j*(slot+2)
            box(sx,y+8,slot,9,["#F3F5F7","#EEF5FA","#FCF3EB"][j],colors[j])
            txt(sx+slot/2,y+12.5,labels[j],8,color=colors[j])
            arrow(sx+slot/2,y+17.5,sx+slot/2,y+20,colors[j],"--" if j==2 else "-")
        txt(x+w/2,y+22,f"Concatenate: p × {3*din}",7.1)
        txt(x+w/2,y+28,f"Linear {3*din} → 48  |  LayerNorm",7.1)
        txt(x+w/2,y+34,"GELU → Dropout (0.20)",7.1)
        txt(x+w/2,y+39,(r"$H_1$: $p\times48$" if din==5 else r"$H_2$: $p\times48$"),7.0,"bold")
    box(3,85,20,15,"#EEF5FA","#9DBBD0")
    txt(13,90,r"$H_0=X$",8.0,"bold")
    txt(13,96,r"$p\times5$",7.4)
    arrow(23.5,92.5,29,92.5)
    block(30,73,62,"Block 1   (k = 0)",5)
    arrow(92.5,92.5,100,92.5)
    block(101,73,66,"Block 2   (k = 1)",48)
    txt(85,119,"12,772 trainable parameters per network; dropout is disabled at inference.",7.0,color=GREY)

    panel("c","Graph readout, classification and three-model ensemble",126)
    widths=[36,43,70];xs=[3,47,97]
    for x,w in zip(xs,widths):
        box(x,133,w,22,"#F1F7F5" if x==3 else "#F4F7FA","#A4BBB4" if x==3 else "#BCC7D1")
    txt(21,138,"Mean + max pooling",7.2,"bold",TEAL)
    txt(21,146,r"Across nodes of $H_2$",6.8)
    txt(21,151,r"$p\times48\ \longrightarrow\ 96$",7.4)
    arrow(39.5,144,46,144)
    txt(68.5,137.5,"Classification head",7.3,"bold")
    txt(68.5,144,"96 → 48 → GELU → Dropout",6.5)
    txt(68.5,150.5,"48 → 4 logits → Softmax",7.0)
    arrow(90.5,144,96.5,144)
    txt(132,138,"Average three probability vectors",7.1,"bold")
    txt(132,144.5,"Independent seed offsets: 11, 29, 47",6.7)
    txt(132,151,r"$s_{\rm anomaly}=1-\bar P(\mathrm{valid})$",8.3,"bold")
    txt(85,160,"Four training states: valid | value-order fault | sign corruption | LD mismatch",7.1)
    box(3,165,164,8,"#FFF8EC","#D6B574",.7)
    txt(85,169,"Warning only: no automatic repair, no eligibility override, no biological inference.",7.1,"bold",color="#75500D")

    for ext in ("pdf","svg","png"):
        fig.savefig(OUT/f"Fig9_LD_GuardNet_architecture.{ext}",dpi=600,facecolor="white")
    plt.close(fig)

audit={"source_file":SOURCE.name,"source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
       "type":"deterministic code-grounded architecture schematic; no performance data",
       "parameter_count":12772,"tensor_shapes":shapes,"operator_formula_check":"passed",
       "width_mm":170,"height_mm":176,"png_dpi":600,"model_or_result_changes":False,
       "journal_guidance":"https://link.springer.com/journal/12864/submission-guidelines",
       "figure_formats":["vector PDF with embedded TrueType fonts","editable SVG text","600 dpi PNG"],
       "files":[]}
for p in sorted(OUT.glob("Fig9*")):
    audit["files"].append({"file":p.name,"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
(OUT/"architecture_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
print(json.dumps(audit,indent=2))
