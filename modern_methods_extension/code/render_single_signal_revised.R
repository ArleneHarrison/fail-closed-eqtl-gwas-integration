# Replot preserved consumer results; correct the historical fault label only.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(svglite)})
a <- commandArgs(TRUE); stopifnot(length(a)==2)
d <- read.delim(a[1]); dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
keys <- c("correct","summary_order_permuted","ten_percent_sign_flipped","EAS_LD_substituted","AFR_LD_substituted")
labels <- c("Correct input","Values permuted\n(IDs unchanged)","10% sign flips","EAS LD\nsubstituted","AFR LD\nsubstituted")
d$condition <- factor(d$condition,levels=keys,labels=labels)
panel <- function(y,label,tag,fill,ink) {
  ggplot(d,aes(x=condition,y=.data[[y]]))+
    geom_boxplot(width=.65,outlier.size=.7,fill=fill,colour=ink)+
    labs(x=NULL,y=label,tag=tag)+theme_classic(base_size=9)+
    theme(axis.text.x=element_text(angle=20,hjust=1),
          plot.tag=element_text(face="bold"),plot.margin=margin(8,8,8,12))
}
fig <- panel("pip_l1_drift","PIP L1 drift versus correct input","a","#BFD7EA","#17324D") /
       panel("causal_pip","PIP of simulated causal variant","b","#D9EAD3","#254117") /
       panel("max_pp_h4","Regional PP.H4 (coloc.abf)","c","#FCE5CD","#6B3A1E")
ggsave(file.path(a[2],"Fig6.png"),fig,width=6.69,height=7.5,dpi=600,bg="white")
ggsave(file.path(a[2],"Fig6.svg"),fig,width=6.69,height=7.5,device=svglite,bg="white")
