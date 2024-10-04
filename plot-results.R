library(magrittr)
library(data.table)
## library(ggpubr)

## * null data - with spike=0
null_prk=readRDS("/rds/project/rds-csoP2nj6Y6Y/wyl37/Rbookdown/GSEApublicdat/nextflowResults/fgseaprerank/GSE81259.logFC.prerankGSEA.rds")
## prerank GSEA results based on true (pheno0) and 100 randomised data (pheno1-pheno100).
## A list of pheno0 - pheno 100. for each pheno, there are "generanks" (input for prerank), "bgpah" (pathways), "gseares" (fgsea::fgsea results)

null_lbl=readRDS("/rds/project/rds-csoP2nj6Y6Y/wyl37/Rbookdown/GSEApublicdat/nextflowResults/flgsea/GSE81259.logFC.flGSEA.rds")

null=lapply(seq_along(null_prk), function(i) {
    list(spiked=character(0),
         lbl=null_lbl[[i]][[3]],
         prk=null_prk[[i]][[3]])})

## * real data
real_bigsea=readRDS("/rds/project/rds-csoP2nj6Y6Y/wyl37/Rbookdown/GSEApublicdat/nextflowResults/BIGSEA/GSE81259.logFC.exprGSEA.rds")[[2]]
real_bigsea[,padj:=p.adjust(pval,method="BH")]


## * results of spike in expt
proj_dir <- "/rds/project/rds-csoP2nj6Y6Y/wyl37/Rbookdown/GSEApublicdat"
source(file.path(proj_dir, "r_functions.R"))
spike100 <- file.path(
proj_dir,
"nextflowResults/SPIKE100/GSE81259.pheno26.spikein_realexpr.Run.rds"
) %>%
    readRDS()

## list of 200
## first 100 = one spike, second 100 = two spikes
## each element is a list
## spikepath is the path spiked in
## fl_gseares - results using permuted gsea - data.table with 50 obs, including pval, padj. one row per pathway
## fprerank_res - results using preranked gsea - data.table as above


## * read in pathways
library(fgsea)
gmtin <- "/rds/project/cew54/rds-cew54-wallace-share/People/Chris/Elena/pathways/h.all.v2023.1.Hs.symbols.gmt"
gmt <- gmtPathways(gmtin)
## gmt is a list of 50 pathways, each pathway is a character vector of gene symbols

library(FastJaccard)
o=overlap_symlist(gmt)

## library(pheatmap)
## pheatmap(o) # mostly small, one example 100, another over 50
j=jaccard_symlist(gmt)
## pheatmap(j)
hist(as.vector(j[lower.tri(j)]),breaks=100)
hist(as.vector(o[lower.tri(o)]))

jv=as.vector(j[lower.tri(j)])
summary(jv)
quantile(jv,seq(0.8,1,by=0.01))
length(jv)
sum(jv > 0.1)
sum(jv > 0.05)
sum(jv > 0.01)

## 95% of pairs have JI 0.453 or less. a few higher, with 14/1225 with JI > 0.1 and 52 with JI > 0.5.

## try evaluating but discarding where JI  > 0.1
diag(j)=0
drop=which(j>0.05,arr.ind=TRUE)
dropc=cbind(names(gmt)[drop[,1]], names(gmt)[drop[,2]])

## evaluate p values for spiked
## library(qvalue)
f=function(l) {
    spiked=l[[1]]
    lbl=l[[2]]
    prk=l[[3]]
    w=which(dropc[,1] %in% spiked)
    if(length(w)) {
        lbl=lbl[ !(pathway %in% dropc[w,2]) ]
        prk=prk[ !(pathway %in% dropc[w,2]) ]
    }
    lbl[,padj:=p.adjust(pval,method="BH")]
    prk[,padj:=p.adjust(pval,method="BH")]
    dt=merge(prk[,.(pathway,pval.prk=pval,padj.prk=padj)], lbl[,.(pathway,pval.lbl=pval,padj.lbl=padj)], by="pathway")
    dt[,spike:=pathway %in% spiked]
    dt[,design:=paste0("spike_",length(spiked))]
    ## ## store padj at the spike first
    ## dt=data.table(spike=FALSE,
    ##            lbl=c(lbl[pathway %in% spiked]$padj, lbl[!(pathway %in% spiked)]$padj),
    ##            prk=c(prk[pathway %in% spiked]$padj, prk[!(pathway %in% spiked)]$padj))
    ## dt[1:length(spiked),spike:=TRUE]
    dt
}
sresults=lapply(spike100, f)  %>% rbindlist()
nresults=lapply(null,f)
nresults[[1]]$design="observed"
nresults %<>% rbindlist()
oresults=nresults[design=="observed"]

results=rbind(sresults, nresults)

results[,.(lbl=mean(pval.lbl<0.05),prk=mean(pval.prk<0.05, na.rm=TRUE)), by=c("spike","design")]

m=melt(results, c("pathway","spike","design"))
## order pathway by observed prk results
oresults=oresults[order(padj.prk)]
m[,pathway:=factor(pathway,levels=oresults$pathway)]

summ=m[,.(mn=mean(value<0.05,na.rm=TRUE),sd=sd(value<0.05,na.rm=TRUE),n=.N),by=c("pathway","spike","design","variable")]
summ[,lci:=mn-1.96*(sd/sqrt(n))][,uci:=mn+1.96*(sd/sqrt(n))]


## overall summary
osumm=m[,.(mn=mean(value<0.05,na.rm=TRUE),sd=sd(value<0.05,na.rm=TRUE),n=.N),by=c("spike","design","variable")]
osumm[,lci:=mn-1.96*(sd/sqrt(n))][,uci:=mn+1.96*(sd/sqrt(n))]
osumm=osumm[variable %in% c("pval.prk","pval.lbl") & design!="observed"]
osumm[variable=="pval.prk",variable:="gene label permutation (fGSEA)"]
osumm[variable=="pval.lbl",variable:="outcome permutation (GSEA)"]
## osumm[,spike:=ifelse(spike==FALSE,"Type 1 error rate","Power")]

theme_set(theme_cowplot(font_size=18))
p1=ggplot(osumm[spike==FALSE & design %in% c("spike_0","spike_2")], aes(x=variable,y=mn,ymin=lci,ymax=uci,col=variable)) +
    geom_pointrange(size=1) +
    background_grid() + 
    geom_hline(yintercept=0.05) +
    scale_x_discrete("analysis",labels=c("fGSEA","GSEA")) +
    scale_y_continuous(limits=c(0,1)) +
    scale_colour_discrete("analysis") +
    facet_grid(. ~ design) +
    theme(legend.position="bottom",axis.title.x=element_blank(),
          strip.background = element_rect(color="black", fill="#ffffff", size=1.5, linetype="solid")) +
    labs(x="analysis",y="value") +
    ggtitle("Type 1 error rate") 
p2=ggplot(osumm[spike==TRUE & design %in% c("spike_0","spike_2")], aes(x=variable,y=mn,ymin=lci,ymax=uci,col=variable)) +
    geom_pointrange(size=1) +
    background_grid() + 
    geom_hline(yintercept=0.05) +
    scale_x_discrete("analysis",labels=c("fGSEA","GSEA")) +
    scale_y_continuous(limits=c(0,1)) +
    scale_colour_discrete("analysis") +
    facet_grid(. ~ design) +
    theme(legend.position="none",axis.title.x=element_blank(),
          strip.background = element_rect(color="black", fill="#ffffff", size=1.5, linetype="solid")) +
    labs(x="analysis",y="value") +
    ggtitle("Power")

library(patchwork)
p1 + p2 + plot_layout(widths=c(2,1),axes="collect")

ggsave("summary_type1_power.png",height=6,width=10,bg="white")



## * read in raw data 
nulldir <- "/rds/project/rds-csoP2nj6Y6Y/wyl37/Rbookdown/GSEApublicdat/nextflowResults/dat4gsea"
    
    ## real expression data with permuated phenotype
    ## 1. objinput (DGEList); read count, gene info, 
    ## 2. qlf (MArrayLM class): DGE analysis based on permutated phenotype
dsnin <- "GSE81259.pheno26.rds"
dsn1 <- readRDS(file.path(nulldir, dsnin))

rc_mat <- dsn1$objinput
identical(rownames(dsn1$objinput),dsn1$objinput@.Data[[3]]$Geneid)
use=!duplicated(dsn1$objinput@.Data[[3]]$gene_name)
rc_mat=rc_mat[use,]
rownames(rc_mat)=dsn1$objinput@.Data[[3]]$gene_name[use]
rc_mat  %<>% as.matrix()  %>% t()

## correlation within pathways
corr=lapply(gmt, function(x) {
    xx=intersect(x,colnames(rc_mat))
    r=cor(rc_mat[,xx])
    r[upper.tri(r)]
})
names(corr)=names(gmt)

corr.dt=data.table(corr=unlist(corr),
                   pathway=rep(names(corr),times=sapply(corr,length)))
medcorr=corr.dt[,.(med=median(corr)),by="pathway"][order(med)]
corr.dt=merge(corr.dt, oresults[,.(pathway,sig=padj.prk<0.05)], by="pathway")
corr.dt[,pathway:=factor(pathway, levels=oresults$pathway)]
summ=merge(summ,oresults[,.(pathway,sig=padj.prk<0.05)], by="pathway")

tmp=summ[design=="spike_0" & variable=="pval.prk"][order(mn)]

corr.dt[,pathway:=as.character(pathway)][,pathway:=factor(pathway,levels=as.character(tmp$pathway))]
corr.dt[,variable:="correlation"]
summ[,pathway:=as.character(pathway)][,pathway:=factor(pathway,levels=as.character(tmp$pathway))]
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
pcorr=ggplot(corr.dt, aes(x=corr^2,y=pathway, fill=sig)) + geom_boxplot() + background_grid() + geom_vline(xintercept=0,lty="dashed",colour="grey") + facet_grid(sig ~ variable, scales="free_y")

m=merge(m, oresults[,.(pathway,sig=padj.prk<0.05)], by="pathway")
pp=ggplot(summ[spike==FALSE & design!="observed" & grepl("pval",variable)],
       aes(x=mn,y=pathway,col=variable)) +
    geom_pointrange(aes(xmin=lci,xmax=uci)) +
    geom_point(aes(x=value),data=m[design=="observed" & grepl("pval",variable)]) +
    facet_grid(sig ~ design, scales="free_y") +
    background_grid(major="y") +
    geom_vline(xintercept=0.05,linetype="dashed",col="grey")

plot_grid(pp + theme(legend.position="top"),
          pcorr + theme(legend.position="top",
                        axis.text.y=element_blank(),
                        axis.title.y=element_blank()),
          rel_widths=c(3,1))

# high correlation, low p
a=rc_mat[,intersect(colnames(rc_mat),gmt$HALLMARK_INTERFERON_ALPHA_RESPONSE)]
# mid correlation, low p
b=rc_mat[,intersect(colnames(rc_mat),gmt$HALLMARK_TNFA_SIGNALING_VIA_NFKB)]

## look at predictors
cs=corr.dt[,.(med=median(corr),mx=max(abs(corr)),q99=quantile(corr,.99),q95=quantile(corr,.95),q90=quantile(corr,.9),n=.N),by="pathway"]
cs=merge(cs,tmp,by="pathway")
lm(mn ~ mx + n.x, data=cs)  %>% summary()

tmp=summ[design=="spike_0" & variable=="pval.prk"][order(mn)]  %>%
    merge(., 
