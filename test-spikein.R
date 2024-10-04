## * load data

nulldir <- "/rds/project/rds-csoP2nj6Y6Y/wyl37/Rbookdown/GSEApublicdat/nextflowResults/dat4gsea"

## real expression data with permuated phenotype
## 1. objinput (DGEList); read count, gene info, 
## 2. qlf (MArrayLM class): DGE analysis based on permutated phenotype
dsnin <- "GSE81259.pheno26.rds"
dsn1 <- readRDS(file.path(nulldir, dsnin))

## dsn1$qlf$design  %<>% rbind(.,.)
## dsn1$objinput$counts  %<>% cbind(.,.)
## dsn1$objinput$samples  %<>% rbind(.,.)


gmtin <- "/rds/project/cew54/rds-cew54-wallace-share/People/Chris/Elena/pathways/h.all.v2023.1.Hs.symbols.gmt"

n_cpus <- 4

source("/rds/project/rds-csoP2nj6Y6Y/wyl37/Rbookdown/GSEApublicdat/foorun.R")

## * look at spikein_realexpr

library(limma)
library(data.table)
library(fgsea)
library(magrittr)
pname=NULL; nrich=1; indep=FALSE; dropout=NULL
## DGE based on permuated phenotype
abc <- topTable(dsn1$qlf,
                coef = ncol(dsn1$qlf$design),
                n = Inf)
## prep_generank function from foorun.R
## 1. rankings: ranked gene based on rancol
## 2. Geneid: indicate which ones to keep
abc.r <- prep_generank(res = abc, rankcol = "logFC")

## read count in the order the pre-rank
rc_mat <- dsn1$objinput[abc.r$Geneid, ]
rownames(rc_mat) <- names(abc.r$rankings)
design <- dsn1$qlf$design
## rc_mat  %<>% cbind(., ., ., .) # double sample size
## design  %<>%  rbind(., ., ., .)

## voom normalisation
vv <- voom(counts = rc_mat, design = design)

## normalisation expression for exprGSEA
vE <- vv$E

## gene sets in as list
gmt <- gmtPathways(gmtin)
## gmt <- copy(pathways)

## pathways to be enriched
doplot=TRUE
get_path_p=function(doplot=FALSE) {
    path_p <- sample(names(gmt), nrich)

    ## genes in the target pathways and in the expression matrix
    opath_pg <- path_pg <- Reduce(f = union, x = gmt[path_p]) %>%
        intersect(., rownames(vE))

    ## find a subset of genes in path_pg that have limited correlation with the remainder
    path_r <- setdiff(rownames(vE), path_pg)
    cr=cor(t(vE)[,path_pg], t(vE)[,path_r])^2
    if(doplot)
        pheatmap(cr)
    mn=rowMeans(cr)
    keep=mn < quantile(mn,.5)
    path_pg=path_pg[keep]
    path_r <- setdiff(rownames(vE), path_pg)
    
    path_r <- setdiff(rownames(vE), path_pg)
    path_order <- c(path_pg,
                    sample(path_r, size = length(path_r)))

    vE <- vE[path_order, ] # log2 CPM

    ## assign beta
    bb <- c(runif(n = length(path_pg), min = 20, max = 50),
            runif(n = length(path_r), min = 0, max = 0))

    eee <- rnorm(ncol(vE), sd = 0.1)

    ## newy <- (bb %*% abs(vE) + eee) %>% as.vector
    newy <- (bb %*% vE + eee) %>% as.vector

    ## ## check: is newy correlated with pathway genes?
    if(doplot) {
        par(mfrow=c(2,1))
        hist(cor(newy, t(vE[opath_pg,])) %>% as.vector(), xlim=c(-1,1),main=path_p)
        hist(cor(newy, t(vE[setdiff(path_r,opath_pg),])) %>% as.vector(), xlim=c(-1,1))
    }
    
    ## dge based on newy
    design[, ncol(design)] <- newy

    vv <- voom(counts = rc_mat, design = design)

    ## normalisation expression for exprGSEA
    vE <- vv$E

    fit <- lmFit(vv, design)
    fit <- eBayes(fit)
    resin <- limma::topTable(fit,
                             coef = ncol(fit$design),
                             n = Inf)

    result <- fgseaLabel(pathways = gmt,
                         mat = vE,
                         labels = newy,
                         nperm = 10000, nproc = n_cpus )[order(pval)]
    w=which(result$pathway==path_p)
    data.table(pathway=path_p,
               raw_p=result$pval[w],
               adj_p=p.adjust(result$pval)[w])
}

pvals2=replicate(20,get_path_p(),simplify=FALSE)
rbindlist(pvals2)
rbindlist(pvals)
