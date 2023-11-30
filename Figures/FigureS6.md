# Figure S6

```R
library(data.table)
load("phewas_power_P0.001.RData")
load("pheWAS_sig_SAIGE_5e8.RData")

load("gwas_catalog_mat.RData")
gwas_mat = gwas_mat[!is.na(gwas_mat$CHR_ID), ]
gwas_mat = gwas_mat[-which(gwas_mat$CHR_ID %in% c("X", "Y")), ]

get_gwas_overlap <- function(gwas_mat, res_sig, width){
    gwas_mat$chrom = paste0("chr", gwas_mat$CHR_ID)
    gwas_mat$start = gwas_mat$CHR_POS - width
    gwas_mat$start[gwas_mat$start<0] = 0
    gwas_mat$end = gwas_mat$CHR_POS + width
    gwas_mat$ID = gsub("chr", "", paste(gwas_mat$chrom, gwas_mat$start, gwas_mat$end, sep=":"))

    gwas_mat_bed = na.omit(unique(gwas_mat[,c("chrom","start","end","ID")]))

    res_sig_gwas = res_sig
    res_sig_gwas$CHR = as.character(res_sig_gwas$CHR)

    res_sig_gwas = res_sig_gwas[,c("CHR","BP")]
    res_sig_gwas$start = res_sig_gwas$BP - 1
    colnames(res_sig_gwas) = c("chrom", "end", "start")
    res_sig_gwas = res_sig_gwas[, c("chrom", "start", "end")]
    res_sig_gwas$chrom = paste0("chr", res_sig_gwas$chrom)
    res_sig_gwas = unique(bedtoolsr::bt.intersect(res_sig_gwas, gwas_mat_bed, wb=T,wa=T))
    res_sig_gwas$ID = gsub("chr", "", paste(res_sig_gwas$V1, res_sig_gwas$V3, sep=":"))
    res_sig_gwas$ID_gwas = gsub("chr", "", paste(res_sig_gwas$V4, res_sig_gwas$V5, res_sig_gwas$V6, sep=":"))
    res_sig_gwas = res_sig_gwas[,c("ID", "ID_gwas")]
    return(res_sig_gwas)
}

res_sig_gwas_50kb = get_gwas_overlap(gwas_mat, res_sig, width=50000)

library(stringr)
library(dplyr)
library(ggplot2)
library(ggsci)
get_hist_density_dat = function(res_sig_gwas_50kb){
    res_sig_gwas_50kb = cbind(res_sig_gwas_50kb, str_split(res_sig_gwas_50kb$ID, pattern=":", simplify=T))
    colnames(res_sig_gwas_50kb)[3:4] = c("ID_chr", "ID_pos")
    res_sig_gwas_50kb$ID_pos = as.numeric(as.character(res_sig_gwas_50kb$ID_pos))
    res_sig_gwas_50kb = cbind(res_sig_gwas_50kb, str_split(res_sig_gwas_50kb$ID_gwas, pattern=":", simplify=T))
    colnames(res_sig_gwas_50kb)[5:7] = c("ID_gwas_chr", "ID_gwas_start", "ID_gwas_end")
    res_sig_gwas_50kb$ID_gwas_start = as.numeric(as.character(res_sig_gwas_50kb$ID_gwas_start))
    res_sig_gwas_50kb$ID_gwas_end = as.numeric(as.character(res_sig_gwas_50kb$ID_gwas_end))
    res_sig_gwas_50kb$ID_gwas_pos = sapply(1:nrow(res_sig_gwas_50kb), function(x){median(c(res_sig_gwas_50kb$ID_gwas_start[x], res_sig_gwas_50kb$ID_gwas_end[x]))})
    res_sig_gwas_50kb$diff = res_sig_gwas_50kb$ID_pos - res_sig_gwas_50kb$ID_gwas_pos
    res_sig_gwas_50kb = res_sig_gwas_50kb %>% group_by(ID) %>% summarise(diff = diff[which.min(abs(diff))])
    return(res_sig_gwas_50kb)
}

res_sig_gwas_50kb_dat = get_hist_density_dat(res_sig_gwas_50kb)

res_sig_gwas_50kb_dat$type = "50kb"

plot_hist_density = function(x, breaks = 10){
    res = density(x)
    hist(x, breaks, prob=T, ylim = c(0, max(res$y)), main = "", xlab = "Position to known risk loci")
    lines(res,lwd=2)
}

pdf("FigureS6.pdf", width = 6, heigh = 6)
plot_hist_density(res_sig_gwas_50kb_dat$diff, breaks = 50)
dev.off()
```
