# Figure 5 and S7

## Figure 5A

```R
load("ICD_FF_status_num10.RData")
ICD_mat_status = ICD_mat_status[, c(colnames(ICD_mat_status)[1:4], names(which(colSums(ICD_mat_status[,-c(1:4)], na.rm = T)>30)))]

ICD_res_sig = read.table("ICD_res_sig.txt", header=T,sep="\t",stringsAsFactors=F)

ldblock_1e3_Rsq8 = lapply("LDblock_1e3_Rsq0.8.RData", function(x){
    load(x)
    return(ldblock)
})
ldblock_1e3_Rsq8 = ldblock_1e3_Rsq8[[1]]
ICD_res_sig$LDblock = ldblock_1e3_Rsq8$L1[match(ICD_res_sig$snp, ldblock_1e3_Rsq8$value)]

load("phewas_power_P0.001.RData")
res_sig_raw = snp_power[which(snp_power$power > 0.8), ]

MC_ld = merge(unique(res_sig_raw[,c("MarkerID","icd", "LDblock")]), unique(ICD_res_sig[,c("snp","icd", "pvalue_norm", "LDblock")]), by.x="LDblock", by.y = "LDblock")

MC_ld_hc_C = hclust(dist(table(unique(MC_ld[,c("icd.y","icd.x")])), method="binary"), method="ward.D2")
MC_ld_hc_M = hclust(dist(table(unique(MC_ld[,c("icd.x","icd.y")])), method="binary"), method="ward.D2")

MC_ld_hc_C_d <- as.dendrogram(MC_ld_hc_C)
MC_ld_hc_M_d <- as.dendrogram(MC_ld_hc_M)

MC_ld5_hc_C = hclust(dist(table(unique(MC_ld5[,c("icd.y","icd.x")])), method="binary"), method="ward.D2")
MC_ld5_hc_M = hclust(dist(table(unique(MC_ld5[,c("icd.x","icd.y")])), method="binary"), method="ward.D2")

MC_ld5_hc_C_d <- as.dendrogram(MC_ld5_hc_C)
MC_ld5_hc_M_d <- as.dendrogram(MC_ld5_hc_M)

# Rectangular lines
library(ggdendro)

MC_ld_dat = reshape2::melt(table(MC_ld[,c("icd.x","icd.y")]))
MC_ld_dat = merge(MC_ld_dat, MC_ld, all.x=T)
MC_ld_dat$P = -log10(MC_ld_dat$pvalue_norm)
MC_ld_dat = MC_ld_dat %>% group_by(icd.x, icd.y, value, LDblock) %>% summarise(pvalue = min(pvalue_norm), P = max(P))
MC_ld_dat$icd.x = as.character(MC_ld_dat$icd.x)
MC_ld_dat$icd.y = as.character(MC_ld_dat$icd.y)
MC_ld_dat_1 = merge(MC_ld_dat, reshape2::melt(table(sort(unique(c(MC_ld_dat$icd.x, MC_ld_dat$icd.y))), sort(unique(c(MC_ld_dat$icd.x, MC_ld_dat$icd.y)))))[,-3], by.x = c("icd.x", "icd.y"), by.y = c("Var1", "Var2"), all.y=T)

MC_ld_dat$icd.x = factor(MC_ld_dat$icd.x, levels = MC_ld_hc_M$labels[MC_ld_hc_M$order])
MC_ld_dat$icd.y = factor(MC_ld_dat$icd.y, levels = MC_ld_hc_C$labels[MC_ld_hc_C$order])

MC_GT_pvalue_files_hist = list.files("MC_ICD_barplot", pattern="*[012].hist.txt", full.names = T)
MC_GT_pvalue_hist = do.call(rbind, lapply(MC_GT_pvalue_files_hist, function(x){
    a = read.table(x, header=T,sep="\t",stringsAsFactors =F)
    a$file=gsub(".txt","",basename(x))
    return(a)
}))
MC_GT_pvalue_hist$ICD_1 = gsub("[\\+\\-]", "", str_split(MC_GT_pvalue_hist$ICD, pattern=",", simplify=T)[,2])
MC_GT_pvalue_hist$ICD_1[MC_GT_pvalue_hist$ICD_1 == ""] = MC_GT_pvalue_hist$ICD[MC_GT_pvalue_hist$ICD_1 == ""]

ICD_icd_paired_01_hist = MC_GT_pvalue_hist %>% filter(p<0.05) %>% select(ICD_1, icd.icd.y) %>% unique
ICD_icd_paired_01_hist = paste(ICD_icd_paired_01_hist$ICD_1, ICD_icd_paired_01_hist$icd.icd.y)

MC_ld_dat$lab_hist = ""
MC_ld_dat$lab_hist[paste(MC_ld_dat$icd.x, MC_ld_dat$icd.y) %in% c(ICD_icd_paired_01_hist)] = "*"

p_MC_ld_0 = ggplot(MC_ld_dat) + geom_tile(aes(icd.x, icd.y, fill = P), colour = "black") +
        geom_text(aes(icd.x, icd.y, label = lab_hist), hjust = 0.5, vjust = 0.8, size=5) +
        labs(x="Maternal",y="Children", fill="-log10(pvalue)") +
        scale_fill_distiller(palette = "Spectral",na.value="white", limits = c(1.3,1.6)) +
        coord_fixed(ratio = 1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
            axis.text.y.right = element_text(),
            legend.position = "right",
            legend.key.width=unit(0.5,"cm"),
            legend.key.height=unit(1,"cm"),
            plot.background = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

pdf("Figure5A.pdf", width = 10, height = 10)
p_MC_ld_0
dev.off()
```

## Figure 5B and S7

```R
MC_GT_plot_hist = function(snp, icd, ICD, GTgroup = NULL, plot = T){
    ref = strsplit(snp, split=":")[[1]][3]
    alt = strsplit(snp, split=":")[[1]][4]
    if(is.null(GTgroup)){
        GTlab = c(paste(ref, ref, sep=":"), paste(ref, alt,sep=":"), paste(alt, alt,sep=":"))
    }else if(GTgroup == "01"){
        GTlab = c(paste(paste(ref, ref, sep=":"), paste(ref, alt,sep=":"), sep="/"), paste(paste(ref, ref, sep=":"), paste(ref, alt,sep=":"), sep="/"), paste(alt, alt,sep=":"))
    }else if(GTgroup == "12"){
        GTlab = c(paste(ref, ref,sep=":"), paste(paste(ref, alt,sep=":"), paste(alt, alt,sep=":"), sep="/"), paste(paste(ref, alt,sep=":"), paste(alt, alt,sep=":"), sep="/"))
    }
    int = intersect(rownames(clin_status), ICD_mat_status$ID)
    dat = child_ICD_NIPT[which(child_ICD_NIPT$ID %in% int),c("ID", "follow_year", icd)]
    colnames(dat)[3] = 'icd'
    dat$ICD = ifelse(clin_status[dat$ID, ICD]==1, paste0(ICD, "+"), paste0(ICD, "-"))
    dat$ICD[dat$ICD==paste0(ICD, "+")] = paste(as.numeric(table(dat$ICD)[paste0(ICD, "+")]), dat$ICD[dat$ICD==paste0(ICD, "+")], sep=",")
    dat$ICD[dat$ICD==paste0(ICD, "-")] = paste(as.numeric(table(dat$ICD)[paste0(ICD, "-")]), dat$ICD[dat$ICD==paste0(ICD, "-")], sep=",")
    dat$ICD = factor(dat$ICD, levels = c(grep("-", dat$ICD,v=T)[1], grep("-", dat$ICD,v=T, inver=T)[1]))
    dat$GT = geno[dat$ID, snp]

    num = reshape2::melt(table(dat[,c("ICD","GT")]))
    num$labs = GTlab[num$GT + 1]
    # num$labs = paste(num$labs, num$value, sep="\n")
    dat = merge(dat, num[,c("ICD","GT","labs")])
    dat$GT = factor(geno[dat$ID, snp], levels = c(0, 1, 2), labels = GTlab)
    dat$icd1 = ifelse(dat$icd>0, paste0(icd, "+"), paste0(icd, "-"))
    # dat$icd1 = factor(dat$icd1, levels = c(paste0(icd, "+"), paste0(icd, "-")))
    a = as.data.frame(do.call(rbind, lapply(unique(as.character(dat$ICD)), function(x){
            return(c(ICD = x, icd = icd, snp = snp, p = fisher.test(table(dat[dat$ICD==x, c("icd1", "GT")]))$p.value))
        })),stringsAsFactors=F)
    if(!plot){
        return(as.data.frame(a, stringsAsFactors=F))
    }
    p = ggplot(dat, aes(GT, ..count../sum(..count..),label=..count.., fill = icd1)) + geom_bar(position = "fill") +
        geom_text(stat = 'count', position='fill', vjust = 1.5, size = 5) +
        labs(y=paste(icd, "Prevalence"), x=snp, fill="") +
        facet_grid(~paste("P =", signif(fisher.test(table(dat[, c("icd1", "GT")]))$p.value,3))) +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_d3() +
        theme_bw() + theme(legend.position="none",aspect.ratio=1)
    dat$ICD1 = paste0(dat$ICD,"\nP = ", signif(as.numeric(a$p[match(dat$ICD, a$ICD)]),3))
    p0 = ggplot(dat, aes(GT, ..count../sum(..count..),label=..count.., fill = icd1)) + geom_bar(position = "fill") +
        geom_text(stat = 'count', position='fill', vjust = 1.5, size = 5) +
        labs(y=paste(icd, "Prevalence"), x=snp, fill="") +
        facet_grid(~ICD1) +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_d3() +
        theme_bw() + theme(aspect.ratio=1)
    return(plot_grid(p, p0, ggtexttable(num %>% arrange(ICD, GT), rows = NULL, theme = ttheme("classic")),
        nrow = 1, rel_widths = c(1.1, 2, 1), axis="h", align="tb"))
}

library(ggpubr)
library(ggsci)
library(cowplot)
library(dplyr)
dat0 = merge(unique(res_sig_raw[,c("MarkerID","icd", "LDblock")]), unique(ICD_res_sig[,c("snp","icd")]), by.x="MarkerID", by.y = "snp")

plot_hist_list = apply(dat0, 1, function(x){
    MC_GT_plot_hist(x["MarkerID"], x["icd.y"], x["icd.x"])
})

pdf("Figure5B_S7.pdf", width = 12, heigh = 4)
plot_hist_list
dev.off()

MC_GT_hist_pvalue = apply(dat0, 1, function(x){
    MC_GT_plot_hist(x["MarkerID"], x["icd.y"], x["icd.x"], plot = F)
})
MC_GT_hist_pvalue = do.call(rbind, MC_GT_hist_pvalue)
write.table(MC_GT_hist_pvalue, file="ICD_icd01_snp012_hist.txt", r=F,c=T,sep="\t",quote=F)
```
