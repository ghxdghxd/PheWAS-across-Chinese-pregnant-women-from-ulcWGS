# Figure 3

## Figure 3A

```R
library(data.table)
library(stringr)
load("phewas_power_P0.001.RData")
res_sig_raw = snp_power[which(snp_power$power > 0.8), ]

bed_raw = unique(res_sig_raw$MarkerID)
bed_raw = data.frame(SNP=bed_raw, stringr::str_split(bed_raw, pattern=":", simplify=T)[,1:2], stringsAsFactors = F)
colnames(bed_raw) = c("SNP", "CHR", "BP")
bed_raw$CHR = as.numeric(bed_raw$CHR)
bed_raw$BP = as.numeric(bed_raw$BP)

res_sig_raw = cbind(res_sig_raw, bed_raw[match(res_sig_raw$MarkerID, bed_raw$SNP), c("CHR", "BP")])
res_sig_raw_bed = res_sig_raw[,c("CHR","BP")]
res_sig_raw_bed$start = res_sig_raw_bed$BP - 1
colnames(res_sig_raw_bed) = c("chrom", "end", "start")
res_sig_raw_bed = res_sig_raw_bed[, c("chrom", "start", "end")]

#gene Symbol
gene_bed = read.table("hg38.RefSeqGene.bed", header=F,sep="\t",stringsAsFactors=F)
colnames(gene_bed) = c("chrom", "start", "end", "gene")
gene_bed$chrom = gsub("chr", "", gene_bed$chrom)

res_sig_raw_bed_gene = unique(bedtoolsr::bt.intersect(res_sig_raw_bed, gene_bed, wb=T)[,c(1,3,7)])

res_sig_raw = merge(res_sig_raw, res_sig_raw_bed_gene, by.x=c("CHR","BP"), by.y=c("V1","V3"), all.x=T)
colnames(res_sig_raw)[ncol(res_sig_raw)] = "gene"
res_sig_raw$gene = as.character(res_sig_raw$gene)

# cytoband
cytoband_bed = read.table("cytoband_hg38.txt", header=F,sep="\t",stringsAsFactors=F)
colnames(cytoband_bed) = c("chrom", "start", "end", "cytoband", "V5")

res_sig_raw_bed_cytoband = unique(bedtoolsr::bt.intersect(res_sig_raw_bed, cytoband_bed, wb=T)[,c(1,3,7)])

res_sig_raw = merge(res_sig_raw, res_sig_raw_bed_cytoband, by.x=c("CHR","BP"), by.y=c("V1","V3"), all.x=T)
colnames(res_sig_raw)[ncol(res_sig_raw)] = "cytoband"
res_sig_raw$cytoband = as.character(res_sig_raw$cytoband)

res_sig_raw$ICDanno = ""
res_sig_raw$ICDanno[grep("[AB]", res_sig_raw$icd)] = "Certain infectious and parasitic diseases"
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24", "D25", "D26"))] = "Neoplasms"
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
res_sig_raw$ICDanno[grep("E", res_sig_raw$icd)] = "Endocrine, nutritional and metabolic diseases"
res_sig_raw$ICDanno[grep("F", res_sig_raw$icd)] = "Mental and behavioural disorders"
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("H00", "H01", "H02", "H04", "H10","H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92","H93"))] = "Ear and mastoid process"
res_sig_raw$ICDanno[grep("I", res_sig_raw$icd)] = "Circulatory system"
res_sig_raw$ICDanno[grep("J", res_sig_raw$icd)] = "Respiratory system"
res_sig_raw$ICDanno[grep("K", res_sig_raw$icd)] = "Digestive system"
res_sig_raw$ICDanno[grep("L", res_sig_raw$icd)] = "Skin and subcutaneous tissue"
res_sig_raw$ICDanno[grep("M", res_sig_raw$icd)] = "Musculoskeletal system and connective tissue"
res_sig_raw$ICDanno[grep("N", res_sig_raw$icd)] = "Genitourinary system"
res_sig_raw$ICDanno[grep("O", res_sig_raw$icd)] = "Pregnancy, childbirth and the puerperium"
res_sig_raw$ICDanno[grep("P", res_sig_raw$icd)] = "Certain conditions originating in the perinatal period"
res_sig_raw$ICDanno[grep("Q", res_sig_raw$icd)] = "Congenital malformations, deformations and chromosomal abnormalities"

load("LDblock_1e3_Rsq0.8.RData")
res_sig_raw$LDblock = ldblock$L1[match(res_sig_raw$MarkerID, ldblock$value)]

res_sig_raw_anno = fread("pheWAS_sig_SAIGE_1e3_power0.8.lite.hg38_multianno.txt", header=T,sep="\t",stringsAsFactors=F)
res_sig_raw_anno$ID = paste(res_sig_raw_anno$Otherinfo4, res_sig_raw_anno$Otherinfo5, res_sig_raw_anno$Otherinfo7, res_sig_raw_anno$Otherinfo8, sep=":")

res_sig_raw$Func_refGene = res_sig_raw_anno$Func.refGene[match(res_sig_raw$MarkerID, res_sig_raw_anno$ID)]
save(res_sig_raw, file="pheWAS_sig_SAIGE_1e3_power0.8.RData")

pdf("Figure3A.pdf", width=15, height=6)
plot_all_icd_snp(res_sig_raw[res_sig_raw$p.value<1e-5, ], title = "res_sig_raw_manhantan_ld_1e5", ldblock = T)
dev.off()
```

## Figure 3B

```R
library(UpSetR)
dat = as.data.frame.matrix(table(unique(res_sig[,c("LDblock", "ICDanno")])))
colnames(dat) = str_split(colnames(dat), pattern=" ", simplify=T)[,1]

pdf("Figure3B.pdf", width= 6, heigh = 6)
# upset(dat, sets = colnames(dat)[-1], order.by = "freq")
upset(dat, sets = names(sort(colSums(dat))), nintersects = 20,
    order.by = "freq", keep.order = T, point.size = 2, line.size = 0.8, mb.ratio = c(0.6, 0.4),
    mainbar.y.label = "Intersection Size", sets.x.label = "Set size")
aplot::plot_list(
    reshape2::melt(table(apply(dat, 1, function(x){paste(colnames(dat)[x>0], collapse=",")}))) %>% arrange(value) %>%
            mutate(Var1 = factor(Var1, levels = rev(Var1))) %>%
            ggplot() + geom_col(aes(Var1, value)) +
            geom_text(aes(Var1, value, label = value), hjust = 0, nudge_y = -0.2, size = 2.5) +
            labs(y = "Intersection Size", x = "") +
            scale_y_break(c(100, 260), scales = 0.3, ticklabels = c(270,320,350), expand = TRUE) +
            scale_y_sqrt(expand = c(0,0)) +
            theme_classic() +
            theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()),
    reshape2::melt(colSums(dat)) %>%
        tibble::rownames_to_column() %>% arrange(value) %>%
        mutate(rowname = factor(rowname, levels = rowname)) %>%
        ggplot() + geom_col(aes(value, rowname)) +
        geom_text(aes(value, rowname, label = value), hjust = 1, nudge_x = 0, size = 2.5) +
        labs(x = "Set Size", y = "") +
        scale_x_break(c(20, 90), scales = 0.2, ticklabels = 95, expand = TRUE) +
        scale_x_break(c(100, 450), scales = 0.3, ticklabels = c(460,800,1000), expand = TRUE) +
        scale_x_sqrt(expand = c(0,0)) +
        theme_classic() +
        theme(axis.ticks.y = element_blank()),
    nrow = 2, heights = c(0.6, 0.4))
dev.off()
```

## Figure 3C

```R
library(ggbreak)
library(dplyr)

plot_snp_levels = function(res_sig){
    ndat_snp_levels = reshape2::melt(table(unique(res_sig[,c("MarkerID","icd","levels")])[,2:3]))
    icds = ndat_snp_levels %>% group_by(icd) %>% summarise(n = sum(value)) %>% arrange(desc(n)) %>% select(icd) %>% unlist %>% as.character
    ndat_snp_levels = ndat_snp_levels[ndat_snp_levels$icd %in% icds, ]
    ndat_snp_levels$icd = factor(ndat_snp_levels$icd, levels = icds)
    # ndat_snp_levels$levels = factor(ndat_snp_levels$levels, levels =c("L4", "L3", "L2", "L1"))
    p_snp_levels = ggplot(ndat_snp_levels) + geom_col(aes(icd, value, fill=levels)) +
        theme_bw() + labs(x="",y="Associations pairs' count") +
        scale_fill_npg() +
        scale_y_sqrt() +
        scale_y_break(c(100,150), scale='free') +
        theme(axis.text.x = element_text(size=10, angle=-90, vjust=0.5, hjust=1))
    return(p_snp_levels)
}

res_sig_ld = res_sig %>% group_by(icd, LDblock) %>% mutate(lab = MarkerID[which.min(p.value)], gwas = any(gwas), AS = any(AS),
    EAS = any(EAS), SAS = any(SAS), CHB = any(CHB), eQTL = any(eQTL), phewas = any(phewas),labwas = any(labwas),geneatlas = any(geneatlas)) %>% filter(MarkerID==lab) %>% mutate(lab=NULL) %>% as.data.frame
p_snp_ld_levels = plot_snp_levels(res_sig_ld)

pdf("pheWAS_snp_icd_count.levels.1.pdf", width=8, height=4)
p_snp_ld_levels
dev.off()
```
