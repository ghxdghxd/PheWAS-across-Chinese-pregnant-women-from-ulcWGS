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
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24","D25", "D26"))] = "Neoplasms"
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
res_sig_raw$ICDanno[grep("E", res_sig_raw$icd)] = "Endocrine, nutritional and metabolic diseases"
res_sig_raw$ICDanno[grep("F", res_sig_raw$icd)] = "Mental and behavioural disorders"
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("H00", "H01", "H02", "H04", "H10", "H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
res_sig_raw$ICDanno[which(res_sig_raw$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92" ))] = "Ear and mastoid process"
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

plot_all_icd_snp = function(res_sig, title, ldblock = F){
    if(ldblock){
        res_sig = res_sig %>% group_by(ICDanno, LDblock) %>%
            summarise(CHR = CHR[which.min(p.value)], BP = BP[which.min(p.value)], OR = OR[which.min(p.value)], BETA = BETA[which.min(p.value)], P = min(p.value), cytoband = cytoband[which.min(p.value)]) %>% as.data.frame(stringsAsFactors = F)
    }else{
        res_sig = res_sig %>% mutate(P = p.value, OR_mean = OR)
    }
    res_sig$OR1 = res_sig$OR
    res_sig$OR1[res_sig$OR >= 15] = ">=15"
    res_sig$OR1[res_sig$OR >= 10 & res_sig$OR < 15] = "[10,15)"
    res_sig$OR1[res_sig$OR >= 5 & res_sig$OR < 10] = "[5,10)"
    res_sig$OR1[res_sig$OR >= 1 & res_sig$OR < 5] = "[1,5)"
    res_sig$OR1[res_sig$OR > 0 & res_sig$OR < 1] = "(0,1)"
    don <- res_sig %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len=max(as.numeric(BP))) %>% 
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        # Add this info to the initial dataset
        left_join(res_sig, ., by=c("CHR"="CHR")) %>%
        # Add a cumulative position of each SNP
        arrange(CHR, BP) %>%
        mutate(BPcum=BP+tot)
    axisdf = don %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2, minBPcum = min(BPcum))
    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
        "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
    don_sub = don # %>% group_by(icd) %>% mutate(ind = 1:length(LDblock)) #%>% filter(ind < 100)
    don_sub$ICDanno = factor(don_sub$ICDanno, levels = rev(names(sort(table(don_sub$ICDanno)))))
    don_sub1 = reshape2::melt(table(unique(don_sub[,c("ICDanno","LDblock")])[,1])) %>% mutate(percent = value/length(unique(don_sub$LDblock))) %>% arrange(desc(value))
    don_sub1$percent = paste0(signif(don_sub1$percent,3)*100, "%")
    don_sub1$Var1 = factor(don_sub1$Var1, levels = don_sub1$Var1)
    don_sub1$percent1 = don_sub1$percent
    p1 = ggplot(don_sub, aes(x=BPcum, y=-log10(P), color=ICDanno, size=OR1)) +
        geom_point(alpha=0.8, position=position_jitter(h = 0, w = 0.5)) +
        geom_hline(yintercept=7.3, linetype = "dashed") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 5))) +
        scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_size_manual(breaks = c("(0,1)", "[1,5)", "[5,10)", "[10,15)", ">=15"), values = c(1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(10, 10), scale='free', space = 0) +
        scale_y_cut(10) +
        theme_bw() +
        theme(
            axis.text = element_text(size=10),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            legend.position= 'bottom',
            legend.direction = "vertical",
            legend.box = "horizontal",
            panel.grid.major.x = element_blank())
    p2 = ggplot(don_sub1 %>% mutate(Var1 = factor(Var1, levels = rev(Var1))), aes(value, Var1, fill=Var1, label = percent)) + geom_col() +
        scale_fill_manual(values=rev(colors37[1:nrow(don_sub1)])) +
        # coord_fixed(ratio = 1/1000) +
        scale_x_sqrt(expand=c(0,0)) +
        geom_text(aes(x=max(don_sub1$value) * 0.75)) +
        theme_classic() +
        theme(axis.title=element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            legend.position= 'none') #axis.text.x = element_text(angle = 60, vjust=1,hjust=1),
    return(list(p1, p2))
    # pdf(paste0(title, ".pdf"), width = 15, heigh = 6)
    # plot(p1)
    # dev.off()
    # pdf(paste0(title, ".count.pdf"), width = 8, heigh = 3)
    # plot(p2)
    # dev.off()
}

pdf("Figure3A.pdf", width=15, height=6)
plot_all_icd_snp(res_sig_raw[res_sig_raw$p.value<1e-6, ], title = "res_sig_raw_manhantan_ld_1e6", ldblock = T)
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

## Figure 3D-F

```shell
export PATH=/data/apps/locuszoom/bin:$PATH
source activate locuszoom

zcat N97.SAIGE.txt.gz|grep "^2:" | cut -f 1,6 | sed 's/:/\t/g' | awk '{print $1":"$2"\t"$5}' > N97.chr2.txt
sed -i 1i"MarkerName\tP-value" N97.chr2.txt
locuszoom --metal N97.chr2.txt \
    --ld-vcf /data/NIPT/1000G/Eagle_Minimac4_GRCh38_positions_Reference_panels/CHB_CHS/ALL.chr2_GRCh38.genotypes.20170504.norm.vcf.gz \
    --build hg38 --chr 2 --start 110000000 --end 113000000 signifLine="7.3" --rundir N97 \
    --refsnp 2:111056516 --add-refsnps 2:110160461

zcat D56.SAIGE.txt.gz|grep ^16: | cut -f 1,6 | sed 's/:/\t/g' | awk '{print $1":"$2"\t"$5}' > D56.chr16.txt
sed -i 1i"MarkerName\tP-value" D56.chr16.txt
locuszoom --metal D56.chr16.txt \
    --ld-vcf /data/NIPT/1000G/Eagle_Minimac4_GRCh38_positions_Reference_panels/CHB_CHS/ALL.chr16_GRCh38.genotypes.20170504.norm.vcf.gz \
    --build hg38 --chr 16 --start 100000 --end 500000 signifLine="7.3" --rundir D56

zcat N70.SAIGE.txt.gz|grep ^2: | cut -f 1,6 | sed 's/:/\t/g' | awk '{print $1":"$2"\t"$5}' > N70.chr2.txt
sed -i 1i"MarkerName\tP-value" N70.chr2.txt
locuszoom --metal N70.chr2.txt \
    --ld-vcf /data/NIPT/1000G/Eagle_Minimac4_GRCh38_positions_Reference_panels/CHB_CHS/ALL.chr2_GRCh38.genotypes.20170504.norm.vcf.gz \
    --build hg38 --chr 2 --start 107000000 --end 110000000 signifLine="7.3" --rundir N70 \
    --refsnp 2:108897145 --add-refsnps 2:109220369

zcat E28.SAIGE.txt.gz|grep ^10: | cut -f 1,6 | sed 's/:/\t/g' | awk '{print $1":"$2"\t"$5}' > E28.chr10.txt
sed -i 1i"MarkerName\tP-value" E28.chr10.txt
locuszoom --metal E28.chr10.txt \
    --ld-vcf /data/NIPT/1000G/Eagle_Minimac4_GRCh38_positions_Reference_panels/CHB_CHS/ALL.chr10_GRCh38.genotypes.20170504.norm.vcf.gz \
    --build hg38 --chr 10 --start 48000000 --end 49500000 signifLine="7.3" --rundir E28

zcat O36.SAIGE.txt.gz|grep ^7: | cut -f 1,6 | sed 's/:/\t/g' | awk '{print $1":"$2"\t"$5}' > O36.chr7.txt
sed -i 1i"MarkerName\tP-value" O36.chr7.txt
locuszoom --metal O36.chr7.txt \
    --ld-vcf /data/NIPT/1000G/Eagle_Minimac4_GRCh38_positions_Reference_panels/CHB_CHS/ALL.chr7_GRCh38.genotypes.20170504.norm.vcf.gz \
    --build hg38 --chr 7 --start 95000000 --end 105000000 signifLine="7.3" --rundir O36
```
