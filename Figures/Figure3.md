# Figure 3

## Figure 3A

```R
plot_all_icd_snp = function(res_sig, title, ldblock = F, h_SNP){
    if(ldblock){
        res_sig1 = res_sig %>% group_by(ICDanno, LDblock) %>%
            summarise(CHR = CHR[which.min(p.value)], BP = BP[which.min(p.value)], OR = OR[which.min(p.value)], BETA = BETA[which.min(p.value)], P = min(p.value), cytoband = cytoband[which.min(p.value)]) %>% as.data.frame(stringsAsFactors = F)
    }else{
        res_sig = res_sig %>% mutate(P = p.value, OR_mean = OR)
    }
    res_sig1$OR1 = res_sig1$OR
    res_sig1$OR1[res_sig1$OR >= 15] = ">=15"
    res_sig1$OR1[res_sig1$OR >= 10 & res_sig1$OR < 15] = "[10,15)"
    res_sig1$OR1[res_sig1$OR >= 5 & res_sig1$OR < 10] = "[5,10)"
    res_sig1$OR1[res_sig1$OR >= 1 & res_sig1$OR < 5] = "[1,5)"
    res_sig1$OR1[res_sig1$OR > 0 & res_sig1$OR < 1] = "(0,1)"
    don <- res_sig1 %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len=max(as.numeric(BP))) %>% 
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        # Add this info to the initial dataset
        left_join(res_sig1, ., by=c("CHR"="CHR")) %>%
        # Add a cumulative position of each SNP
        arrange(CHR, BP) %>%
        mutate(BPcum=BP+tot)
    axisdf = don %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2, minBPcum = min(BPcum))
    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c",#"#a06fda","#df462a",
    "#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
    don_sub = don # %>% group_by(icd) %>% mutate(ind = 1:length(LDblock)) #%>% filter(ind < 100)
    don_sub$ICDanno = factor(don_sub$ICDanno,
                            levels = c("Genitourinary system", "Pregnancy, childbirth and the puerperium", 
                                "Endocrine, nutritional and metabolic diseases", 
                                "Digestive system", "Respiratory system", 
                                "Musculoskeletal system and connective tissue", 
                                "Skin and subcutaneous tissue", 
                                "Blood and blood-forming organs and certain disorders involving the immune mechanism",
                                "Certain infectious and parasitic diseases", "Ear and mastoid process", "Neoplasms", "Circulatory system", 
                                "Eye and adnexa"))
    don_sub1 = reshape2::melt(table(unique(don_sub[,c("ICDanno","LDblock")])[,1])) %>% mutate(percent = value/length(unique(don_sub$LDblock))) %>% arrange(desc(value))
    don_sub1$percent = paste0(signif(don_sub1$percent,3)*100, "%")
    don_sub1$Var1 = factor(don_sub1$Var1, levels = don_sub1$Var1)
    don_sub1$percent1 = don_sub1$percent
    h_SNP = cbind(do.call(rbind, lapply(1:nrow(h_SNP), function(x){
        m = res_sig[which(res_sig$MarkerID==h_SNP$V1[x] & res_sig$icd==h_SNP$V2[x]), ]
        mm = don_sub[which(don_sub$LDblock %in% m$LDblock & don_sub$BP %in% m$BP & don_sub$ICDanno %in% m$ICDanno),]
        return(mm)
    })), h_SNP)
    h_SNP$lab = paste(h_SNP$V1, h_SNP$V2)
    don_sub = left_join(don_sub, h_SNP[,c("LDblock","ICDanno","lab")], by = c("ICDanno", "LDblock"))
    p0 = ggplot(don_sub, aes(x=BPcum, y=-log10(P), color=ICDanno, size=OR1, label=lab)) +
        geom_point(alpha=0.8, position=position_jitter(h = 0, w = 0.5)) +
        geom_hline(yintercept=7.3, linetype = "dashed") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 5))) +
        scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_size_manual(breaks = c("(0,1)", "[1,5)", "[5,10)", "[10,15)", ">=15"), values = c(1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(10, 10), scale='free', space = 0) +
        geom_text_repel(force=20, color="grey20", size=3, point.padding = 1, hjust = 1, vjust = 1,
            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
            segment.color="grey20", segment.size=0.2, segment.alpha=1, nudge_y=1, nudge_x = 1) +
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
    p1 = ggplot(don_sub, aes(x=BPcum, y=-log10(P), color=ICDanno, size=OR1, label=lab)) +
        geom_point(alpha=0.8, position=position_jitter(h = 0, w = 0.5)) +
        geom_hline(yintercept=7.3, linetype = "dashed") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 5))) +
        scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_size_manual(breaks = c("(0,1)", "[1,5)", "[5,10)", "[10,15)", ">=15"), values = c(1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(10, 10), scale='free', space = 0) +
        geom_text_repel(force=20, color="grey20", size=3, point.padding = 1, hjust = 1, vjust = 1,
            arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
            segment.color="grey20", segment.size=0.2, segment.alpha=1, nudge_y=1, nudge_x = 1) +
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
    return(list(p0, p1, p2))
}

load("pheWAS_sig_SAIGE_1e6_power0.8_filtered.RData")
# res_sig_1e6_filter, 5e-8 <= Pvalue < 1e-6
res_sig_1e6 = rbind(res_sig_1e6_filter, res_sig_filter) # get all SNPs with Pvalue < 1e-6
h_SNP = c("1:147155826:C:G", "N97", "2:108999920:T:C", "N97", "6:61238993:A:G", "N81", 
    "7:67431233:CA:C", "E61",
    "10:48576703:C:G", "O02", "10:48551720:G:A", "E28", "12:80511378:C:T", "K13", 
    "16:268762:G:A","D56", "16:60232728:C:T", "O06")
h_SNP = as.data.frame(matrix(h_SNP, ncol=2, byrow=T), stringsAsFactors=F)

pdf("Figure3A.pdf", width=15, height=6)
plot_all_icd_snp(res_sig_1e6, title = "res_sig_raw_manhantan_ld_1e6", ldblock = T, h_SNP)
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
