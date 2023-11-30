# Figure4 and S6

```R
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggbreak)
library(ggrepel)

plot_hist_density = function(x, breaks = 10, xlab = ""){
    res = density(x)
    hist(x, breaks, prob=T, ylim = c(0, max(res$y)), main = "", xlab=xlab)
    lines(res,lwd=2)
    # text(res$x[which.max(res$y)], max(res$y)*1.02, labels = signif(res$x[which.max(res$y)],3))
}

library(dplyr)
load("bed_cnv.RData")
bed_cnv$region = paste(bed_cnv$chr, bed_cnv$start, bed_cnv$end,sep=":")

bed = unique(bed_cnv[,c(1:3,8)])

# board
cytoband_bed = read.table("cytoband_hg38.txt", header=F,sep="\t",stringsAsFactors=F)
colnames(cytoband_bed) = c("chrom", "start", "end", "cytoband", "V5")
board_bed = cytoband_bed
board_bed$board = ""
for(i in c(paste0(c(1:22, "X","Y"), "p"), paste0(c(1:22, "X","Y"), "q"))){
    board_bed$board[grep(paste0("^", i), board_bed$cytoband)] = i
}
board_bed = board_bed %>% group_by(board) %>% summarise(chrom=unique(chrom), start = min(start), end = max(end))
board_bed = board_bed[, c(2:4,1)]
bed = bed_cnv[,c(1:3,6:8)] %>% group_by(chr, start, end, region, type) %>% summarise(n = length(unique(ID)))

bed_board = unique(bedtoolsr::bt.intersect(bed, board_bed, wb=T))
colnames(bed_board) = c("chr", 'start', 'end', 'region', "type", "n", "chrom", "board_start","board_end","board", "CHR")
bed_board$fraction = c(bed_board$end-bed_board$start+1)/c(bed_board$board_end-bed_board$board_start)
bed_board = merge(bed_cnv, bed_board, by=c("chr", 'start', 'end', 'region', "type"))

load("clin.RData")

pdf("FigureS7A.pdf", width=4,height=4)
plot_hist_density(bed_board$fraction, breaks=100, xlab="fraction")
dev.off()


board_bed$CHR = as.numeric(gsub("Y", "24", gsub("X", 23, board_bed$chrom)))
board_bed = board_bed %>% arrange(CHR, start, end)

res = read.table("pheWAS_status_PC_cnv_board.txt", header=T,sep="\t",stringsAsFactors=F)
res = merge(res, board_bed[,1:4], by.x = "cyto", by.y = "board")
res$OR = exp(res$Est)
res$cnv = as.numeric(res$cnv)
res$ICDanno = ""
res$ICDanno[grep("[AB]", res$icd)] = "Certain infectious and parasitic diseases"
res$ICDanno[which(res$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24", "D25", "D26"))] = "Neoplasms"
res$ICDanno[which(res$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
res$ICDanno[grep("E", res$icd)] = "Endocrine, nutritional and metabolic diseases"
res$ICDanno[grep("F", res$icd)] = "Mental and behavioural disorders"
res$ICDanno[which(res$icd %in% c("H00", "H01", "H02", "H04", "H10","H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
res$ICDanno[which(res$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92","H93"))] = "Ear and mastoid process"
res$ICDanno[grep("I", res$icd)] = "Circulatory system"
res$ICDanno[grep("J", res$icd)] = "Respiratory system"
res$ICDanno[grep("K", res$icd)] = "Digestive system"
res$ICDanno[grep("L", res$icd)] = "Skin and subcutaneous tissue"
res$ICDanno[grep("M", res$icd)] = "Musculoskeletal system and connective tissue"
res$ICDanno[grep("N", res$icd)] = "Genitourinary system"
res$ICDanno[grep("O", res$icd)] = "Pregnancy, childbirth and the puerperium"
res$ICDanno[grep("P", res$icd)] = "Certain conditions originating in the perinatal period"
res$ICDanno[grep("Q", res$icd)] = "Congenital malformations, deformations and chromosomal abnormalities"

res$ICDanno = factor(res$ICDanno,
                            levels = c("Genitourinary system", "Pregnancy, childbirth and the puerperium", "Endocrine, nutritional and metabolic diseases", "Digestive system", "Respiratory system", "Musculoskeletal system and connective tissue", "Skin and subcutaneous tissue", 
                                "Blood and blood-forming organs and certain disorders involving the immune mechanism",
                                "Certain infectious and parasitic diseases", "Ear and mastoid process", "Neoplasms", "Circulatory system", 
                                "Congenital malformations, deformations and chromosomal abnormalities", "Mental and behavioural disorders", "Eye and adnexa"))


p_cnv_board_filtered = bed_board %>% filter(fraction>0.8, board %in% names(which(colSums(mat_board!=0)>10))) %>%
    mutate(board = factor(board, levels = unlist(lapply(1:22, function(x){paste0(x, c("p","q"))}))),
        type = factor(type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain"))) %>%
    arrange(ratio) %>% mutate(index = 1:length(ratio)) %>%
    ggplot(aes(x=index, y=ratio, color = type)) +
        geom_point() + labs(x="", y="log2(Ratio)", color="") +
        facet_grid(~board, scale="free") +
        scale_color_d3(guide = guide_legend(override.aes = list(size = 3))) +
        theme_bw(base_size=15) +
        theme(axis.text.x=element_blank(),
            axis.text = element_text(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(),
            axis.title.x=element_blank(),
            legend.title = element_text(),
            legend.text = element_text(),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = c(0.2, 0.2),
            legend.direction = "horizontal",
            legend.background = element_blank(),
            legend.box = "horizontal",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())

pdf("FigureS7B.pdf", width=12, height=2.285714)
p_cnv_board_filtered
dev.off()

plot_all_icd_cnv_filterd = function(res_sig, title, maxP = 5e-8, legend_pos = "bottom"){
    res_sig$OR1 = res_sig$OR
    res_sig$OR1[res_sig$OR >= 30] = ">=30"
    res_sig$OR1[res_sig$OR >= 15 & res_sig$OR < 30] = "[15,30)"
    res_sig$OR1[res_sig$OR >= 10 & res_sig$OR < 15] = "[10,15)"
    res_sig$OR1[res_sig$OR >= 5 & res_sig$OR < 10] = "[5,10)"
    res_sig$OR1[res_sig$OR > 1 & res_sig$OR < 5] = "(1,5)"
    res_sig$OR1[res_sig$OR > 0 & res_sig$OR < 1] = "(0,1)"
    res_sig$cyto = factor(res_sig$cyto, levels = board_bed$board)
    res_sig$pos = rowMeans(res_sig[,c("start","end")])
    res_sig$CHR = factor(as.numeric(gsub("Y", "24", gsub("X", 23, res_sig$chrom))), levels = 1:24)
    res_sig$logP = -log10(res_sig$P) * res_sig$cnv
    res_sig$logP[res_sig$logP > 10] = 10
    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
        "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
    res_sig$icd[res_sig$P >= maxP] = NA
    res_sig$lab = res_sig$icd
    res_sig$lab[is.na(res_sig$lab)] = ""
    p = res_sig %>% group_by(cyto) %>% arrange(logP) %>% mutate(index=1:length(logP)) %>%
        ggplot(aes(x=index, y=logP, color=ICDanno, size = OR1)) +
        geom_point(alpha=0.7) +
        facet_grid(~cyto, scale="free")+
        geom_text_repel(aes(label = lab), size = 5) +
        labs(color = "", size = "Odds ratio", y = "Deletion-Amplication\n\n-log10(Pvalue)") +
        geom_hline(yintercept = 0) +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 3)), na.value = "grey") +
        scale_y_continuous(breaks= c(-5,0,5,10), labels = c(5, 0, 5, ">10")) +
        scale_size_manual(breaks = c("(0,1)", "(1,5)", "(5,10)", "[10,15)", "[15,30)", ">=30"), values = c(0.5, 1, 1.5, 2, 2.5, 3)) +
        theme_bw(base_size=15) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(),
            axis.title.x=element_blank(),
            legend.title = element_text(),
            legend.text = element_text(),
            legend.key.size = unit(1,"line"),
            legend.position = legend_pos,
            legend.direction = "horizontal",
            legend.box = "vertical",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())
    return(p)
}

p_cnv_icd_filtered = plot_all_icd_cnv_filterd(res[res$P<0.05 & res$cyto %in% names(which(colSums(mat_board!=0)>10)), ], 
    title = "CNA-ICD10", maxP = 1e-4, legend_pos = "bottom")

pdf("Figure4.pdf", width=12, height=5.714286)
p_cnv_icd_filtered
dev.off()
```
