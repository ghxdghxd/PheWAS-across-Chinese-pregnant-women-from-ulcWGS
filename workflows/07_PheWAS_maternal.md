# 7. PheWAS for mother

## 7.1 PheWAS for maternal variants

### 7.1.1 plink files for SAIGE

```R
library(data.table)
library(stringr)
library(dplyr)
load("all.1kg.geno.RData")
geno_1000G = as.data.frame(geno_1000G, stringsAsFactors=F)

load("all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.com0.1.callRate0.99.conf95.GT1.RData")
int = intersect(colnames(geno), colnames(geno_1000G))
geno = geno[, ..int]
geno = geno[!duplicated(sample_ID),]
sample_ID = sample_ID[!duplicated(sample_ID)]
load("/data/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge_Model/icd_status_PC_all.RData")
geno = geno[which(sample_ID %in% rownames(clin_status)), ]
sample_ID = sample_ID[which(sample_ID %in% rownames(clin_status))]
geno = t(geno)
colnames(geno) = sample_ID
geno[is.na(geno)] = "./."
geno[geno==0] = "0/0"
geno[geno==1] = "0/1"
geno[geno==2] = "1/1"
vcf = cbind(str_split(rownames(geno), pattern = ":", simplify = T), geno)
colnames(vcf)[1:4] = c("#CHROM", "POS", "REF", "ALT")
vcf = as.data.frame(vcf, stringsAsFactors = F)
vcf = cbind(vcf[,1:2], ID = rownames(vcf), vcf[,3:4], QUAL = ".", FILTER="PASS", INFO = ".", FORMAT = "GT", vcf[,-c(1:4)])
rownames(vcf) = NULL
vcf$`#CHROM` = as.numeric(vcf$`#CHROM`)
vcf$POS = as.numeric(vcf$POS)
vcf = vcf %>% arrange(`#CHROM`, POS)

samples = read.table("samples_25639.txt",header=F,sep=" ",stringsAsFactors=F)
samples = samples$V1
names(samples) = gsub("[AB]$", "", samples)
colnames(vcf)[-c(1:9)] = samples[colnames(vcf)[-c(1:9)]]

write.table('##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=file:///data/hg38/hg38.fa.gz
##contig=<ID=1,length=248956422>
##contig=<ID=2,length=242193529>
##contig=<ID=3,length=198295559>
##contig=<ID=4,length=190214555>
##contig=<ID=5,length=181538259>
##contig=<ID=6,length=170805979>
##contig=<ID=7,length=159345973>
##contig=<ID=8,length=145138636>
##contig=<ID=9,length=138394717>
##contig=<ID=10,length=133797422>
##contig=<ID=11,length=135086622>
##contig=<ID=12,length=133275309>
##contig=<ID=13,length=114364328>
##contig=<ID=14,length=107043718>
##contig=<ID=15,length=101991189>
##contig=<ID=16,length=90338345>
##contig=<ID=17,length=83257441>
##contig=<ID=18,length=80373285>
##contig=<ID=19,length=58617616>
##contig=<ID=20,length=64444167>
##contig=<ID=21,length=46709983>
##contig=<ID=22,length=50818468>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file="for_SAIGE.vcf", r=F,c=F,sep="\t",quote=F)

write.table(vcf, file="for_SAIGE.vcf", r=F,c=T,sep="\t", quote=F, append =T)
```

```sh
plink2 --vcf for_SAIGE.vcf.gz --make-bed --out for_SAIGE

awk -F "\t" '{$5=2;OFS="\t";print $0}' for_SAIGE.fam > for_SAIGE.fam.1
mv for_SAIGE.fam.1 for_SAIGE.fam
```

### 7.1.2 run SAIGE

```sh
step1_fitNULLGLMM.R     \
    --plinkFile=for_SAIGE \
    --phenoFile=clin_status.gender.txt \
    --phenoCol=$icd \
    --covarColList=age \
    --sampleIDColinphenoFile=ID \
    --traitType=binary \
    --outputPrefix=./output/$icd \
    --nThreads=30 \
    --LOCO=TRUE \
    --IsOverwriteVarianceRatioFile=TRUE

step2_SPAtests.R \
    --vcfFile=$j \
    --vcfFileIndex=$j.tbi \
    --vcfField=GT \
    --chrom=$chr \
    --minMAF=0.01 \
    --minMAC=0.5 \
    --sampleFile=samples_25639.txt \
    --GMMATmodelFile=./output1/$i.rda \
    --varianceRatioFile=./output1/$i.varianceRatio.txt \
    --SAIGEOutputFile=./output2/$i.SAIGE.genotype.$chr.txt \
    --LOCO=TRUE

```

### 7.1.3 combine SAIGE results

```R
files = list.files("output2", pattern="*.txt", full.names=T)
library(data.table)
imputed_info = fread("snp_imputed_R0.5_MAF0.01.txt",header=T,sep="\t",stringsAsFactor=F)

load("icd_status_PC_all_BMI.RData")
clin_status = clin_status[,names(which(colSums(clin_status[,-c(1:7)], na.rm = T)>30))]

library(parallel)

for(icd in colnames(clin_status)){
    res_SAIGE = mclapply(grep(icd, files,v=T), function(x){
        print(x)
        m = fread(x, header=T,sep="\t", stringsAsFactors = F)
        m = m[m$MarkerID %in% imputed_info$SNP, c("MarkerID", "BETA", "SE", "Tstat", "var", "p.value", "Is.SPA", "AF_case", "AF_ctrl")]
        return(m)
    }, mc.cores = 22)
    res_SAIGE = do.call(rbind, res_SAIGE)
    fwrite(res_SAIGE, file = paste0(icd, ".SAIGE.txt.gz"), sep="\t", row.names = F, col.names = T, quote = F)
}
```

### 7.1.4 calculating statistical power

```R
load("icd_status_PC_all_BMI.RData")
clin_status = clin_status[,names(which(colSums(clin_status[,-c(1:7)], na.rm = T)>30))]

# observed MAF
load("all_R2_0.5_MAF_in_motherICD_childICD.RData")

library(genpwr)

case_num = colSums(clin_status, na.rm=T)
case_control_num = colSums(!is.na(clin_status), na.rm=T)
case_rate = case_num/case_control_num

library(data.table)
library(stringr)
library(parallel)

for(x in colnames(clin_status)[3:267]){
    res = fread(paste0(x, ".SAIGE.txt.gz"), header=T, sep="\t", stringsAsFactors=F)
    res = res[res$p.value < 1e-3, ]
    res$OR = exp(res$BETA)
    gc()
    pow = mclapply(1:nrow(res), function(i){
        print(i)
        res = try(genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                N=case_control_num[x],
                Case.Rate=case_rate[x],
                k=NULL,
                MAF=res_MAF_M[res$MarkerID[i],"MAF"],
                OR=res$OR[i],
                Alpha=0.05,
                True.Model=c("Additive"),
                Test.Model=c("Additive")))
        if(class(res) == "try-error"){
            return(NA)
        }
        return(res$Power_at_Alpha_0.05)
    }, mc.cores = 30)
    res$power = unlist(pow)
    fwrite(res, file=paste0(x, ".SAIGE.power.txt.gz"), row.names = F, col.names = T, sep="\t",quote = F)
}

snp_power = lapply(list.files(".", pattern="*.txt.gz", full.names=T), function(x){
    m = fread(x, header=T,sep='\t',stringsAsFactors = F)
    m$icd = gsub(".SAIGE.power.txt.gz", "", basename(x))
    return(m)
})
snp_power = do.call(rbind, snp_power)
save(snp_power, file="phewas_power.RData")
```

## 7.2 PheWAS for maternal CNVs

### 7.2.1 CNVs of chromosome boards

## cnv_cytoband, cnv_board and cnv_gene

```R
plot_hist_density = function(x, breaks = 10, xlab = ""){
    res = density(x)
    hist(x, breaks, prob=T, ylim = c(0, max(res$y)), main = "", xlab=xlab)
    lines(res,lwd=2)
    # text(res$x[which.max(res$y)], max(res$y)*1.02, labels = signif(res$x[which.max(res$y)],3))
}

# data6
library(dplyr)
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/bed_cnv.RData")
bed_cnv$region = paste(bed_cnv$chr, bed_cnv$start, bed_cnv$end,sep=":")

# filter ratio
# bed_cnv = bed_cnv[bed_cnv$ratio > 0.3 | bed_cnv$ratio < -0.5, ]

bed = unique(bed_cnv[,c(1:3,8)])

# cytoband
cytoband_bed = read.table("/data/NIPT/cytoband_hg38.txt", header=F,sep="\t",stringsAsFactors=F)
colnames(cytoband_bed) = c("chrom", "start", "end", "cytoband", "V5")

bed_cytoband = unique(bedtoolsr::bt.intersect(bed, cytoband_bed, wb=T))
bed_cytoband = bed_cytoband[,c(4,8)]

bed_cytoband = merge(bed_cnv, bed_cytoband, by.x="region", by.y="V4")
bed_cytoband_1 = bed_cytoband %>% group_by(ID,V8) %>% summarise(type=paste(sort(unique(type)), collapse=","), ratio = mean(ratio))

mat_cytoband = reshape2::acast(bed_cytoband_1[,1:3], ID ~ V8)

mat_cytoband[mat_cytoband=="gain"] = 1
mat_cytoband[mat_cytoband=="loss"] = -1
mat_cytoband[is.na(mat_cytoband)] = 0
mat_cytoband[mat_cytoband=="gain,loss"] = NA
save(mat_cytoband, bed_cytoband_1, file="/data/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_cytoband.RData")

# filter ratio
#save(mat_cytoband, bed_cytoband_1, file="/data/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_cytoband_filtered.RData")

# board
cytoband_bed = read.table("/share/data3/NIPT/cytoband_hg38.txt", header=F,sep="\t",stringsAsFactors=F)
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
colnames(bed_board) = c("chr", 'start', 'end', 'region', "type", "n", "chrom", "board_start","board_end","board")
bed_board$fraction = c(bed_board$end-bed_board$start+1)/c(bed_board$board_end-bed_board$board_start)
bed_board = merge(bed_cnv, bed_board, by=c("chr", 'start', 'end', 'region', "type"))

load("../pheWAS_sig_merge/for_SAIGE/fq_ICD_clin.fix_NA.fixColumn.RData")

# cbind(bed_board[match(intersect(ICD_clin_mother$ID, bed_board$ID), bed_board$ID), ],
#     ICD_clin_mother[match(intersect(ICD_clin_mother$ID, bed_board$ID), ICD_clin_mother$ID), grep("唐筛", colnames(ICD_clin_mother))]) %>%
#     filter(fraction>0.8, chrom %in% c(13, 18, 21)) %>% select("board", grep("唐筛", colnames(ICD_clin_mother),v=T)) %>% rowwise %>%
#     mutate(anno= paste(`唐筛结果.筛查模式`, `唐筛结果.早期母血清筛查风险率`, `唐筛结果.早期18三体综合征.1.`, `唐筛结果.早期21三体综合征.1.`, `唐筛结果.中期母血清筛查风险率`, `唐筛结果.中期18三体综合征.1.`, `唐筛结果.中期21三体综合征.1.`, `唐筛结果`, collapse = "_")) %>%
#     select(anno) %>% 

pdf("board_fraction.pdf", width=4,height=4)
plot_hist_density(bed_board$fraction, breaks=100, xlab="fraction")
dev.off()

bed_board_1 = bed_board %>% filter(fraction > 0.8, !grepl("X", board)) %>%
    group_by(ID, board) %>% summarise(type=paste(sort(unique(type)), collapse=","), ratio = mean(ratio), zscore = mean(zscore))

#load("/share/data3/NIPT/clin/fq_clin.fixColumn.RData")

mat_board = reshape2::acast(bed_board_1[,1:3], ID ~ board)

mat_board[mat_board=="gain"] = 1
mat_board[mat_board=="loss"] = -1
mat_board[is.na(mat_board)] = 0
save(mat_board, bed_board_1, file="/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_board.RData")

# gene
gene_bed = read.table("/share/data3/NIPT/bam_vcf/hg38.RefSeqGene.bed", header=F,sep="\t",stringsAsFactors=F)
colnames(gene_bed) = c("chrom", "start", "end", "gene")
gene_bed$chrom = gsub("chr", "", gene_bed$chrom)
gene_bed = gene_bed[-grep("_", gene_bed$chrom), ]
gene_bed = gene_bed %>% group_by(gene) %>% summarise(chrom=unique(chrom), start=min(start), end = max(end))
gene_bed = gene_bed[, c(2:4, 1)]
bed_gene = bedtoolsr::bt.intersect(bed, gene_bed, wb=T)
bed_gene = unique(bed_gene[,c(4,8)])

bed_gene = merge(bed_cnv, bed_gene, by.x="region", by.y="V4")
bed_gene_1 = bed_gene %>% group_by(ID,V8) %>% summarise(type=paste(sort(unique(type)), collapse=","))

mat_gene = reshape2::acast(bed_gene_1, ID ~ V8)

mat_gene[mat_gene=="gain"] = 1
mat_gene[mat_gene=="loss"] = -1
mat_gene[is.na(mat_gene)] = 0
mat_gene[mat_gene=="gain,loss"] = NA
save(mat_gene, file="/data/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_gene.RData")

# filter ratio
#save(mat_gene, file="/data/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_gene_filtered.RData")
```

## cnv_ICD10

### cytoband

```R
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_cytoband.RData")
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge_Model/icd_status_PC_all_BMI.RData")
clin_status$E66 = NULL
clin_status$E66 = ifelse(clin_status$BMI > 24, 1, 0)
clin_status$BMI = NULL
clin_status = clin_status[, grep("\\.1", colnames(clin_status),invert=T)]
clin_status = clin_status[, grep("R", colnames(clin_status),invert=T)]
int = intersect(rownames(clin_status), rownames(mat_cytoband))
mat_cytoband_sub = as.data.frame(mat_cytoband[int, ], stringsAsFactors = F)
clin_status_sub = clin_status[int, ]
clin_status_sub = clin_status_sub[, c(colnames(clin_status_sub)[1:7], names(which(colSums(clin_status_sub[,-c(1:7)], na.rm = T)>30)))]

# cnv_count = sapply(names(mat_cytoband_sub), function(x){
#     m = table(mat_cytoband_sub[, x])
#     m = c(m["0"], m["1"], m["-1"])
#     names(m) = c("Neu", "Amp", "Del")
#     return(m)
# })

library(parallel)
res = mclapply(colnames(mat_cytoband_sub), function(x){
    res = lapply(colnames(clin_status_sub)[-c(1:7)], function(y){
        print(paste(x, y))
        dat = na.omit(cbind(clin_status_sub[,c(y, "age","PC1","PC2","PC3","PC4","PC5")], cnv = mat_cytoband_sub[,x]))
        if(length(unique(dat$cnv))>1){
            colnames(dat)[1]='icd'
            res_glm = coef(summary(glm(icd ~ cnv + age + PC1 + PC2 + PC3 + PC4 + PC5, data = dat, family="binomial")))[, c(1,4)]
            res_glm = as.data.frame(res_glm[c(1, grep("cnv", rownames(res_glm))),], stringsAsFactors=F)
            colnames(res_glm) = c("Est", "P")
            res_glm$icd = y
            res_glm$cnv = rownames(res_glm)
            return(res_glm[-1, ])
        }else{
            return(NA)
        }
    })
    res = do.call(rbind, res[which(sapply(res, length)>1)])
    if(!is.null(res)){
        res$cyto = x
        res$cnv = gsub("cnv", "", res$cnv)
    }else{
        res = NA
    }
    return(res)
}, mc.cores = 35)
res = do.call(rbind, res)
res$Est = as.numeric(res$Est)
res$P = as.numeric(res$P)
rownames(res) = NULL
res = na.omit(res)
write.table(res, file="pheWAS_status_PC_cnv_cytoband.txt", r=F,c=T,sep="\t",quote=F)
# filter ratio
#write.table(res, file="pheWAS_status_PC_cnv_filtered_cytoband.txt", r=F,c=T,sep="\t",quote=F)
###############################

library(dplyr)
library(ggplot2)
library(ggsci)
library(ggbreak)
library(ggrepel)

cytoband_bed = read.table("/share/data3/NIPT/cytoband_hg38.txt", header=F,sep="\t",stringsAsFactors=F)
colnames(cytoband_bed) = c("chrom", "start", "end", "cytoband", "V5")
cytoband_bed$CHR = as.numeric(gsub("Y", "24", gsub("X", 23, cytoband_bed$chrom)))
cytoband_bed = cytoband_bed %>% arrange(CHR, start, end)

res = read.table("pheWAS_status_PC_cnv_cytoband.txt", header=T,sep="\t",stringsAsFactors=F)
res = merge(res, cytoband_bed[,1:4], by.x = "cyto", by.y = c("cytoband"))
res$OR = exp(res$Est)

res_filterd = read.table("pheWAS_status_PC_cnv_filtered_cytoband.txt", header=T,sep="\t",stringsAsFactors=F)
res_filterd = merge(res_filterd, cytoband_bed[,1:4], by.x = "cyto", by.y = c("cytoband"))
res_filterd$OR = exp(res_filterd$Est)

res_filterd$ICDanno = ""
res_filterd$ICDanno[grep("[AB]", res_filterd$icd)] = "Certain infectious and parasitic diseases"
res_filterd$ICDanno[which(res_filterd$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24", "D25", "D26"))] = "Neoplasms"
res_filterd$ICDanno[which(res_filterd$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
res_filterd$ICDanno[grep("E", res_filterd$icd)] = "Endocrine, nutritional and metabolic diseases"
res_filterd$ICDanno[grep("F", res_filterd$icd)] = "Mental and behavioural disorders"
res_filterd$ICDanno[which(res_filterd$icd %in% c("H00", "H01", "H02", "H04", "H10","H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
res_filterd$ICDanno[which(res_filterd$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92","H93"))] = "Ear and mastoid process"
res_filterd$ICDanno[grep("I", res_filterd$icd)] = "Circulatory system"
res_filterd$ICDanno[grep("J", res_filterd$icd)] = "Respiratory system"
res_filterd$ICDanno[grep("K", res_filterd$icd)] = "Digestive system"
res_filterd$ICDanno[grep("L", res_filterd$icd)] = "Skin and subcutaneous tissue"
res_filterd$ICDanno[grep("M", res_filterd$icd)] = "Musculoskeletal system and connective tissue"
res_filterd$ICDanno[grep("N", res_filterd$icd)] = "Genitourinary system"
res_filterd$ICDanno[grep("O", res_filterd$icd)] = "Pregnancy, childbirth and the puerperium"
res_filterd$ICDanno[grep("P", res_filterd$icd)] = "Certain conditions originating in the perinatal period"
res_filterd$ICDanno[grep("Q", res_filterd$icd)] = "Congenital malformations, deformations and chromosomal abnormalities"

res_filterd$ICDanno = factor(res_filterd$ICDanno,
                            levels = c("Pregnancy, childbirth and the puerperium", "Genitourinary system", "Endocrine, nutritional and metabolic diseases", "Digestive system", "Respiratory system", "Skin and subcutaneous tissue", "Musculoskeletal system and connective tissue",
                                "Certain infectious and parasitic diseases", "Eye and adnexa", "Ear and mastoid process", "Neoplasms", "Circulatory system", "Mental and behavioural disorders", "Blood and blood-forming organs and certain disorders involving the immune mechanism",
                                "Congenital malformations, deformations and chromosomal abnormalities"))

load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/bed_cnv.RData")
bed_cnv$region = paste(bed_cnv$chr, bed_cnv$start, bed_cnv$end,sep=":")
cnv = bed_cnv
cnv = cnv[cnv$chr!="X", ]
cnv$chr = factor(cnv$chr, levels = unique(cnv$chr))
cnv$pos = rowMeans(cnv[,c("start","end")])
cnv$type = factor(cnv$type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain"))

p_cnv = ggplot(cnv, aes(x=pos, y=ratio, color = type)) +
        geom_point(alpha=0.8) +
        facet_grid(~chr, scale = "free", space = "free") +
        scale_color_d3(guide = guide_legend(override.aes = list(size = 5))) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = c(0.2, 0.2),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())

p_cnv_filtered = ggplot(cnv[cnv$ratio < -0.5 | cnv$ratio > 0.3, ], aes(x=pos, y=ratio, color = type)) +
        geom_point(alpha=0.8) +
        facet_grid(~chr, scale = "free", space = "free") +
        scale_color_d3(guide = guide_legend(override.aes = list(size = 5))) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = c(0.2, 0.2),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())

load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_cytoband.RData")
bed_cytoband_1 = cbind(as.data.frame(bed_cytoband_1), cytoband_bed[match(bed_cytoband_1$V8, cytoband_bed$cytoband),])
bed_cytoband_1 = bed_cytoband_1[bed_cytoband_1$chrom!="X", ]
bed_cytoband_1$chrom = factor(bed_cytoband_1$chrom, levels = 1:22)
bed_cytoband_1$pos = rowMeans(bed_cytoband_1[,c("start", "end")])
bed_cytoband_1$type = factor(bed_cytoband_1$type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain"))

bed_cytoband_1_filterd = lapply(1, function(x){
    load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_cytoband_filtered.RData")
    return(bed_cytoband_1)
    })
bed_cytoband_1_filterd = bed_cytoband_1_filterd[[1]]
bed_cytoband_1_filterd = cbind(as.data.frame(bed_cytoband_1_filterd), cytoband_bed[match(bed_cytoband_1_filterd$V8, cytoband_bed$cytoband),])
bed_cytoband_1_filterd = bed_cytoband_1_filterd[bed_cytoband_1_filterd$chrom!="X", ]
bed_cytoband_1_filterd$chrom = factor(bed_cytoband_1_filterd$chrom, levels = 1:22)
bed_cytoband_1_filterd$pos = rowMeans(bed_cytoband_1_filterd[,c("start", "end")])
bed_cytoband_1_filterd$type = factor(bed_cytoband_1_filterd$type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain"))

p_cnv_cytoband = ggplot(bed_cytoband_1, aes(x=pos, y=ratio, color = type)) +
        geom_point(alpha=0.8) +
        labs(color = "") +
        facet_grid(~chrom, scale = "free", space = "free") +
        scale_color_d3(guide = guide_legend(override.aes = list(size = 5))) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = c(0.2, 0.2),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())

p_cnv_cytoband_filtered = ggplot(bed_cytoband_1_filterd, aes(x=pos, y=ratio, color = type)) +
        geom_point(alpha=0.8) +
        labs(color = "") +
        facet_grid(~chrom, scale = "free", space = "free") +
        scale_color_d3(guide = guide_legend(override.aes = list(size = 5))) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = c(0.2, 0.2),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())

p_cnv_icd = plot_all_icd_cnv(res[!grepl("X", res$cyto),], title = "CNA-ICD10", maxP = 5e-8, legend_pos = c(0.3, 0.8))

pdf("man_cytoband_ICD10.pdf", width=15, height=6)
cowplot::plot_grid(p_cnv, p_cnv_icd, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.8, 1))
cowplot::plot_grid(p_cnv_cytoband, p_cnv_icd, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.8, 1))
dev.off()

p_cnv_filterd_icd = plot_all_icd_cnv_filterd(res_filterd[!grepl("X", res_filterd$cyto),], title = "CNA-ICD10", maxP = 1e-5, legend_pos = "bottom")

pdf("man_cytoband_filtered_ICD10.pdf", width=15, height=6)
cowplot::plot_grid(p_cnv_filtered, p_cnv_filterd_icd, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.3, 1))
cowplot::plot_grid(p_cnv_cytoband_filtered, p_cnv_filterd_icd, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.3, 1))
dev.off()

plot_all_icd_cnv = function(res_sig, title, maxP = 5e-8, legend_pos = "bottom"){
    res_sig$OR1 = res_sig$OR
    res_sig$OR1[res_sig$OR >= 60] = ">=60"
    res_sig$OR1[res_sig$OR >= 30 & res_sig$OR < 60] = "[30,60)"
    res_sig$OR1[res_sig$OR >= 15 & res_sig$OR < 30] = "[15,30)"
    res_sig$OR1[res_sig$OR >= 10 & res_sig$OR < 15] = "[10,15)"
    res_sig$OR1[res_sig$OR >= 5 & res_sig$OR < 10] = "[5,10)"
    res_sig$OR1[res_sig$OR > 1 & res_sig$OR < 5] = "(1,5)"
    res_sig$OR1[res_sig$OR > 0 & res_sig$OR < 1] = "(0,1)"
    res_sig$cyto = factor(res_sig$cyto, levels = cytoband_bed$cytoband)
    res_sig$pos = rowMeans(res_sig[,c("start","end")])
    res_sig$CHR = factor(as.numeric(gsub("Y", "24", gsub("X", 23, res_sig$chrom))), levels = 1:24)
    res_sig$logP = -log10(res_sig$P) * res_sig$cnv
    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
        "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
    res_sig$icd[res_sig$P >= maxP] = NA
    p = ggplot(res_sig, aes(x=pos, y=logP, color=icd, size = OR1)) +
        geom_point(alpha=0.8) +
        labs(color = "ICD10", size = "Odds ratio") +
        geom_hline(yintercept = -log10(maxP), linetype = "dashed") +
        geom_hline(yintercept = log10(maxP), linetype = "dashed") +
        facet_grid(~CHR, scale = "free", space = "free") +
        # scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 5)), na.value = "grey") +
        scale_color_d3(guide = guide_legend(nrow = 1, override.aes = list(size = 5)), na.value = "grey") +
        # scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_size_manual(breaks = c("(0,1)", "(1,5)", "(5,10)", "[10,15)", "[15,30)", "[30,60)", ">=60"), values = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
        # scale_y_break(c(20.5, 41), scale=0.3, space = 0) +
        # scale_y_cut(20) +
        # scale_y_sqrt() +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = legend_pos,
            legend.direction = "horizontal",
            legend.box = "horizontal",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())
    return(p)
}

plot_all_icd_cnv_filterd = function(res_sig, title, maxP = 5e-8, legend_pos = "bottom"){
    res_sig$OR1 = res_sig$OR
    res_sig$OR1[res_sig$OR >= 30] = ">=30"
    res_sig$OR1[res_sig$OR >= 15 & res_sig$OR < 30] = "[15,30)"
    res_sig$OR1[res_sig$OR >= 10 & res_sig$OR < 15] = "[10,15)"
    res_sig$OR1[res_sig$OR >= 5 & res_sig$OR < 10] = "[5,10)"
    res_sig$OR1[res_sig$OR > 1 & res_sig$OR < 5] = "(1,5)"
    res_sig$OR1[res_sig$OR > 0 & res_sig$OR < 1] = "(0,1)"
    res_sig$cyto = factor(res_sig$cyto, levels = cytoband_bed$cytoband)
    res_sig$pos = rowMeans(res_sig[,c("start","end")])
    res_sig$CHR = factor(as.numeric(gsub("Y", "24", gsub("X", 23, res_sig$chrom))), levels = 1:24)
    res_sig$logP = -log10(res_sig$P) * res_sig$cnv
    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
        "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
    res_sig$icd[res_sig$P >= maxP] = NA
    res_sig$lab = res_sig$icd
    res_sig$lab[is.na(res_sig$lab)] = ""
    # res_sig$ICDanno[res_sig$P >= maxP] = NA
    p = ggplot(res_sig, aes(x=pos, y=logP, color=ICDanno, size = OR1)) +
        geom_point(alpha=0.7) +
        geom_text_repel(aes(label = lab), size = 3) +
        labs(color = "ICD10", size = "Odds ratio") +
        geom_hline(yintercept = -log10(maxP), linetype = "dashed") +
        geom_hline(yintercept = log10(maxP), linetype = "dashed") +
        facet_grid(~CHR, scale = "free", space = "free") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 4, override.aes = list(size = 5)), na.value = "grey") +
        # scale_color_d3(guide = guide_legend(nrow = 1, override.aes = list(size = 5)), na.value = "grey") +
        # scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_y_continuous(breaks= c(-7,-5,-3,3,5,7)) +
        scale_size_manual(breaks = c("(0,1)", "(1,5)", "(5,10)", "[10,15)", "[15,30)", ">=30"), values = c(0.5, 1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(20.5, 41), scale=0.3, space = 0) +
        # scale_y_cut(20) +
        # scale_y_sqrt() +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = legend_pos,
            legend.direction = "horizontal",
            legend.box = "vertical",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())
    return(p)
}
```

### board

```R
#data6
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_board.RData")
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge_Model/icd_status_PC_all_BMI.RData")
clin_status$E66 = NULL
clin_status$E66 = ifelse(clin_status$BMI > 24, 1, 0)
clin_status$BMI = NULL
clin_status = clin_status[, grep("\\.1", colnames(clin_status),invert=T)]
clin_status = clin_status[, grep("R", colnames(clin_status),invert=T)]
clin_status = clin_status[, c(colnames(clin_status)[1:7], names(which(colSums(clin_status[,-c(1:7)], na.rm = T)>30)))]

mat_board = mat_board[match(rownames(clin_status), rownames(mat_board)),]
rownames(mat_board) = rownames(clin_status)
mat_board[is.na(mat_board)] = 0

library(parallel)
res = mclapply(colnames(mat_board), function(x){
    res = lapply(colnames(clin_status)[-c(1:7)], function(y){
        print(paste(x, y))
        dat = na.omit(cbind(clin_status[,c(y, "age","PC1","PC2","PC3","PC4","PC5")], cnv = mat_board[,x]))
        if(length(unique(dat$cnv))>1){
            colnames(dat)[1]='icd'
            res_glm = coef(summary(glm(icd ~ cnv + age + PC1 + PC2 + PC3 + PC4 + PC5, data = dat, family="binomial")))[, c(1,4)]
            res_glm = as.data.frame(res_glm[c(1, grep("cnv", rownames(res_glm))),], stringsAsFactors=F)
            colnames(res_glm) = c("Est", "P")
            res_glm$icd = y
            res_glm$cnv = rownames(res_glm)
            return(res_glm[-1, ])
        }else{
            return(NA)
        }
    })
    res = do.call(rbind, res[which(sapply(res, length)>1)])
    if(!is.null(res)){
        res$cyto = x
        res$cnv = gsub("cnv", "", res$cnv)
    }else{
        res = NA
    }
    return(res)
}, mc.cores = 35)
res = do.call(rbind, res)
res$Est = as.numeric(res$Est)
res$P = as.numeric(res$P)
rownames(res) = NULL
res = na.omit(res)
write.table(res, file="pheWAS_status_PC_cnv_board.txt", r=F,c=T,sep="\t",quote=F)

###############################

library(dplyr)
library(ggplot2)
library(ggsci)
library(ggbreak)
library(ggrepel)

cytoband_bed = read.table("/share/data3/NIPT/cytoband_hg38.txt", header=F,sep="\t",stringsAsFactors=F)
colnames(cytoband_bed) = c("chrom", "start", "end", "cytoband", "V5")
board_bed = cytoband_bed
board_bed$board = ""
for(i in c(paste0(c(1:22, "X","Y"), "p"), paste0(c(1:22, "X","Y"), "q"))){
    board_bed$board[grep(paste0("^", i), board_bed$cytoband)] = i
}
board_bed = board_bed %>% group_by(board) %>% summarise(chrom=unique(chrom), start = min(start), end = max(end))
board_bed = board_bed[, c(2:4,1)]

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
                            levels = c("Pregnancy, childbirth and the puerperium", "Genitourinary system", "Endocrine, nutritional and metabolic diseases", "Digestive system", "Respiratory system", "Skin and subcutaneous tissue", "Musculoskeletal system and connective tissue",
                                "Certain infectious and parasitic diseases", "Eye and adnexa", "Ear and mastoid process", "Neoplasms", "Circulatory system", "Mental and behavioural disorders", "Blood and blood-forming organs and certain disorders involving the immune mechanism",
                                "Congenital malformations, deformations and chromosomal abnormalities"))


p_cnv_board = bed_board %>% filter(fraction>0.8, !grepl("X", board)) %>% mutate(board = factor(board, levels = unlist(lapply(1:22, function(x){paste0(x, c("p","q"))})))) %>%
    mutate(type = factor(type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain"))) %>%
    ggplot(aes(x=board, y=abs(ratio), color = type)) +
        geom_boxplot() +
        scale_color_d3(guide = guide_legend(override.aes = list(size = 5))) +
        theme_bw() +
        theme(
            axis.text = element_text(size=10),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = c(0.2, 0.8),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())

p_cnv_icd = plot_all_icd_cnv(res[res$P<0.05 & res$icd!="O04", ], title = "CNA-ICD10", maxP = 1e-4, legend_pos = "bottom")

pdf("man_board_ICD10.pdf", width=15, height=6)
cowplot::plot_grid(p_cnv_board, p_cnv_icd, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.8, 2))
dev.off()

#  bed_board %>% filter(fraction>0.8) %>%
#     mutate(board = factor(board, levels = unlist(lapply(1:22, function(x){paste0(x, c("p","q"))}))),
#         type = factor(type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain"))) 

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

p_cnv_icd_filtered = plot_all_icd_cnv_filterd(res[res$P<0.05 & res$icd !="O04" & res$cyto %in% names(which(colSums(mat_board!=0)>10)), ], 
    title = "CNA-ICD10", maxP = 1e-4, legend_pos = "bottom")

# p_cnv_icd_filtered_2 = plot_all_icd_cnv_filterd(res[res$P<0.05 & res$cyto %in% names(which(colSums(mat_board!=0)>10)), ], 
#     title = "CNA-ICD10", maxP = 1e-4, legend_pos = "bottom")

pdf("man_board_filtered_ICD10.pdf", width=12, height=8)
cowplot::plot_grid(p_cnv_board_filtered, p_cnv_icd_filtered, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.8, 2))
dev.off()

plot_all_icd_cnv = function(res_sig, title, maxP = 5e-8, legend_pos = "bottom"){
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
    # res_sig$ICDanno[res_sig$P >= maxP] = NA
    p = res_sig %>% group_by(cyto) %>% arrange(logP) %>% mutate(index=1:length(logP)) %>%
        ggplot(aes(x=index, y=logP, color=ICDanno, size = OR1)) +
        geom_point(alpha=0.7) +
        facet_grid(~cyto, scale="free")+
        geom_text_repel(aes(label = lab), size = 3) +
        labs(color = "ICD10", size = "Odds ratio", y = "Deletion --- Amplication\n\n-log10(Pvalue)") +
        geom_hline(yintercept = 0) +
        # geom_hline(yintercept = log10(maxP), linetype = "dashed") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 4, override.aes = list(size = 5)), na.value = "grey") +
        # scale_color_d3(guide = guide_legend(nrow = 1, override.aes = list(size = 5)), na.value = "grey") +
        # scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_y_continuous(breaks= c(-5,0,5,10), labels = c(5, 0, 5, ">10")) +
        scale_size_manual(breaks = c("(0,1)", "(1,5)", "(5,10)", "[10,15)", "[15,30)", ">=30"), values = c(0.5, 1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(23.5, 53), scale=0.3, space = 0) +
        # scale_y_cut(20) +
        # scale_y_sqrt() +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = legend_pos,
            legend.direction = "horizontal",
            legend.box = "vertical",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())
    return(p)
}

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
    # res_sig$ICDanno[res_sig$P >= maxP] = NA
    p = res_sig %>% group_by(cyto) %>% arrange(logP) %>% mutate(index=1:length(logP)) %>%
        ggplot(aes(x=index, y=logP, color=ICDanno, size = OR1)) +
        geom_point(alpha=0.7) +
        facet_grid(~cyto, scale="free")+
        geom_text_repel(aes(label = lab), size = 5) +
        labs(color = "", size = "Odds ratio", y = "Deletion-Amplication\n\n-log10(Pvalue)") +
        geom_hline(yintercept = 0) +
        # geom_hline(yintercept = log10(maxP), linetype = "dashed") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 3)), na.value = "grey") +
        # scale_color_d3(guide = guide_legend(nrow = 1, override.aes = list(size = 5)), na.value = "grey") +
        # scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_y_continuous(breaks= c(-5,0,5,10), labels = c(5, 0, 5, ">10")) +
        scale_size_manual(breaks = c("(0,1)", "(1,5)", "(5,10)", "[10,15)", "[15,30)", ">=30"), values = c(0.5, 1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(23.5, 53), scale=0.3, space = 0) +
        # scale_y_cut(20) +
        # scale_y_sqrt() +
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
            # legend.position= 'bottom',
            legend.position = legend_pos,
            legend.direction = "horizontal",
            legend.box = "vertical",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())
    return(p)
}

######## 长CNV

cnv_1 = bed_board %>% filter(fraction > 0.8, !grepl("X", board)) %>% select(chr, start, end, ratio)
cnv_1$chr = paste0("chr", cnv_1$chr)

library(circlize)
library(ggsci)
library(Cairo)

# CairoPNG("test_1.png", width=600,height=600)
pdf("cnv_long.pdf", width = 6, heigh = 6)
circos.par(gap.after=c(rep(1,21), 20), "start.degree" = 90, track.height = 0.3)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), plotType = c("ideogram", "labels"), species='hg38')
circos.genomicTrackPlotRegion(cnv_1, ylim = c(-4, 3), panel.fun = function(region, value, ...) {
    col = ifelse(value[[1]] > 0, "#FF7F0EFF", "#1F77B4FF")
    circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16) 
}, track.height = 0.5)
circos.yaxis(side = "left", at = seq(-4, 3, by = 0.5),
             sector.index = get.all.sector.index()[1], labels.cex = 0.5)
circos.clear()
text(-0.15, 0.65, "ratio", cex = 1, srt = 95)
legend("center", title = "", pch = 16, cex = 1, bty = "n",
       col = c("#FF7F0EFF", "#1F77B4FF"),
       legend = c("Amplication", "Deletion"))
dev.off()

######## 短CNV

bed_board_2 = bed_board %>% filter(fraction < 0.2, !grepl("X", board), n > 10)

cnv_2 = bed_board_2[, c("chr", "start", "end", "ratio")]
cnv_2$chr = paste0("chr", cnv_2$chr)

library(circlize)
library(ggsci)
library(Cairo)

# CairoPNG("test_1.png", width=600,height=600)
pdf("cnv_short_n10.pdf", width = 6, heigh = 6)
circos.par(gap.after=c(rep(1,21), 20), "start.degree" = 90, track.height = 0.3)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), plotType = c("ideogram", "labels"), species='hg38')
circos.genomicTrackPlotRegion(cnv_2, ylim = c(-4, 3), panel.fun = function(region, value, ...) {
    col = ifelse(value[[1]] > 0, "#FF7F0EFF", "#1F77B4FF")
    circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16) 
}, track.height = 0.5)
circos.yaxis(side = "left", at = seq(-4, 3, by = 0.5),
             sector.index = get.all.sector.index()[1], labels.cex = 0.5)
circos.clear()
text(-0.15, 0.65, "ratio", cex = 1, srt = 95)
legend("center", title = "", pch = 16, cex = 1, bty = "n",
       col = c("#FF7F0EFF", "#1F77B4FF"),
       legend = c("Amplication", "Deletion"))
dev.off()

bed_board_2 %>% group_by(region) %>%
    mutate(pos = median(c(start, end)),
        board = factor(board, levels = unlist(lapply(1:22, function(x){paste0(x, c("p","q"))}))),
        type = factor(type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain")))

p_cnv_board_filtered_short = bed_board_2 %>% mutate(
        pos = median(c(start, end)),
        board = factor(board, levels = unlist(lapply(1:22, function(x){paste0(x, c("p","q"))}))),
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

bed_board_2 %>% group_by(region) %>%
    mutate(pos = median(c(start, end)),
        board = factor(board, levels = unlist(lapply(1:22, function(x){paste0(x, c("p","q"))}))),
        type = factor(type, labels = c("Deletion", "Amplication"), levels = c("loss", "gain"))) %>%
    arrange(ratio) %>% mutate(index = 1:length(ratio)) %>%
    ggplot(aes(x=pos, y=ratio, color = type)) +
        geom_point() + labs(x="", y="log2(Ratio)", color="") +
        facet_grid(~chrom, scale="free") +
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

p_cnv_icd_filtered = plot_all_icd_cnv_filterd(res[res$P<0.05 & res$icd !="O04" & res$cyto %in% names(which(colSums(mat_board!=0)>10)), ], 
    title = "CNA-ICD10", maxP = 1e-4, legend_pos = "bottom")

pdf("man_board_filtered_ICD10.1.pdf", width=12, height=8)
cowplot::plot_grid(p_cnv_board_filtered, p_cnv_icd_filtered, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.8, 2))
dev.off()
```

## child

```R
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_cnv/mat_cytoband_filtered.RData")
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge/avsnp150_ChainMAP_genos/all/com_0.01/SNPs_RD_AD_cohort1_cohort2/ICD_FF_status_num10.RData")
load("/share/data3/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge_Model/icd_status_PC_all_BMI.RData")

int = intersect(ICD_mat_status$ID, rownames(mat_cytoband))
ICD_mat_status = ICD_mat_status[ICD_mat_status$ID %in% int, ]
mat_cytoband = mat_cytoband[ICD_mat_status$ID, ]

ICD_mat_status = ICD_mat_status[, c(colnames(ICD_mat_status)[1:4], names(which(colSums(ICD_mat_status[,-c(1:4)], na.rm = T)>30)))]

library(parallel)
res = mclapply(colnames(mat_cytoband), function(x){
    res = lapply(colnames(ICD_mat_status)[-c(1:7)], function(y){
        print(paste(x, y))
        dat = na.omit(cbind(ICD_mat_status[,c(y, "follow_year","gender","FF")], cnv = mat_cytoband[,x]))
        if(length(unique(dat$cnv))>1){
            colnames(dat)[1]='icd'
            res_glm = coef(summary(glm(icd ~ cnv + follow_year + gender + FF, data = dat, family="binomial")))[, c(1,4)]
            res_glm = as.data.frame(res_glm[c(1, grep("cnv", rownames(res_glm))),], stringsAsFactors=F)
            colnames(res_glm) = c("Est", "P")
            res_glm$icd = y
            res_glm$cnv = rownames(res_glm)
            return(res_glm[-1, ])
        }else{
            return(NA)
        }
    })
    res = do.call(rbind, res[which(sapply(res, length)>1)])
    if(!is.null(res)){
        res$cyto = x
        res$cnv = gsub("cnv", "", res$cnv)
    }else{
        res = NA
    }
    return(res)
}, mc.cores = 35)
res = do.call(rbind, res)
res$Est = as.numeric(res$Est)
res$P = as.numeric(res$P)
rownames(res) = NULL
res = na.omit(res)
write.table(res, file="child/pheWAS_status_PC_cnv_filtered_cytoband.txt", r=F, c=T, sep="\t", quote=F)


res_filterd = res
res_filterd = merge(res_filterd, cytoband_bed[,1:4], by.x = "cyto", by.y = c("cytoband"))
res_filterd$OR = exp(res_filterd$Est)
res_filterd$cnv = as.numeric(res_filterd$cnv)
res_filterd$ICDanno = ""
res_filterd$ICDanno[grep("[AB]", res_filterd$icd)] = "Certain infectious and parasitic diseases"
res_filterd$ICDanno[which(res_filterd$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24", "D25", "D26"))] = "Neoplasms"
res_filterd$ICDanno[which(res_filterd$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
res_filterd$ICDanno[grep("E", res_filterd$icd)] = "Endocrine, nutritional and metabolic diseases"
res_filterd$ICDanno[grep("F", res_filterd$icd)] = "Mental and behavioural disorders"
res_filterd$ICDanno[which(res_filterd$icd %in% c("H00", "H01", "H02", "H04", "H10","H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
res_filterd$ICDanno[which(res_filterd$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92","H93"))] = "Ear and mastoid process"
res_filterd$ICDanno[grep("I", res_filterd$icd)] = "Circulatory system"
res_filterd$ICDanno[grep("J", res_filterd$icd)] = "Respiratory system"
res_filterd$ICDanno[grep("K", res_filterd$icd)] = "Digestive system"
res_filterd$ICDanno[grep("L", res_filterd$icd)] = "Skin and subcutaneous tissue"
res_filterd$ICDanno[grep("M", res_filterd$icd)] = "Musculoskeletal system and connective tissue"
res_filterd$ICDanno[grep("N", res_filterd$icd)] = "Genitourinary system"
res_filterd$ICDanno[grep("O", res_filterd$icd)] = "Pregnancy, childbirth and the puerperium"
res_filterd$ICDanno[grep("P", res_filterd$icd)] = "Certain conditions originating in the perinatal period"
res_filterd$ICDanno[grep("Q", res_filterd$icd)] = "Congenital malformations, deformations and chromosomal abnormalities"
res_filterd = res_filterd[res_filterd$ICDanno !="",]
res_filterd$ICDanno = factor(res_filterd$ICDanno,
                            levels = c("Pregnancy, childbirth and the puerperium", "Genitourinary system", "Endocrine, nutritional and metabolic diseases", "Digestive system", "Respiratory system", "Skin and subcutaneous tissue", "Musculoskeletal system and connective tissue",
                                "Certain infectious and parasitic diseases", "Eye and adnexa", "Ear and mastoid process", "Neoplasms", "Circulatory system", "Mental and behavioural disorders", "Blood and blood-forming organs and certain disorders involving the immune mechanism",
                                "Congenital malformations, deformations and chromosomal abnormalities", "Certain conditions originating in the perinatal period"))

p_cnv_filterd_icd = plot_all_icd_cnv_filterd(res_filterd[!grepl("X", res_filterd$cyto),], title = "child_CNA-ICD10", maxP = 1e-3, legend_pos = "bottom")

pdf("man_cytoband_filtered_ICD10_child.pdf", width=15, height=6)
p_cnv_filterd_icd
# cowplot::plot_grid(p_cnv_filtered, p_cnv_filterd_icd, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.3, 1))
# cowplot::plot_grid(p_cnv_cytoband_filtered, p_cnv_filterd_icd, ncol = 1, axis = 'lr', align = "v", rel_heights=c(0.3, 1))
dev.off()

plot_all_icd_cnv_filterd = function(res_sig, title, maxP = 5e-8, legend_pos = "bottom"){
    res_sig$OR1 = res_sig$OR
    res_sig$OR1[res_sig$OR >= 30] = ">=30"
    res_sig$OR1[res_sig$OR >= 15 & res_sig$OR < 30] = "[15,30)"
    res_sig$OR1[res_sig$OR >= 10 & res_sig$OR < 15] = "[10,15)"
    res_sig$OR1[res_sig$OR >= 5 & res_sig$OR < 10] = "[5,10)"
    res_sig$OR1[res_sig$OR > 1 & res_sig$OR < 5] = "(1,5)"
    res_sig$OR1[res_sig$OR > 0 & res_sig$OR < 1] = "(0,1)"
    res_sig$cyto = factor(res_sig$cyto, levels = cytoband_bed$cytoband)
    res_sig$pos = rowMeans(res_sig[,c("start","end")])
    res_sig$CHR = factor(as.numeric(gsub("Y", "24", gsub("X", 23, res_sig$chrom))), levels = 1:24)
    res_sig$logP = -log10(res_sig$P) * res_sig$cnv
    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
        "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
    res_sig$icd[res_sig$P >= maxP] = NA
    res_sig$lab = res_sig$icd
    res_sig$lab[is.na(res_sig$lab)] = ""
    # res_sig$ICDanno[res_sig$P >= maxP] = NA
    p = ggplot(res_sig, aes(x=pos, y=logP, color=ICDanno, size = OR1)) +
        geom_point(alpha=0.7) +
        geom_text_repel(aes(label = lab), size = 3) +
        labs(color = "ICD10", size = "Odds ratio") +
        geom_hline(yintercept = -log10(maxP), linetype = "dashed") +
        geom_hline(yintercept = log10(maxP), linetype = "dashed") +
        facet_grid(~CHR, scale = "free", space = "free") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 4, override.aes = list(size = 5)), na.value = "grey") +
        # scale_color_d3(guide = guide_legend(nrow = 1, override.aes = list(size = 5)), na.value = "grey") +
        # scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_y_continuous(breaks= c(-5,-3,-1,1,3,5)) +
        scale_size_manual(breaks = c("(0,1)", "(1,5)", "(5,10)", "[10,15)", "[15,30)", ">=30"), values = c(0.5, 1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(20.5, 41), scale=0.3, space = 0) +
        # scale_y_cut(20) +
        # scale_y_sqrt() +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=10),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            # legend.position= 'bottom',
            legend.position = legend_pos,
            legend.direction = "horizontal",
            legend.box = "vertical",
            panel.spacing=unit(0,"line"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.grid.major.x = element_blank())
    return(p)
}
```
