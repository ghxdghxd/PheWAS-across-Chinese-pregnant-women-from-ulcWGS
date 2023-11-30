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
load("icd_status_PC_all.RData")
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
    --covarColList=age,PC1,PC2,PC3,PC4,PC5 \
    --sampleIDColinphenoFile=ID \
    --traitType=binary \
    --outputPrefix=./output1/$icd \
    --nThreads=30 \
    --LOCO=TRUE \
    --maxMissingRateforGRM=0.95 \
    --minMAFforGRM=0 \
    --IsOverwriteVarianceRatioFile=TRUE

step2_SPAtests.R \
    --vcfFile=$j \
    --vcfFileIndex=$j.csi \
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

load("icd_status_PC_all.RData")
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
load("icd_status_PC_all.RData")
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
save(snp_power, file="phewas_power_P0.001.RData")
```

## 7.2 PheWAS for maternal CNVs

### 7.2.1 CNVs of chromosome boards

```R
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
colnames(bed_board) = c("chr", 'start', 'end', 'region', "type", "n", "chrom", "board_start","board_end","board")
bed_board$fraction = c(bed_board$end-bed_board$start+1)/c(bed_board$board_end-bed_board$board_start)
bed_board = merge(bed_cnv, bed_board, by=c("chr", 'start', 'end', 'region', "type"))

load("clin.RData")

bed_board_1 = bed_board %>% filter(fraction > 0.8, !grepl("X", board)) %>%
    group_by(ID, board) %>% summarise(type=paste(sort(unique(type)), collapse=","), ratio = mean(ratio), zscore = mean(zscore))

mat_board = reshape2::acast(bed_board_1[,1:3], ID ~ board)

mat_board[mat_board=="gain"] = 1
mat_board[mat_board=="loss"] = -1
mat_board[is.na(mat_board)] = 0
save(mat_board, bed_board_1, file="mat_board.RData")
```

### 7.2.1 CNVs associated with maternal phenotypes

```R
#data6
load("mat_board.RData")
load("icd_status_PC_all.RData")
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
```
