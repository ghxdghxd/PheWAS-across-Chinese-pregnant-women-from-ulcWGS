# Genotyping

## WildType germline variants of each sample

```sh
bcftools query -i 'ALT=="." & GT == "0/0"' -f '%CHROM\t%POS\t%REF\t%DP4{0}\t%DP4{1}\n' $sampleID.b.vcf.gz | \
    awk 'BEGIN{DP=0;s["+"]="-";s["-"]="+";ss="-"}
        {OFS="\t";if(DP!=$4+$5){DP=$4+$5;ss=s[ss]};print $1,$2-1,$2,"",DP,ss}' | \
    bedtools merge -s -c 5 -o distinct | bgzip > $sampleID.b.wt.gz

intersectBed -wb -a all.avsnp150.ChinaMAP.AF0.rm_um35.lite.vcf.gz -b $sampleID.b.wt.gz |\
    cut -f 1-2,4-5,12|bgzip -@2 > $sampleID.avsnp150.ChinaMAP.wt.gz
```

## Mutant germline varinats of each sample

```sh
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP4{0}\t%DP4{1}\t%DP4{2}\t%DP4{3}\t[%GT]\n' $sampleID.sorted.dedup.bqsr.b.var.norm.vcf.gz | \
    awk '{split($9,a,"/");OFS="\t";print $1,$2,$3,$4,$5+$6,$7+$8,a[1]+a[2]}' > $sampleID.avsnp150.ChinaMAP.var
bgzip $sampleID.avsnp150.ChinaMAP.var
```

## combining all genotypes

### for each sample

```R
library(data.table)
get_geno = function(x){
    wt = fread(wt_files[x],header=F,sep="\t",stringsAsFactors=F)
    wt = wt[, min(V5), by = c("V1","V2","V3","V4")]
    colnames(wt)[5] = "V5"
    setkey(wt, V1, V2, V3, V4)
    wt$V6=0
    wt$V7=0
    var = fread(var_files[x],header=F,sep="\t",stringsAsFactors=F)
    setkey(var, V1, V2, V3, V4)
    geno = rbind(wt[!var], var)
    colnames(geno) = c("chr","pos","ref","alt","RD", "AD", "GT")
    setkey(geno, chr, pos, ref, alt)
    fwrite(bed[geno, nomatch=0], file = out_files[x], row.names = F, col.names = T, sep="\t",quote=F)
}

wt_files = list.files("wt_genos", pattern="*.gz", full.names = T)
var_files = gsub("wt", "var", wt_files)
out_files = gsub("wt", "geno", basename(wt_files))

bed = fread("all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.lite.bed.gz", header = F, sep="\t", stringsAsFactors = F)
colnames(bed) = c("chr","pos","ref","alt","cytoband")
setkey(bed, chr,pos,ref,alt)
bed$cytoband=NULL

library(foreach)
library(doParallel)

registerDoParallel(20)
foreach(a = 1:length(wt_files)) %dopar% {
    get_geno(a)
}
```

### for all samples

```R
library(data.table)
bed = fread("all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.lite.bed.gz", header = T, sep = "\t", stringsAsFactors=F)
bed = bed[,c(1:5)]
colnames(bed) = c("chr", "pos", "ref", "alt", "cytoband")
setkey(bed, chr,pos,ref,alt)
# 1p31.1

files = list.files("all_wt_var_genos", pattern = "*.geno.gz", full.names=T)

library(stringr)
library(parallel)

cytoband = unique(bed$cytoband)
index_1 = names(which(table(bed$cytoband)==1))

cytoband_list = lapply(unique(bed$cytoband), function(x){
    if(x == "9q11"){
        return(c("9q11", "9q12"))
    }else if(x %in% c(index_1)){
        return(NULL)
    }else{
        return(c(x))
    }
})
cytoband_list = cytoband_list[!sapply(cytoband_list, is.null)]

get_geno = function(fs){
    geno = mclapply(fs, function(f){
        geno <- fread(f, header=T,sep="\t",stringsAsFactors=F)
        geno = unique(geno)
        setkey(geno, chr,pos,ref,alt)
        geno = geno[bed]
        gc()
        return(geno[, 5:7])
    }, mc.cores = 10)
    RD = do.call(rbind, lapply(geno, function(x){
        t(data.table(x$RD))
    }))
    AD = do.call(rbind, lapply(geno, function(x){
        t(data.table(x$AD))
    }))
    rm(geno)
    gc()
    write.table(gsub(".avsnp150.ChinaMAP.geno.gz","", basename(fs)), file = "sample_ID1.txt", row.names = F, col.names = F, append = T, sep = "\t", quote = F)
    res = mclapply(cytoband_list, function(j){
        fwrite(RD[, which(bed$cytoband %in% j)], file = paste0("RD/", paste(j, collapse="_"), ".RD.txt.gz"), row.names = F, col.names = F, append = T, sep = "\t", quote = F)
        fwrite(AD[, which(bed$cytoband %in% j)], file = paste0("AD/", paste(j, collapse="_"), ".AD.txt.gz"), row.names = F, col.names = F, append = T, sep = "\t", quote = F)
    }, mc.cores = 38)
}

indexs = parallel::splitIndices(length(files), length(files)/50)
for(i in indexs){
    get_geno(files[i])
    gc()
}
```

#### genotyping for all samples

```R
library(data.table)
library(dplyr)
library(SNPassoc)
library(parallel)

bed <- fread("all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.lite.bed.gz", header=T, sep="\t",stringsAsFactors=F)

files = list.files("/data/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge/avsnp150_ChainMAP_genos/all/raw/GT", pattern = "*.gz", full.names = T)

get_AF_callNum_HWE_GT = function(f){
    bed_sub = bed[bed$cytoband %in% strsplit(gsub(".GT.txt.gz", "", basename(f)), split="_")[[1]], ]
    index = which(bed_sub$AF > 0.01 & bed_sub$AF < 0.99)
    if(length(index)==0){
        return(NULL)
    }
    bed_sub = bed_sub[index]
    RD = fread(gsub("GT","RD",f), header = F, sep = "\t", stringsAsFactors = F)
    RD = RD[, ..index]
    AD = fread(gsub("GT","AD",f), header = F, sep = "\t", stringsAsFactors = F)
    AD = AD[, ..index]
    GT = AD
    GT[RD>0 & AD==0] = 0
    GT[RD>0 & AD>0] = 1
    GT[RD==0 & AD>0] = 2
    rm(RD)
    rm(AD)
    gc()
    bed_sub$callNum = colSums(!is.na(GT))
    # Binomial distribution
    bed_sub$AF_sd = sqrt((bed_sub$AF*(1-bed_sub$AF))/bed_sub$callNum)
    AF_callNum = foreach(n = 1:ncol(GT)) %dopar% {
        res = summary(snp(unlist(GT[,..n]), name.genotypes=c(0,1,2)))
        return(c(obAF = res$allele.freq["B","percentage"]))
    }
    AF_callNum = as.data.frame(do.call(rbind, AF_callNum), stringsAsFactors = F)
    AF_callNum[, "obAF"] = AF_callNum[,"obAF"]/100
    bed_sub = cbind(bed_sub, AF_callNum)
    fwrite(GT, file=paste0("GT/", basename(f)), row.names=F, col.names=F, sep="\t")
    fwrite(bed_sub, file=paste0("SNPs/", gsub("GT", "SNPs", basename(f))), row.names = F, col.names = T, sep = "\t", quote = F)
}


for(i in rev(files)){
    if(!file.exists(paste0("GT/", basename(i)))){
        print(i)
        get_AF_callNum_HWE_GT(i)
    }
}
```
