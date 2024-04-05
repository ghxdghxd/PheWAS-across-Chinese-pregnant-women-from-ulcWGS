# 4 Ancestry analysis

## 4.1 filter germline variants

```R
library(data.table)
library(parallel)
bed_files = list.files("SNPs", pattern="*.gz", full.names=T)
bed = mclapply(bed_files, function(x){
    bed <- fread(x, header=T,sep="\t",stringsAsFactors=F)
    bed$cytoband = gsub(".SNPs.txt.gz", "", basename(x))
    return(bed)
}, mc.cores = 20)
bed = do.call(rbind, bed)

min_CallNum = quantile(bed$callNum, probs = 0.99)
bed = bed[which((bed$obAF >= bed$AF - 1.96*bed$AF_sd & bed$obAF <= bed$AF + 1.96*bed$AF_sd) & bed$AF > 0.1 & bed$AF < 0.9 & bed$callNum > min_CallNum)]
write.table(bed$id, file="all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.com0.1.callRate0.99.conf95.snpID.txt", row.names = F, col.names = F, sep="\t",quote=F)

geno = mclapply(unique(bed$cytoband), function(x){
    SNPs = bed[bed$cytoband == x]
    GT = fread(paste0("GT/", x, ".GT.txt.gz"), header=F,sep="\t",stringsAsFactors=F)
    return(GT)
}, mc.cores = 20)

geno = do.call(cbind, geno)
sample_ID <- read.table("all_sample_ID.txt", header=F,sep="\t",stringsAsFactors=F)

save(sample_ID, geno, file="all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.com0.1.callRate0.99.conf95.GT1.RData")
```

## 4.2 get genotypes in 1000 Genomes Project

```sh
for i in {1..22};
do
{
    vcftools --gzvcf 1000G/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr$i.filtered.shapeit2-duohmm-phased.vcf.gz \
        --snps all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.com0.1.callRate0.99.conf95.snpID.txt --recode --out $i
} &
done

grep -v "##" 1.recode.vcf |sed 's/#//' > all.1kg.vcf
for i in {2..22};do grep -v "#" $i.recode.vcf >> all.1kg.vcf;done
rm *recode.vcf *log
```

```R
library(data.table)
geno_1000G = fread("all.1kg.vcf", header=T,sep="\t",stringsAsFactors=F)
geno_1000G = geno_1000G[, -c(1:2,4:9)]
geno_1000G[geno_1000G=="0|0"]=0
geno_1000G[geno_1000G=="0|1"]=1
geno_1000G[geno_1000G=="1|0"]=1
geno_1000G[geno_1000G=="1|1"]=2
geno_1000G = as.data.frame(geno_1000G, stringsAsFactors=F)
rownames(geno_1000G) = geno_1000G$ID
geno_1000G$ID=NULL
geno_1000G = t(geno_1000G)
ID_1000G = read.csv("1000G/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/igsr_samples.tsv", header=T,sep="\t",stringsAsFactors=F)
save(geno_1000G, ID_1000G, file="all.1kg.geno.RData")
```

## 4.3 Principal component analysis

```R
library(data.table)
load("all.1kg.geno.RData")
geno_1000G = as.data.frame(geno_1000G, stringsAsFactors=F)

load("all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.com0.1.callRate0.99.conf95.GT1.RData")
int = intersect(colnames(geno), colnames(geno_1000G))
geno = geno[, ..int]
geno = geno[!duplicated(sample_ID),]

library(SYNCSA)
g2 = rbind(geno, geno_1000G[ID_1000G$Sample.name[ID_1000G$Population.code %in% c("CHB","CHS","CDX", "JPT")],])
g2 = data.matrix(g2)
g2 = g2[, names(which(apply(g2, 2, function(x){sd(x, na.rm=T)}) > 0.1))]
res = pca(g2)
save(res, file="Pop_pca_NIPTadd1kgCHB_JPT_com0.1_callRate0.99_conf95_SYNCSA.RData")

a = as.data.frame(res$individuals[,1:10], stringsAsFactors = F)
a$ID = c(sample_ID[!duplicated(sample_ID)], ID_1000G$Sample.name[ID_1000G$Population.code %in% c("CHB","CHS","CDX", "JPT")])
a$id = c(rep("Case", 48734), ID_1000G$Superpopulation.code[match(ID_1000G$Sample.name[ID_1000G$Population.code %in% c("CHB","CHS","CDX", "JPT")], ID_1000G$Sample.name)])
a$id1 = c(rep("Case", 48734), ID_1000G$Population.code[match(ID_1000G$Sample.name[ID_1000G$Population.code %in% c("CHB","CHS","CDX", "JPT")], ID_1000G$Sample.name)])
a$id2 = a$id1
a$id2[grep("CH",a$id2)] = "HAN"

write.table(a, file="Pop_pca_NIPTadd1kgCHB_JPT_com0.1_callRate0.99_conf95_SYNCSA.txt", r=F,c=T,sep="\t",quote=F)
```

## 4.4 predicted the ancestry using SVM

```R
a <- read.table("Pop_pca_NIPTadd1kgCHB_JPT_com0.1_callRate0.99_conf95_SYNCSA.txt", header=T,sep="\t",stringsAsFactors=F)

library(e1071) 
train_data <- a[a$id2!="Case", c(1:2)]
test_data <- a[a$id2=="Case", 1:2]

train_labels <- factor(a$id2[a$id2!="Case"], labels=1:3, levels = c("HAN","CDX","JPT"))

svm_model <- svm(train_data, train_labels)  

predicted_labels <- predict(svm_model, test_data)

table(predicted_labels)
```
