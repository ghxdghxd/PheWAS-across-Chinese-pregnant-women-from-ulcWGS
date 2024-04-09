# 5 Genotypes imputation

## 5.1 eagle + minimac4

```sh
for chr in {1..22};
do
{
qctool_v2.0.7 \
    -filetype gen \
    -g $chr.GT.MAF0.01.ICDsample.gz \
    -s all_samples_ICDint.txt \
    -og - | \
    awk '{OFS="\t";
        if($0!~/#/){
            $1="'$chr'";$9="GT";split($3,a,";");$3=a[1];
            gsub("0,0,0","./.",$0);gsub("1,0,0","0/0",$0);gsub("0,1,0","0/1",$0);gsub("0,0,1","1/1",$0);
            print $0}
        else{if(NR==2){print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\""}else{print $0}}}'| \
    bgzip -@5 > $chr.MAF0.01.ICDsample.vcf.gz
} &
done

function run_eagle_chr(){
    chr=$1
    cpus=$2
    Eagle/eagle \
    --numThreads $cpus \
    --chrom $chr \
    --vcfRef 1000G/Eagle_Minimac4_GRCh38_positions_Reference_panels/ALL.chr$chr\_GRCh38.genotypes.20170504.norm.bcf \
    --vcfTarget $chr.MAF0.01.ICDsample.vcf.gz \
    --geneticMapFile Eagle/GRCh38_map/eagle_chr$chr\_b38.map \
    --noImpMissing \
    --allowRefAltSwap \
    --vcfOutFormat z \
    --outputUnphased \
    --outPrefix $chr.MAF0.01.ICDsample.phased
}

function run_minimac4_chr(){
    chr=$1
    cpus=$2
    end=`grep -w ^$chr cytoband_hg38.txt |tail -1|cut -f 3`
    minimac4 \
    --cpus $cpus \
    --chr $chr \
    --start 1 \
    --end $end \
    --minRatio 0.000001 \
    --window 500000 \
    --refhaps 1000G/Eagle_Minimac4_GRCh38_positions_Reference_panels/ALL.chr$chr\_GRCh38.genotypes.20170504.norm.m3vcf.gz \
    --haps $chr.MAF0.01.ICDsample.phased.vcf.gz \
    --noPhoneHome \
    --allTypedSites \
    --format GT,DS,GP \
    --prefix $chr.MAF0.01.ICDsample.imputed \
    --vcfBuffer 10392
}
```

## 5.2 filter imputed genotypes

```sh
for i in *.imputed.dose.vcf.gz;
do
{
out=`echo ${i}|cut -d '.' -f 1-5`
bcftools query -i 'R2>0.5' -f '%ID\t%R2\t[%GT\t]\n' $i |sed 's/0|0/0/g;s/1|1/2/g;s/0|1/1/g;s/1|0/1/g'|bgzip -@10 > $out.R2_0.5.geno.gz
} &
done

for i in *.R2_0.5.geno.gz;
do
{
zcat $i|awk '{if($2>=0.8){print > "'split_R2_0.5/${i%_*}_0.8_1.geno'"}
            else if($2>=0.7 && $2<0.8){print > "'split_R2_0.5/${i%_*}_0.7_0.8.geno'"}
            else if($2>=0.6 && $2<0.7){print > "'split_R2_0.5/${i%_*}_0.6_0.7.geno'"}
            else if($2>0.5 && $2<0.6){print > "'split_R2_0.5/${i%_*}_0.5_0.6.geno'"}}'
} &
done

for j in split_R2_0.5/*.geno;do bgzip -@5 $j;done
```

```R
files = list.files(".", pattern="*MAF0.01.*info", full.names = T)

info = lapply(files, function(x){
    info = fread(x, header=T,sep="\t",stringsAsFactors = F)
    return(info)
})
info = do.call(rbind, info)
info$Rsq = as.numeric(info$Rsq)

write.table(info[info$Rsq>0.5 & (info$MAF>0.01 & info$MAF<0.99),], file="snp_imputed_R0.5_MAF0.01.txt", r=F,c=T,sep="\t",quote=F)
```

## 5.3 SNPs' LDproxy based on CHX of 1000G(hg38)

```sh
################ plink2 LDproxy
for i in {1..22} X;
do
vcftools --gzvcf ../ALL.chr${i}_GRCh38.genotypes.20170504.norm.vcf.gz --plink-tped --out CHB_CHS.chr${i}_GRCH38
zgrep -v "^#" ../ALL.chr${i}_GRCh38.genotypes.20170504.norm.vcf.gz | \
    awk '{print $1":"$2":"$4":"$5}' | paste -d "\t" - CHB_CHS.chr${i}_GRCH38.tped | cut -f 1-2,4- | awk '{OFS="\t";a=$1;$1=$2;$2=a;print $0}' > CHB_CHS.chr${i}_GRCH38.tped.1
mv CHB_CHS.chr${i}_GRCH38.tped.1 CHB_CHS.chr${i}_GRCH38.tped
plink --tfile CHB_CHS.chr${i}_GRCH38 --make-bed --out CHB_CHS.chr${i}_GRCH38
done

for i in {22..1};
do
grep -P "^$i:" imputed_SNP.txt > $i.txt
plink --bfile CHB_CHS.chr${i}_GRCH38 --ld-snp-list $i.txt --ld-window-kb 500000 --ld-window 99999 --r2 --ld-window-r2 0.2 --out ${i}_ld_results --threads 40

awk '{OFS="\t";print $3,$6,$7}' ${i}_ld_results.ld |bgzip -@30 > ${i}_ld_results.txt.gz
rm ${i}.txt
rm ${i}_ld_results.[ln]*
done
```

```R
library(data.table)
library(dplyr)
library(stringr)

files = list.files(".", pattern="*.txt.gz")

filter_ldproxy = function(f){
    ldproxy <- fread(f, header=T, sep="\t", stringsAsFactors=F)
    index <- which(ldproxy$SNP_A == ldproxy$SNP_B)
    ldproxy_self <- ldproxy[index]
    ldproxy <- ldproxy[-index]

    index <- which(ldproxy$SNP_A > ldproxy$SNP_B)

    ldproxy_T <- ldproxy[index] # nolint
    ldproxy_F <- ldproxy[-index] # nolint

    setkey(ldproxy_T, SNP_A, SNP_B)
    colnames(ldproxy_F) <- c("SNP_B", "SNP_A", "R2")
    setkey(ldproxy_F, SNP_A, SNP_B)
    ldproxy_F <- ldproxy_F[!ldproxy_T] # nolint
    colnames(ldproxy_F) <- c("SNP_A", "SNP_B", "R2")

    ldproxy <- rbind(ldproxy_self, ldproxy_T, ldproxy_F)
    save(ldproxy, file=gsub("txt.gz", "RData", f))
}

files = list.files(".", pattern="*_ld_results.RData")

ldproxy = do.call(rbind, mclapply(files, function(x){
    load(x)
    ldproxy$SNP_B1 = gsub(":[atcgATCG]+","", ldproxy$SNP_B)
    return(ldproxy)
}, mc.cores = 22))

snp <- read.table("snp.txt", header = F, sep = "\t", stringsAsFactors = F)
# snp <- snp[grep("^22:", snp$V1), ]
a <- data.table(SNP_A = setdiff(snp$V1, ldproxy$SNP_A))
a$SNP_B <- a$SNP_A # nolint
a$R2 <- 1 # nolint
ldproxy <- rbind(ldproxy, a)

save(ldproxy, file="ld_results.RData")
```

## 5.4 SNPs' LDproxy in ulcWGS

```sh
################################ LD-R2 of imputed variants (MAF>0.01 and imputation score R2>0.5)
for chr in {1..22};
do
{
    bcftools view -i 'R2>0.5' $chr.MAF0.01.ICDsample.imputed.dose.vcf.gz -Oz -o $chr.MAF0.01.ICDsample.imputed.R2_0.5.vcf.gz
    vcftools --gzvcf $chr.MAF0.01.ICDsample.imputed.R2_0.5.vcf.gz --plink-tped --out $chr.MAF0.01.ICDsample.imputed.R2_0.5
    plink --tfile $chr.MAF0.01.ICDsample.imputed.R2_0.5 --make-bed --out $chr.MAF0.01.ICDsample.imputed.R2_0.5
    plink --bfile $chr.MAF0.01.ICDsample.imputed.R2_0.5 --ld-window-kb 500000 --ld-window 99999 --r2 --ld-window-r2 0.2 \
    --out $chr.MAF0.01.ICDsample.imputed.R2_0.5.ld_results --threads 2
} &
done


for i in *.ld;
do
out=`echo $i|sed 's/ld$/txt.gz/g'`
awk '{OFS="\t";print $3,$6,$7}' ${i} |bgzip -@30 > ${out}
# rm ${i}.txt
# rm ${i}_ld_results.[ln]*
done
```

```R
# imputed R2>0.5
files = list.files(".", pattern="*.ld_results.txt.gz$", full.names = T)

minimac4_ld = lapply(files, function(x){
    a <- fread(x, header=T, sep="\t", stringsAsFactors=F)
    return(a)
})
minimac4_ld = do.call(rbind, minimac4_ld)

ldproxy = fread("ld_results.txt.gz", header=T,sep="\t",stringsAsFactors=F)
colnames(ldproxy)[3] = "R2_1kg"

dat = rbind(merge(minimac4_ld, ldproxy, by = c("SNP_A", "SNP_B")), 
    merge(minimac4_ld, ldproxy, by.x = c("SNP_A", "SNP_B"), by.y = c("SNP_B", "SNP_A")))
dat$R2_diff = dat$R2_1kg - dat$R2

dat_filter = dat[abs(dat$R2_diff) < 0.01, ]
write.table(snp_filter, file="snp_filter.txt", r=F,c=F,sep="\t",quote=F)
