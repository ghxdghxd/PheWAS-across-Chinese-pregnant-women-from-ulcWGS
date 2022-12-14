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
