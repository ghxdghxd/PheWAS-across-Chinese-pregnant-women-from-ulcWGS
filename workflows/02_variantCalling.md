# 2 Variants calling

## 2.1 calling variants of each sample

```sh
bcftools mpileup -Ou --adjust-MQ 50 --min-BQ 20 --min-MQ 30 --per-sample-mF -a "DP,AD" -f hg38.fa.gz \
    $sampleID.sorted.dedup.bqsr.bam | \
    bcftools call --multiallelic-caller --keep-alts -Ou| \
    bcftools filter -s LowQual -e 'MQ="." || MQ<30 || DP<2' --set-GTs . \
    -Oz > $sampleID.sorted.dedup.bqsr.b.vcf.gz > bcftools.nohup.out 2>&1 &
```

## 2.2 normalized variants of each sample

```sh
bcftools view -e 'ALT =="."'$sampleID.sorted.dedup.bqsr.b.vcf.gz -O z -o $sampleID.sorted.dedup.bqsr.b.var.vcf.gz
bcftools norm -m -both -f hg38.fa.gz $sampleID.sorted.dedup.bqsr.b.var.vcf.gz -N -Oz -o $sampleID.sorted.dedup.bqsr.b.var.norm.vcf.gz
```

## 2.3 all unique varinats

```sh
bcftools query -i 'ALT!="."' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t.\t.\t.\n' $sampleID.sorted.dedup.bqsr.b.var.norm.vcf.gz >> $sampleID.lite.vcf
awk '!visited[$0]++' $sampleID.lite.vcf > $sampleID.lite.vcf.1
mv $sampleID.lite.vcf.1 $sampleID.lite.vcf
```

# 2.4 annotated unique variants

```sh
# combine
cat *.lite.vcf| awk '!visited[$0]++' > all.lite.vcf
# annovar
/usr/bin/time -av -o chr$i.lite.vcf.annovar.log \
perl annovar/table_annovar.pl \
    all.lite.vcf \
    -vcfinput annovar/humandb \
    -buildver hg38 \
    -out all.lite \
    -otherinfo -remove \
    -protocol refGene,avsnp150,gnomad_genome \
    -operation g,f,f -nastring .
cut -f 1-19,23-24,26-27 all.lite.hg38_multianno.txt | bgzip -@2 > all.lite.hg38_multianno.txt.gz
rm all.lite.hg38_multianno.txt
```

## 2.5 intersected with ChinaMAP

```sh
bcftools isec all.lite.vcf ChinaMAP/mbiobank_ChinaMAP.phase1.vcf.gz -O z -p .
```

## 2.6 variants selecting

```R
nipt <- fread("0000.vcf.gz", header=T, sep="\t", stringsAsFactors=F)
colnames(nipt)[1] = "CHROM"
nipt = nipt[,1:6]
colnames(nipt)[6] = "AF_ChinaMAP"

nipt_chinamap <- fread("0003.vcf.gz", header=T, sep="\t", stringsAsFactors=F)
colnames(nipt_chinamap)[1] = "CHROM"
library(stringr)
nipt_chinamap = cbind(nipt_chinamap, as.data.frame(str_split(nipt_chinamap$INFO, pattern=";|=", simplify=T)[,c(2,4,6)], stringsAsFactors = F))
colnames(nipt_chinamap)[9:11] = c("AC", "AF", "AN")
nipt_chinamap = nipt_chinamap[, c(1:5, 10)]
colnames(nipt_chinamap)[6] = "AF_ChinaMAP"

nipt = rbind(nipt, nipt_chinamap)
setkey(nipt, CHROM, POS, REF, ALT)

anno = fread("all.lite.hg38_multianno.txt.gz", header=T,sep="\t",stringsAsFactors=F)
setkey(anno, Otherinfo4, Otherinfo5, Otherinfo7, Otherinfo8)
anno = anno[nipt, nomatch=0]
fwrite(anno, file="all.lite.hg38_multianno.ChinaMAP.txt.gz", row.names = F, col.names = T, sep = "\t", quote=F)
anno[, c("GeneDetail.refGene", "gnomAD_genome_ALL", "gnomAD_genome_AFR", "gnomAD_genome_AMR", "gnomAD_genome_ASJ", "gnomAD_genome_FIN", "gnomAD_genome_NFE", "gnomAD_genome_OTH")] = NULL
germline = which(anno$gnomAD_genome_EAS != "." | anno$AF_ChinaMAP != "." | grepl("rs",anno$avsnp150) | grepl("rs",anno$ID))
anno = anno[germline, c(1:5,12:15,11,17)]
fwrite(anno, file="all.avsnp150.ChinaMAP.bed.gz", row.names = F, col.names = F, sep = "\t", quote=F)
```

## 2.7 variants filtering

We filtered the variants:
    *located in the 35-kmer problematic alignmant region and the regionis with a mappability uniqueness score of less than one
    * AF == 0

```sh
zcat all.avsnp150.ChinaMAP.bed.gz|cut -f 6-9 | awk -F '\t' '$10>0||$11>0{OFS="\t";print $1,$2,".",$3,$4,".",".","."}' > all.avsnp150.ChinaMAP.AF0.lite.vcf
intersectBed -a all.avsnp150.ChinaMAP.AF0.lite.vcf -b wgEncodeDukeMapabilityUniqueness35bp.one.rm_um35.hg38.bed > all.avsnp150.ChinaMAP.AF0.rm_um35.lite.vcf
intersectBed -wb -a all.avsnp150.ChinaMAP.AF0.rm_um35.lite.vcf -b cytoband_hg38_chr.txt|cut -f 1-2,4-5,12|bgzip -@10 > all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.lite.bed.gz
```

```R
library(data.table)
bed = fread("all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.lite.bed.gz", header=F,sep="\t",stringsAsFactors=F)
setkey(bed, V1, V2, V3, V4)

all_bed = fread("all.avsnp150.ChinaMAP.bed.gz", header= F, sep = "\t", stringsAsFactors = F)
setkey(all_bed, V6,V7,V8,V9)
all_bed = bed[b, nomatch=0]
all_bed = all_bed[,c(1:5, 11:12)]
colnames(all_bed) = c("chr", "pos", "ref", "alt", "cytoband", "AF_EAS", "AF_ChinaMAP")
all_bed$type = ""
all_bed$type[all_bed$AF_EAS == "." & all_bed$AF_ChinaMAP != "."] = "C"  # ChinaMAP only
all_bed$type[all_bed$AF_EAS != "." & all_bed$AF_ChinaMAP != "."] = "B"  # Both ChinaMAP and EAS
all_bed$type[all_bed$AF_EAS != "." & all_bed$AF_ChinaMAP == "."] = "E"  # EAS only
all_bed$AF = all_bed$AF_ChinaMAP
all_bed$AF[all_bed$type == "E"] = all_bed$AF_EAS[all_bed$type == "E"]
all_bed = unique(all_bed[, c(1:5,8:9)])
fwrite(all_bed, file = "all.avsnp150.ChinaMAP.AF0.rm_um35.cytoband.AF.lite.bed.gz", row.names = F, col.names = T, sep = "\t", quote = F)
```
