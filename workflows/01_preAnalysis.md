# 1. preAnalysis

## 1.1 alignment

```sh
bwa aln -t 35 hg38.fa.gz $NIPT_single_end.fq | bwa samse -r "@RG\tID:$sampleID\tLB:HG38\tPL:XXX\tPU:XXX\tSM:$sampleID\tCREATE_INDEX:True" hg38.fa.gz - $NIPT_single_end.fq | \
    samtools view -Su -@4 - | samtools sort -@4 - -o $sampleID.sorted.bam
```

## 1.2 remove duplicated reads

```sh
/data/NIPT/script/MarkDuplicates.py \
    --path2GATK gatk-4.1.4.1/gatk \
    -i $sampleID.sorted.bam \
    -o $sampleID -t 10 -m 10g \
    --validation_stringency LENIENT \
    --remove_duplicates
```

## 1.3 Base Quality Score Recalibration

```sh
/data/NIPT/script/BQSR.py \
    --path2GATK gatk-4.1.4.1/gatk \
    -r hg38.fa.gz \
    -i $sampleID.sorted.dedup.bam \
    -o $sampleID \
    -l GATK_reaource/hg38_v0/wgs_calling_regions.hg38.interval_list \
    --known-sites GATK_reaource/hg38_v0/dbsnp_146.hg38.vcf.gz \
    --known-sites GATK_reaource/hg38_v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites GATK_reaource/hg38_v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
```
