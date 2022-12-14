# 6. CNV and Fetal fraction

## 6.1 CNV

### 6.1.1 Convert bam to npz

```sh
WisecondorX convert $ampleID.bqsr.bam $sampleID.npz
```

### 6.1.2 Create reference

```sh
WisecondorX newref samples_all_npz_ref/*.npz \
    reference.100kb.npz --nipt \
    --binsize 100000 --cpus 38
```

### 6.1.3 Predict copy number alterations

```sh
WisecondorX predict $sampleID.npz reference.100kb.npz $sampleID/$sampleID --bed --plot
```

## 6.2 Fetal fraction

### 6.2.1 FFY

```sh
bedtools intersect -a $sampleID.vcf.gz -b wgEncodeDukeMapabilityUniqueness35bp.one.rm_um35.hg38.bed | \
    cut -f 1|uniq -c|awk 'BEGIN{sum=0}{sum+=$1;if($2=="chrY"){y=$1};}END{print "'$sampleID'",y,sum}' >> samples_Y_readCount_var_35bp.txt
```

```R
chrYAM = 0.00173 #the median of the Y chromosome percentage from three adult males,https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4956272/
# chrYAM_Optimized = 0.000374
Ycount <- read.table("FFY/samples_Y_readCount_var_35bp.txt", header=F,sep=" ", stringsAsFactors=F)
colnames(Ycount) = c("ID","chry","chrall")
Ycount$chryP_Original = Ycount$chry/Ycount$chrall
chryF_Original = mean(Ycount$chryP_Original[which(Ycount$gender=="F")])

Ycount$FFY = (Ycount$chryP_Original - chryF_Original)/(chrYAM - chryF_Original)
Ycount$FFY[-which(Ycount$gender=='M')] = 0
Ycount$FFY[which(Ycount$FFY < 0)] = 0

write.table(Ycount, file="FFY/samples_Y_readCount_var_35bp_conf.txt", r=F,c=T,quote=F)
```

### 6.2.2 PREFACE Model training and Predicting

```sh
Rscript PREFACE-master/PREFACE.R train \
    --config samples_Y_readCount_var_35bp_conf.txt \
    --outdir PREFACE_model_training \
    --cpus 39

Rscript ./PREFACE-master/PREFACE.R predict \
    --infile $sampleID_bins.bed \
    --model PREFACE_model_training/model.RData \
    --json fetal_fraction/$sampleID.json
```

```R
library(rjson)
library(parallel)
files = list.files("fetal_fraction", pattern="*.json",full.names=T)
mat = mclapply(files, function(x){
    a = as.data.frame(fromJSON(file = x))
    a$ID = gsub(".json", "", basename(x))
    return(a)
}, mc.cores = 48)
mat = do.call(rbind, mat)

Ycount = read.table("samples_Y_readCount_var_35bp_conf.txt", header=T,sep=" ",stringsAsFactors=F)
mat = merge(mat, Ycount, by = "ID")

write.table(mat, file="samples_Y_readCount_var_35bp_conf_addPREFACE.txt", r=F,c=T,sep="\t",quote=F)
```
