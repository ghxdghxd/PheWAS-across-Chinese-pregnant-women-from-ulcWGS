# filter GWAS_catalog

```R
library(data.table)
library(stringr)
gwas <- fread("/share/data3/NIPT/bam_vcf/gwas_catalog_v1.0.2-associations_e105_r2022-04-07.tsv", header=T,sep="\t",stringsAsFactors=F,quote="") # GWAS为hg38
gwas = unique(gwas)
gwas_mat = gwas[, c("STRONGEST SNP-RISK ALLELE","SNPS","CHR_ID","CHR_POS","DISEASE/TRAIT", "MAPPED_GENE", "INITIAL SAMPLE SIZE")]
colnames(gwas_mat) = c("SNPS_ALLELE","SNPS","CHR_ID","CHR_POS","DISEASE_TRAIT", "MAPPED_GENE", "INITIAL_SAMPLE_SIZE")
gwas_mat = unique(gwas_mat)

colnames(gwas) = gsub(" ", "_", colnames(gwas))

a = data.frame(SNPS=grep("[cC]hr:*", unique(gwas_mat$SNPS), v=T), stringsAsFactors=F)
a$snp = gsub("_", ":", a$SNPS)
a$snp = gsub("\\.", ":", a$snp)
a$snp = gsub(" ", "", a$snp)
a$snp = gsub("-", ":", a$snp)
a$snp = gsub("[cC]hr:*", "", a$snp)
a = a[grep("hg18", a$snp, invert=T), ]
a$chr = str_split(a$snp, pattern=":", simplify=T)[,1]
a$pos = str_split(a$snp, pattern=":", simplify=T)[,2]
a$chr[a$chr=="23"] = "X"
a_chr = a$chr
names(a_chr) = a$SNPS
a_pos = a$pos
names(a_pos) = a$SNPS

gwas_mat$CHR_ID[which(gwas_mat$SNPS %in% a$SNPS)] = a_chr[gwas_mat$SNPS[which(gwas_mat$SNPS %in% a$SNPS)]]
gwas_mat$CHR_POS[which(gwas_mat$SNPS %in% a$SNPS)] = a_pos[gwas_mat$SNPS[which(gwas_mat$SNPS %in% a$SNPS)]]
gwas_mat = unique(gwas_mat)

gwas_mat = rbind(gwas_mat[grep(";", gwas_mat$CHR_ID, invert=T),],
    do.call(rbind, lapply(grep(";", gwas_mat$CHR_ID), function(x){
        res = cbind(t(str_split(t(gsub(" ","",gwas_mat[x,1:4])), pattern=";", simplify=T)), gwas_mat[x,5:7])
        colnames(res) = colnames(gwas_mat)
        return(res)
})))

gwas_mat = rbind(gwas_mat[grep("x", gwas_mat$CHR_ID, invert=T),],
    do.call(rbind, lapply(grep("x", gwas_mat$CHR_ID), function(x){
        res = cbind(t(str_split(t(gsub(" ","",gwas_mat[x,1:4])), pattern="x", simplify=T)), gwas_mat[x,5:7])
        colnames(res) = colnames(gwas_mat)
        return(res)
})))
gwas_mat = unique(gwas_mat)

gwas_mat = rbind(gwas_mat[grep(" x ", gwas_mat$SNPS, invert=T),],
    do.call(rbind, lapply(grep(" x ", gwas_mat$SNPS), function(x){
        res = cbind(t(str_split(gsub(" ","",gwas_mat[x,1:2]), pattern="x", simplify=T)), gwas_mat[x,3:7])
        colnames(res) = colnames(gwas_mat)
        return(res)
})))
gwas_mat = unique(gwas_mat)

gwas_mat = rbind(gwas_mat[grep("\\|", gwas_mat$SNPS, invert=T),],
    do.call(rbind, lapply(grep("\\|", gwas_mat$SNPS), function(x){
        res = cbind(t(str_split(gsub(" ","",gwas_mat[x,1:2]), pattern="\\|", simplify=T)), gwas_mat[x,3:7])
        colnames(res) = colnames(gwas_mat)
        return(res)
})))
gwas_mat = unique(gwas_mat)

gwas_mat = rbind(gwas_mat[grep("; ", gwas_mat$SNPS, invert=T),],
    do.call(rbind, lapply(grep("; ", gwas_mat$SNPS), function(x){
        res = cbind(t(str_split(gsub(" ","",gwas_mat[x,1:2]), pattern="; ", simplify=T)), gwas_mat[x,3:7])
        colnames(res) = colnames(gwas_mat)
        return(res)
})))
gwas_mat = unique(gwas_mat)

gwas_mat = rbind(gwas_mat[grep(";", gwas_mat$SNPS, invert=T),],
    do.call(rbind, lapply(grep(";", gwas_mat$SNPS), function(x){
        res = cbind(t(str_split(gsub(" ","",gwas_mat[x,1:2]), pattern=";", simplify=T)), gwas_mat[x,3:7])
        colnames(res) = colnames(gwas_mat)
        return(res)
})))
gwas_mat = unique(gwas_mat)

gwas_mat = rbind(gwas_mat[grep(",", gwas_mat$SNPS, invert=T),],
    do.call(rbind, lapply(grep(",", gwas_mat$SNPS), function(x){
        res = cbind(t(str_split(gsub(" ","",gwas_mat[x,1:2]), pattern=",", simplify=T)), gwas_mat[x,3:7])
        colnames(res) = colnames(gwas_mat)
        return(res)
})))
gwas_mat = unique(gwas_mat)

gwas_mat = gwas_mat[grep("esv", gwas_mat$SNPS, invert=T),]
gwas_mat = gwas_mat[grep("kgp", gwas_mat$SNPS, invert=T),]
gwas_mat = gwas_mat[grep("^Affx", gwas_mat$SNPS, invert=T),]
gwas_mat = gwas_mat[grep("\\*", gwas_mat$SNPS, invert=T),]

gwas_mat$SNPS = gsub("^psy_", "", gwas_mat$SNPS)
gwas_mat$SNPS = gsub("^ ", "", gwas_mat$SNPS)
gwas_mat$SNPS = gsub("^imm_", "", gwas_mat$SNPS)
gwas_mat$SNPS_ALLELE = gsub("^psy_", "", gwas_mat$SNPS_ALLELE)
gwas_mat$SNPS_ALLELE = gsub("^ ", "", gwas_mat$SNPS_ALLELE)
gwas_mat$SNPS_ALLELE = gsub("^imm_", "", gwas_mat$SNPS_ALLELE)
gwas_mat = gwas_mat[grep("^im:", gwas_mat$SNPS, invert=T),]

a = data.frame(SNPS=gwas_mat$SNPS[which(gwas_mat$CHR_ID=="" & !grepl("rs", gwas_mat$SNPS) & grepl(":", gwas_mat$SNPS))], stringsAsFactors=F)
a$snp = gsub("ch", "", a$SNPS)
a$snp = gsub("\\]r", ":", a$snp)
a$snp = gsub("^del-", "", a$snp)
a$snp = gsub("^e", "", a$snp)
a$snp = gsub("^:", "", a$snp)
a$snp = gsub("_", ":", a$snp)
a = a[grep("hg18", a$snp, invert=T), ]
a$chr = str_split(a$snp, pattern=":", simplify=T)[,1]
a$pos = str_split(a$snp, pattern=":", simplify=T)[,2]
a_chr = a$chr
names(a_chr) = a$SNPS
a_pos = a$pos
names(a_pos) = a$SNPS

gwas_mat$CHR_ID[gwas_mat$SNPS %in% a$SNPS] = a_chr[gwas_mat$SNPS[gwas_mat$SNPS %in% a$SNPS]]
gwas_mat$CHR_POS[gwas_mat$SNPS %in% a$SNPS] = a_pos[gwas_mat$SNPS[gwas_mat$SNPS %in% a$SNPS]]
gwas_mat = unique(gwas_mat)
gwas_mat$CHR_POS = as.numeric(gwas_mat$CHR_POS)

gwas_mat$SNPS = gsub(":IND", "", gwas_mat$SNPS)
gwas_mat$SNPS = gsub("[ †:a-z]$", "", gwas_mat$SNPS)

gwas_mat$SNPS[which(gwas_mat$CHR_ID=="" & grepl("rs", gwas_mat$SNPS))] = str_split(gwas_mat$SNPS[which(gwas_mat$CHR_ID=="" & grepl("rs", gwas_mat$SNPS))], pattern="[_-]", simplify=T)[,1]
gwas_mat$SNPS[gwas_mat$SNPS=="s2736898"] = "rs2736898"
gwas_mat$SNPS[gwas_mat$SNPS=="r12465214"] = "rs12465214"
gwas_mat$SNPS[gwas_mat$SNPS=="ss1388061783"] = "rs1388061783"

gwas_pop = data.frame(INITIAL_SAMPLE_SIZE = unique(gwas_mat$INITIAL_SAMPLE_SIZE), stringsAsFactors=F)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("[0-9]+,[0-9]+ ", "", gwas_pop$INITIAL_SAMPLE_SIZE)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("[0-9]+ ", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("individuals *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("ancestry *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("cases *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("controls *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("female *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("male *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("women *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub("men *", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)
gwas_pop$INITIAL_SAMPLE_SIZE_1 = gsub(" +$", "", gwas_pop$INITIAL_SAMPLE_SIZE_1)


gwas_pop_1 = data.frame(INITIAL_SAMPLE_SIZE = unique(gwas_pop$INITIAL_SAMPLE_SIZE_1), stringsAsFactors=F)

# AFR, African
# AMR, Ad Mixed American
# EAS, East Asian
# EUR, European
# SAS, South Asian

index = grepl("European", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("American", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("African", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("British", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Peruvian", gwas_pop_1$INITIAL_SAMPLE_SIZE)
gwas_pop_1$OTH = ifelse(index, "OTH", "")
gwas_pop_1$AS = ifelse(grepl("Asian", gwas_pop_1$INITIAL_SAMPLE_SIZE), "AS", "")

index = grepl("East Asian", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Japanese", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Thai", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Korean", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Vietnamese", gwas_pop_1$INITIAL_SAMPLE_SIZE)
gwas_pop_1$EAS = ifelse(index, "EAS", "")
gwas_pop_1$CHB = ifelse(grepl("Chinese", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Taiwanese", gwas_pop_1$INITIAL_SAMPLE_SIZE), "CHB", "")
gwas_pop_1$SAS = ifelse(grepl("South Asian", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Indian", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Bangladeshi", gwas_pop_1$INITIAL_SAMPLE_SIZE) | grepl("Indonesian", gwas_pop_1$INITIAL_SAMPLE_SIZE), "SAS", "")
gwas_pop_1$OTH[which(gwas_pop_1$OTH == "" & gwas_pop_1$AS=="" & gwas_pop_1$EAS=="" & gwas_pop_1$CHB=="" & gwas_pop_1$SAS=="")] = "OTH"

gwas_pop = cbind(gwas_pop, gwas_pop_1[match(gwas_pop$INITIAL_SAMPLE_SIZE_1, gwas_pop_1$INITIAL_SAMPLE_SIZE),-1])

gwas_mat = cbind(gwas_mat, gwas_pop[match(gwas_mat$INITIAL_SAMPLE_SIZE, gwas_pop$INITIAL_SAMPLE_SIZE),-c(1:2)])

bed <- fread("/share/data3/apps/annovar/humandb/hg38_avsnp150.txt", header=F, sep="\t", stringsAsFactors=F)
colnames(bed) = c("chr","start","end","ref","alt","rsid")
bed = unique(bed[, c(1,3,6)])
save(bed, file="hg38_avsnp150.RData")

gwas_bed = bed[match(gwas_mat$SNPS, bed$rsid), ]
gwas_mat[which(gwas_mat$CHR_ID!="" & gwas_mat$CHR_ID != gwas_bed$chr), c("CHR_ID", "CHR_POS")] = gwas_bed[which(gwas_mat$CHR_ID!="" & gwas_mat$CHR_ID != gwas_bed$chr), c("chr","end")]
gwas_mat[which(gwas_mat$CHR_ID==""), c("CHR_ID", "CHR_POS")] = gwas_bed[which(gwas_mat$CHR_ID==""), c("chr","end")]

gwas_mat$EAS[which(gwas_mat$CHB == "CHB")] = "EAS"
gwas_mat$AS[which(gwas_mat$CHB == "CHB" | gwas_mat$EAS == "EAS" | gwas_mat$SAS == "SAS")] = "AS"
save(gwas_mat, file='gwas_catalog_mat_1.RData')


```
