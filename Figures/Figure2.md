# Figure2

## Figure2A

```R
load("icd_status_PC_all.RData")
clin_status = clin_status[,names(which(colSums(clin_status[,-c(1:7)], na.rm = T)>30))]

load("gwas_catalog_mat.RData")
gwas_mat$id = paste(gwas_mat$CHR_ID, gwas_mat$CHR_POS, sep=":")
gwas_mat = as.data.frame(gwas_mat, stringsAsFactors=F)
gwas_trait = read.csv("gwas_catalog_mat_trait.txt", header=F,sep="\t",stringsAsFactors=F)
gwas_mat$DISEASE_TRAIT_zh = gwas_trait$V2[match(gwas_mat$`DISEASE/TRAIT`, gwas_trait$V1)]
colnames(gwas_mat)[4] = "DISEASE_TRAIT"
gwas_catalog_trait_ICD10 = read.table("gwas_catalog_mat_trait_ICD10.tsv",header=T,sep="\t",stringsAsFactors=F)
gwas_mat$ICD10 = gwas_catalog_trait_ICD10$ICD10[match(gwas_mat$DISEASE_TRAIT, gwas_catalog_trait_ICD10$Disease.trait)]

library(data.table)
library(stringr)
imputed_info = fread("snp_imputed_R0.5_MAF0.01_anno.txt",header=T,sep="\t",stringsAsFactor=F)
imputed_info$EAS_AF = as.numeric(imputed_info$EAS_AF)
imputed_info$AF_sd = sqrt((imputed_info$EAS_AF * (1-imputed_info$EAS_AF))/25639)

imputed_info = cbind(imputed_info, str_split(imputed_info$SNP, pattern=":", simplify=T))

gwas_mat_sub = gwas_mat[which(gwas_mat$id %in% intersect(paste(imputed_info$V1, imputed_info$V2, sep=":"), gwas_mat$id)), ]
gwas_mat_sub = gwas_mat_sub[gwas_mat_sub$ICD10 %in% names(clin_status), ]

imputed_info_sub = imputed_info[which(paste(imputed_info$V1, imputed_info$V2, sep = ":") %in% gwas_mat_sub$id), ]
imputed_info_sub$id = paste(imputed_info_sub$V1, imputed_info_sub$V2, sep = ":")

res_gwas = mclapply(unique(gwas_mat_sub$ICD10), function(x){
    res = fread(paste0(x, ".SAIGE.txt.gz"), header=T, sep="\t", stringsAsFactors=F)
    res = res[res$MarkerID %in% imputed_info_sub$SNP[which(imputed_info_sub$id %in% gwas_mat_sub$id[gwas_mat_sub$ICD10==x])], ]
    res$icd = x
    return(res)
}, mc.cores = 20)
res_gwas = do.call(rbind, res_gwas)

res_gwas$snp_id = imputed_info_sub$id[match(res_gwas$MarkerID, imputed_info_sub$SNP)]
res_gwas$OR = exp(res_gwas$BETA)

# observed MAF, res_MAF_M for mother and res_MAF_C for child
load("all_R2_0.5_MAF_in_motherICD_childICD.RData")

library(genpwr)

res_gwas$obMAF = res_MAF_M[res_gwas$MarkerID, "MAF"]
case_num = colSums(clin_status, na.rm=T)
case_control_num = colSums(!is.na(clin_status), na.rm=T)
case_rate = case_num/case_control_num

pow = mclapply(1:nrow(res_gwas), function(x){
    return(genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
            N=case_control_num[res_gwas$icd][x],
            Case.Rate=case_rate[res_gwas$icd][x],
            k=NULL,
            MAF=res_MAF_M[res_gwas$MarkerID,"MAF"][x],
            OR=res_gwas$OR[x],
            Alpha=0.05,
            True.Model=c("Additive"),
            Test.Model=c("Additive")))
}, mc.cores = 10)
pow = do.call(rbind, pow)
res_gwas$power = unlist(pow)

get_fraction2 = function(min_obMAF, min_power){
    SNPs = intersect(rownames(res_MAF_M)[res_MAF_M$MAF > min_obMAF], res_gwas$MarkerID[which(res_gwas$power > min_power)])
    ids = imputed_info_sub$id[imputed_info_sub$SNP %in% SNPs]

    ## number of snp-ICD10 pairs in GWAS
    a = length(unique(paste(gwas_mat_sub$id, gwas_mat_sub$ICD10, sep=":")[gwas_mat_sub$id %in% ids]))
    ## number of snp-ICD10 pairs replicated in study
    b = length(intersect(paste(res_gwas$snp_id, res_gwas$icd, sep=":")[res_gwas$MarkerID %in% SNPs & res_gwas$p.value<0.05],
        paste(gwas_mat_sub$id, gwas_mat_sub$ICD10, sep=":")[gwas_mat_sub$id %in% ids]))
    return(c(min_obMAF = min_obMAF, min_power = min_power, snpICD10_all_GWAS = a, snpICD10_in_GWAS = b))
}

res_gwas_fraction2 = do.call(rbind, lapply(c(0, 0.01, 0.05, 0.1), function(min_obMAF){
    return(do.call(rbind, lapply(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), function(min_pow){
        return(get_fraction2(min_obMAF, min_pow))
    })))
}))
res_gwas_fraction2 = as.data.frame(res_gwas_fraction2, stringsAsFactors=F)
res_gwas_fraction2$min_obMAF = paste("MAF>", res_gwas_fraction2$min_obMAF)
res_gwas_fraction2$snpICD10_all_GWAS = as.numeric(res_gwas_fraction2$snpICD10_all_GWAS)
res_gwas_fraction2$snpICD10_in_GWAS = as.numeric(res_gwas_fraction2$snpICD10_in_GWAS)
res_gwas_fraction2$fraction = signif(res_gwas_fraction2$snpICD10_in_GWAS/res_gwas_fraction2$snpICD10_all_GWAS, 3)*100
res_gwas_fraction2 = res_gwas_fraction2 %>% filter(min_obMAF!="MAF>0")
res_gwas_fraction2$pop = "ALL"

#replication in different populations
res_gwas_pop = lapply(1:nrow(res_gwas), function(x){
    return(gwas_mat[which(gwas_mat$id==res_gwas$snp_id[x] & gwas_mat$ICD10==res_gwas$icd[x]),c("id","ICD10","OTH", "AS","EAS","CHB","SAS")])
})
res_gwas_pop = do.call(rbind, res_gwas_pop)
res_gwas_pop = merge(res_gwas, res_gwas_pop, by.x=c("snp_id","icd"), by.y=c("id","ICD10"))

res_gwas_pop$gene = as.character(res_gwas$gene[match(res_gwas_pop$MarkerID, res_gwas$MarkerID)])

get_fraction2_pop = function(min_obMAF, min_power, pop){
    SNPs = intersect(rownames(res_MAF_M)[res_MAF_M$MAF > min_obMAF], res_gwas_pop$MarkerID[which(res_gwas_pop$power > min_power & res_gwas_pop[, ..pop] == pop)])
    ids = imputed_info_sub$id[imputed_info_sub$SNP %in% SNPs]

    ## number of snp-ICD10 pairs in GWAS
    a = length(unique(paste(gwas_mat_sub$id, gwas_mat_sub$ICD10, sep=":")[gwas_mat_sub$id %in% ids & gwas_mat_sub[, pop] == pop]))
    ## number of snp-ICD10 pairs replicated
    b = length(intersect(paste(res_gwas_pop$snp_id, res_gwas_pop$icd, sep=":")[res_gwas_pop$MarkerID %in% SNPs & res_gwas_pop$p.value<0.05 & res_gwas_pop[, ..pop] == pop],
        paste(gwas_mat_sub$id, gwas_mat_sub$ICD10, sep=":")[gwas_mat_sub$id %in% ids & gwas_mat_sub[, pop] == pop]))
    return(c(min_obMAF = min_obMAF, min_power = min_power, pop = pop, snpICD10_all_GWAS = a, snpICD10_in_GWAS = b))
}
###
res_gwas_fraction2_pop = do.call(rbind, lapply(c(0, 0.01, 0.05, 0.1), function(min_obMAF){
    return(do.call(rbind, lapply(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), function(min_pow){
        return(do.call(rbind, lapply(c("AS","EAS","CHB","SAS","OTH"), function(pop){
            return(get_fraction2_pop(min_obMAF, min_pow, pop))
        })))
    })))
}))

res_gwas_fraction2_pop = as.data.frame(res_gwas_fraction2_pop, stringsAsFactors=F)
res_gwas_fraction2_pop$min_obMAF = paste("MAF>", res_gwas_fraction2_pop$min_obMAF)
res_gwas_fraction2_pop$snpICD10_all_GWAS = as.numeric(res_gwas_fraction2_pop$snpICD10_all_GWAS)
res_gwas_fraction2_pop$snpICD10_in_GWAS = as.numeric(res_gwas_fraction2_pop$snpICD10_in_GWAS)
res_gwas_fraction2_pop$fraction = signif(res_gwas_fraction2_pop$snpICD10_in_GWAS/res_gwas_fraction2_pop$snpICD10_all_GWAS, 3)*100

res_gwas_fraction2_pop = rbind(res_gwas_fraction2, res_gwas_fraction2_pop)

save(res_gwas_pop, res_gwas_fraction2_pop, file="res_gwas_fraction_pop.RData")

##
res_gwas_pop$ICDanno = ""
res_gwas_pop$ICDanno[grep("[AB]", res_gwas_pop$icd)] = "Certain infectious and parasitic diseases"
res_gwas_pop$ICDanno[which(res_gwas_pop$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24", "D25", "D26"))] = "Neoplasms"
res_gwas_pop$ICDanno[which(res_gwas_pop$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
res_gwas_pop$ICDanno[grep("E", res_gwas_pop$icd)] = "Endocrine, nutritional and metabolic diseases"
res_gwas_pop$ICDanno[grep("F", res_gwas_pop$icd)] = "Mental and behavioural disorders"
res_gwas_pop$ICDanno[which(res_gwas_pop$icd %in% c("H00", "H01", "H02", "H04", "H10","H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
res_gwas_pop$ICDanno[which(res_gwas_pop$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92","H93"))] = "Ear and mastoid process"
res_gwas_pop$ICDanno[grep("I", res_gwas_pop$icd)] = "Circulatory system"
res_gwas_pop$ICDanno[grep("J", res_gwas_pop$icd)] = "Respiratory system"
res_gwas_pop$ICDanno[grep("K", res_gwas_pop$icd)] = "Digestive system"
res_gwas_pop$ICDanno[grep("L", res_gwas_pop$icd)] = "Skin and subcutaneous tissue"
res_gwas_pop$ICDanno[grep("M", res_gwas_pop$icd)] = "Musculoskeletal system and connective tissue"
res_gwas_pop$ICDanno[grep("N", res_gwas_pop$icd)] = "Genitourinary system"
res_gwas_pop$ICDanno[grep("O", res_gwas_pop$icd)] = "Pregnancy, childbirth and the puerperium"
res_gwas_pop$ICDanno[grep("P", res_gwas_pop$icd)] = "Certain conditions originating in the perinatal period"
res_gwas_pop$ICDanno[grep("Q", res_gwas_pop$icd)] = "Congenital malformations, deformations and chromosomal abnormalities"
res_gwas_pop$ICDanno = factor(res_gwas_pop$ICDanno, levels = sort(unique(res_gwas_pop$ICDanno))[c(12,8,6,4,13,14,10,2,7,5,11,3,9,1)])
res_gwas_pop$pop = ""
res_gwas_pop$pop[res_gwas_pop$CHB=="CHB"]="CHB"
res_gwas_pop$pop[res_gwas_pop$EAS=="EAS" & res_gwas_pop$pop == ""]="EAS"
res_gwas_pop$pop[res_gwas_pop$SAS=="SAS" & res_gwas_pop$pop == ""]="SAS"
res_gwas_pop$pop[res_gwas_pop$AS=="AS" & res_gwas_pop$pop == ""]="AS"
res_gwas_pop$pop[res_gwas_pop$OTH=="OTH" & res_gwas_pop$pop == ""]="ALL"

colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
    "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
    "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")


pdf("Figure2A.pdf",width=10, height=2.4)
ggarrange(newpage = F,
    ggplot(res_gwas_pop) + geom_point(aes(power, -log10(p.value), color=ICDanno), size=0.6) +
        geom_hline(aes(yintercept=-log10(0.05)), linetype = "dashed") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 1, override.aes = list(size = 2))) +
        scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
        scale_y_continuous(position='right',limits = c(0,10), breaks = seq(0,10,2.5)) +
        theme_bw() +
        theme(aspect.ratio= 0.8, legend.key.size = unit(0, 'lines')),
    res_gwas_fraction2_pop %>% filter(min_obMAF=="MAF>0.01", pop %in% c("CHB", "EAS", "AS", "ALL")) %>%
        mutate(pop = factor(pop, levels = c("CHB", "EAS", "AS", "ALL"))) %>%
        ggplot(aes(as.numeric(min_power), fraction, color = pop, group = pop)) +
            geom_point(size=2, shape = 15) +
            geom_smooth(method="lm", se = FALSE) +
            scale_color_manual(values=pal_d3("category10")(5)[c(2:4,1)]) +
            labs(x="minPower in PheWAS", y="Replicated(%)", fill="") +
            scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
            scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
            theme_bw() +
            theme(legend.position="right",
                legend.title = element_blank(),
                legend.background = element_blank(),
                axis.text.x = element_text(),
                aspect.ratio = 0.8),
    ncol=1)
dev.off()
```

## Figure2B

```R
load("icd_status_PC_all.RData")
clin_status = clin_status[, names(which(colSums(clin_status[,-c(1:7)], na.rm = T) > 30))]

load("gwas_catalog_mat.RData")
gwas_mat$id = paste(gwas_mat$CHR_ID, gwas_mat$CHR_POS, sep=":")
gwas_mat = as.data.frame(gwas_mat, stringsAsFactors=F)
gwas_trait = read.csv("gwas_catalog_mat_trait.txt", header=F,sep="\t",stringsAsFactors=F)
gwas_mat$DISEASE_TRAIT_zh = gwas_trait$V2[match(gwas_mat$`DISEASE/TRAIT`, gwas_trait$V1)]
colnames(gwas_mat)[4] = "DISEASE_TRAIT"
gwas_catalog_trait_ICD10 = read.table("gwas_catalog_trait_ICD10.tsv",header=T,sep="\t",stringsAsFactors=F)
gwas_mat$ICD10 = gwas_catalog_trait_ICD10$ICD10[match(gwas_mat$DISEASE_TRAIT, gwas_catalog_trait_ICD10$Disease.trait)]

library(data.table)
library(stringr)
imputed_info = fread("snp_imputed_R0.5_MAF0.01.txt",header=T,sep="\t",stringsAsFactor=F)
imputed_info = cbind(imputed_info, str_split(imputed_info$SNP, pattern=":", simplify=T))


gwas_mat_sub = gwas_mat[which(gwas_mat$ICD10 %in% names(clin_status)), ]

lod("ld_results.RData")
ldproxy$SNP_B1 = sub("([0-9]+):([0-9]+):([ATCG]*):([ATCG]*)", "\\1:\\2", ldproxy$SNP_B)

ldproxy_sub = ldproxy[which(ldproxy$SNP_A %in% imputed_info$SNP & ldproxy$SNP_B1 %in% gwas_mat_sub$id), ]
gwas_mat_sub = gwas_mat_sub[which(gwas_mat_sub$id %in% ldproxy_sub$SNP_B1), ]

res_gwas_ld = mclapply(unique(gwas_mat_sub$ICD10), function(x){
    res = fread(paste0(x, ".SAIGE.txt.gz"), header=T, sep="\t", stringsAsFactors=F)
    res = res[res$MarkerID %in% ldproxy_sub$SNP_A[which(ldproxy_sub$SNP_B1 %in% gwas_mat_sub$id[gwas_mat_sub$ICD10==x])], ]
    res = merge(res, ldproxy_sub, by.x="MarkerID",by.y="SNP_A")
    res = res[which(res$SNP_B1 %in% gwas_mat_sub$id[gwas_mat_sub$ICD10==x]),]
    res$icd = x
    return(res)
}, mc.cores = 20)
res_gwas_ld = do.call(rbind, res_gwas_ld)

res_gwas_ld$OR = exp(res_gwas_ld$BETA)

save(res_gwas_ld, file="res_gwas_ld_fraction.RData")

res_gwas_pop_ld = apply(unique(res_gwas_ld[,c("SNP_B1", "icd")]), 1, function(x){
    print(x)
    return(gwas_mat_sub[which(gwas_mat_sub$id==x["SNP_B1"] & gwas_mat_sub$ICD10==x["icd"]),c("id","ICD10","OTH", "AS","EAS","CHB","SAS")])
})
res_gwas_pop_ld = unique(do.call(rbind, res_gwas_pop_ld))
library(dplyr)
res_gwas_pop_ld = res_gwas_pop_ld %>% group_by(id, ICD10) %>%
    summarise(OTH = paste(sort(unique(OTH)), collapse = ""),
        AS = paste(sort(unique(AS)), collapse = ""),
        EAS = paste(sort(unique(EAS)), collapse = ""),
        CHB = paste(sort(unique(CHB)), collapse = ""),
        SAS = paste(sort(unique(SAS)), collapse = ""))

res_gwas_pop_ld = merge(res_gwas_ld, res_gwas_pop_ld, by.x=c("SNP_B1","icd"), by.y=c("id","ICD10"))
save(res_gwas_pop_ld, file="res_gwas_ld_fraction_pop.RData")

# observed MAF, res_MAF_M for mother and res_MAF_C for child
load("all_R2_0.5_MAF_in_motherICD_childICD.RData")

library(genpwr)

case_num = colSums(clin_status, na.rm=T)
case_control_num = colSums(!is.na(clin_status), na.rm=T)
case_rate = case_num/case_control_num

library(data.table)
library(stringr)
library(parallel)

pow = lapply(unique(res_gwas_pop_ld$icd), function(x){
    res = unique(res_gwas_pop_ld[res_gwas_pop_ld$icd==x,c("MarkerID", "OR")])
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
    res$icd = x
    res$power = unlist(pow)
    return(res)
})
pow = do.call(rbind, pow)

save(res_gwas_pop_ld, pow, file="res_gwas_ld_fraction_pop.RData")

res_gwas_pop_ld = merge(res_gwas_pop_ld, pow[,-2], by = c("MarkerID", "icd"))

get_fraction2_pop_ld = function(min_obMAF, min_power, pop, min_Rsq){
    if(pop == "ALL"){
        res = res_gwas_pop_ld[which(res_gwas_pop_ld$MarkerID %in% rownames(res_MAF_M)[res_MAF_M$MAF > min_obMAF] & res_gwas_pop_ld$power > min_power & res_gwas_pop_ld$R2 > min_Rsq), ]
    }else{
        res = res_gwas_pop_ld[which(res_gwas_pop_ld$MarkerID %in% rownames(res_MAF_M)[res_MAF_M$MAF > min_obMAF] & res_gwas_pop_ld$power > min_power & res_gwas_pop_ld$R2 > min_Rsq & res_gwas_pop_ld[, ..pop] == pop), ]
    }
    ## number of snp-ICD10 pairs in GWAS
    a = nrow(unique(res[,1:2]))
    ## number of snp-ICD10 pairs replicated in study
    b = nrow(unique(res[res$p.value<0.05,1:2]))
    return(c(min_obMAF = min_obMAF, min_power = min_power, pop = pop, min_Rsq = min_Rsq, snpICD10_all_GWAS = a, snpICD10_in_GWAS = b))
}

res_gwas_fraction2_pop_ld = do.call(rbind, lapply(c(0, 0.01, 0.05, 0.1), function(min_obMAF){
    return(do.call(rbind, lapply(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), function(min_pow){
        return(do.call(rbind, lapply(c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), function(min_Rsq){
            return(do.call(rbind, lapply(c("AS","EAS","CHB","SAS","OTH", "ALL"), function(pop){
                return(get_fraction2_pop_ld(min_obMAF, min_pow, pop, min_Rsq))
            })))
        })))
    })))
}))

save(res_gwas_pop_ld, res_gwas_fraction2_pop_ld, file="res_gwas_ld_fraction_pop.RData") #data3

res_gwas_fraction2_pop_ld = as.data.frame(res_gwas_fraction2_pop_ld, stringsAsFactors=F)
res_gwas_fraction2_pop_ld$min_obMAF = paste0("obMAF>", res_gwas_fraction2_pop_ld$min_obMAF)
res_gwas_fraction2_pop_ld$min_power = res_gwas_fraction2_pop_ld$min_power
res_gwas_fraction2_pop_ld$min_Rsq = res_gwas_fraction2_pop_ld$min_Rsq
res_gwas_fraction2_pop_ld$snpICD10_all_GWAS = as.numeric(res_gwas_fraction2_pop_ld$snpICD10_all_GWAS)
res_gwas_fraction2_pop_ld$snpICD10_in_GWAS = as.numeric(res_gwas_fraction2_pop_ld$snpICD10_in_GWAS)
res_gwas_fraction2_pop_ld$fraction = (res_gwas_fraction2_pop_ld$snpICD10_in_GWAS/res_gwas_fraction2_pop_ld$snpICD10_all_GWAS)*100

plot_hist_density = function(x, breaks = 10){
    res = density(x)
    hist(x, breaks, prob=T, ylim = c(0, max(res$y)), main = "", xlab = "Position to known risk loci")
    lines(res,lwd=2)
    # text(res$x[which.max(res$y)], max(res$y)*1.02, labels = signif(res$x[which.max(res$y)],3))
}

res_gwas_pop_ld$power_levels = "0-0.2"
res_gwas_pop_ld$power_levels[which(res_gwas_pop_ld$power>0.2)] = "0.2-0.4"
res_gwas_pop_ld$power_levels[which(res_gwas_pop_ld$power>0.4)] = "0.4-0.6"
res_gwas_pop_ld$power_levels[which(res_gwas_pop_ld$power>0.6)] = "0.6-0.8"
res_gwas_pop_ld$power_levels[which(res_gwas_pop_ld$power>0.8)] = "0.8-1"
# res_gwas_pop_ld$power_levels = factor(res_gwas_pop_ld$power_levels, levels = c(">0", ">0.2", ">0.4", ">0.6", ">0.8"))
res_gwas_pop_ld$R2_levels = "0-0.2"
res_gwas_pop_ld$R2_levels[res_gwas_pop_ld$R2>0.2] = "0.2-0.4"
res_gwas_pop_ld$R2_levels[res_gwas_pop_ld$R2>0.4] = "0.4-0.6"
res_gwas_pop_ld$R2_levels[res_gwas_pop_ld$R2>0.6] = "0.6-0.8"
res_gwas_pop_ld$R2_levels[res_gwas_pop_ld$R2>0.8] = "0.8-1"
# res_gwas_pop_ld$R2_levels = factor(res_gwas_pop_ld$R2_levels, levels = c(">0", ">0.2", ">0.4", ">0.6", ">0.8"))

res_gwas_pop_ld$replicated = ifelse(res_gwas_pop_ld$p.value < 0.05, "Replicated","Not replicated")
res_gwas_pop_ld$replicated_pop = res_gwas_pop_ld$replicated
res_gwas_pop_ld$replicated_pop[res_gwas_pop_ld$CHB == "CHB" & res_gwas_pop_ld$replicated_pop == "Replicated"] = "CHB"
res_gwas_pop_ld$replicated_pop[res_gwas_pop_ld$EAS == "EAS" & res_gwas_pop_ld$replicated_pop == "Replicated"] = "EAS"
res_gwas_pop_ld$replicated_pop[res_gwas_pop_ld$SAS == "SAS" & res_gwas_pop_ld$replicated_pop == "Replicated"] = "SAS"
res_gwas_pop_ld$replicated_pop[res_gwas_pop_ld$AS == "AS" & res_gwas_pop_ld$replicated_pop == "Replicated"] = "AS"
res_gwas_pop_ld$replicated_pop[res_gwas_pop_ld$OTH == "OTH" & res_gwas_pop_ld$replicated_pop == "Replicated"] = "OTH"
res_gwas_pop_ld$P = -log10(res_gwas_pop_ld$p.value)


res_gwas_pop_ld$power_levels1 = rowMeans(apply(str_split(res_gwas_pop_ld$power_levels, pattern="-",simplify=T),2,as.numeric))
res_gwas_pop_ld$ICDanno = ""
res_gwas_pop_ld$ICDanno[grep("[AB]", res_gwas_pop_ld$icd)] = "Certain infectious and parasitic diseases"
res_gwas_pop_ld$ICDanno[which(res_gwas_pop_ld$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24", "D25", "D26"))] = "Neoplasms"
res_gwas_pop_ld$ICDanno[which(res_gwas_pop_ld$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
res_gwas_pop_ld$ICDanno[grep("E", res_gwas_pop_ld$icd)] = "Endocrine, nutritional and metabolic diseases"
res_gwas_pop_ld$ICDanno[grep("F", res_gwas_pop_ld$icd)] = "Mental and behavioural disorders"
res_gwas_pop_ld$ICDanno[which(res_gwas_pop_ld$icd %in% c("H00", "H01", "H02", "H04", "H10","H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
res_gwas_pop_ld$ICDanno[which(res_gwas_pop_ld$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92","H93"))] = "Ear and mastoid process"
res_gwas_pop_ld$ICDanno[grep("I", res_gwas_pop_ld$icd)] = "Circulatory system"
res_gwas_pop_ld$ICDanno[grep("J", res_gwas_pop_ld$icd)] = "Respiratory system"
res_gwas_pop_ld$ICDanno[grep("K", res_gwas_pop_ld$icd)] = "Digestive system"
res_gwas_pop_ld$ICDanno[grep("L", res_gwas_pop_ld$icd)] = "Skin and subcutaneous tissue"
res_gwas_pop_ld$ICDanno[grep("M", res_gwas_pop_ld$icd)] = "Musculoskeletal system and connective tissue"
res_gwas_pop_ld$ICDanno[grep("N", res_gwas_pop_ld$icd)] = "Genitourinary system"
res_gwas_pop_ld$ICDanno[grep("O", res_gwas_pop_ld$icd)] = "Pregnancy, childbirth and the puerperium"
res_gwas_pop_ld$ICDanno[grep("P", res_gwas_pop_ld$icd)] = "Certain conditions originating in the perinatal period"
res_gwas_pop_ld$ICDanno[grep("Q", res_gwas_pop_ld$icd)] = "Congenital malformations, deformations and chromosomal abnormalities"
colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
    "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
    "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

res_gwas_pop_ld$R2_levels1 = factor(ifelse(res_gwas_pop_ld$R2>0.8, "high", "low"), levels = c("low", "high"))

pdf("Figure2B.pdf", width = 10, height = 2.4)
ggplot(res_gwas_pop_ld, aes(power_levels, -log10(p.value), fill=R2_levels1)) +
    geom_boxplot(outlier.size = 0.3, outlier.shape = NA, width = 0.6) +
    scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2.5)) +
    facet_grid(~replicated_pop) +
    scale_fill_manual(values = rev(pal_ucscgb("default")(5))[c(1,3)]) + theme_bw() +
    geom_hline(aes(yintercept=-log10(5e-8)), linetype = "dashed") +
    geom_hline(aes(yintercept=-log10(1e-3)), linetype = "dashed") +
    geom_hline(aes(yintercept=-log10(0.05)), linetype = "dashed") +
    theme(aspect.ratio = 1)
dev.off()
```

## Figure2C

```R
library(data.table)
load("phewas_power_P0.001.RData")
res_sig = snp_power[which(snp_power$power > 0.8), ]

all_SAIGE_snp = fread("all_SAIGE_snp.txt.gz", header=F, sep="\t",stringsAsFactors=F)
imputed_info = fread("snp_imputed_R0.5_MAF0.01.txt",header=T,sep="\t",stringsAsFactor=F)

######## ldblock
ldproxy <- fread("ld_results.txt.gz", header=T,sep="\t",stringsAsFactors=F)
ldproxy = ldproxy[which(ldproxy$SNP_A %in% intersect(imputed_info$SNP, all_SAIGE_snp$V1)), ]
ldproxy$SNP_A1 = gsub(":[atcgATCG]+", "", ldproxy$SNP_A)
ldproxy$SNP_B1 = gsub(":[atcgATCG]+", "", ldproxy$SNP_B)

ldproxy_sig = ldproxy[ldproxy$SNP_A %in% res_sig$MarkerID & ldproxy$R2 > 0.8, ]
save(ldproxy_sig, file="ld_results_1e3_Rsq0.8.RData")

load("gwas_catalog_mat.RData")
gwas_mat = gwas_mat[-which(gwas_mat$CHR_ID %in% c("MT","X","Y")), ]
gwas_mat$id = paste(gwas_mat$CHR_ID, gwas_mat$CHR_POS, sep=":")
gwas_mat = as.data.frame(gwas_mat, stringsAsFactors=F)

res_OR = lapply(c(0.2,0.3,0.4,0.5,0.6,0.7,0.8), function(Rsq){
    all_num = 5901875   #intersect(all_SAIGE_snp$V1, imputed_info$SNP)
    if(Rsq > 0.2){
        ldproxy_sub = ldproxy[ldproxy$R2 >= Rsq, ]
    }else{
        ldproxy_sub = ldproxy
    }
    res_pop = do.call(rbind, lapply(c("CHB", "EAS", "AS", "ALL"), function(pop){
        if(pop=="ALL"){
            id = unique(gwas_mat$id)
        }else{
            id = unique(gwas_mat$id[gwas_mat[,pop] == pop])
        }
        all_in_gwas_snp = unique(ldproxy_sub$SNP_A[which(ldproxy_sub$SNP_B1 %in% id)])
        all_in_gwas_num = length(all_in_gwas_snp)

        res_P = do.call(rbind, lapply(c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9), function(pvalue){
                sig_snp = unique(res_sig$MarkerID[res_sig$p.value < pvalue])
                sig_num = length(sig_snp)

                sig_in_gwas_num = length(intersect(all_in_gwas_snp, sig_snp))

                print(paste(pop, pvalue, Rsq))
                FoE_RR = (sig_in_gwas_num/sig_num)/(all_in_gwas_num/all_num)    #relative risk
                Var_log_FoE_RR=1/sig_in_gwas_num - 1/sig_num + 1/(all_in_gwas_num - sig_in_gwas_num) - 1/(all_num - sig_num)
                FoE_RR_confintMin = exp(log(FoE_RR) - 1.96 * Var_log_FoE_RR)
                FoE_RR_confintMax = exp(log(FoE_RR) + 1.96 * Var_log_FoE_RR)
                res = fisher.test(matrix(c(sig_in_gwas_num, sig_num - sig_in_gwas_num, all_in_gwas_num - sig_in_gwas_num, all_num - sig_num - all_in_gwas_num + sig_in_gwas_num), nrow=2, byrow=T))
                return(c(pop = pop, maxP = pvalue, minRsq = Rsq, FoE_RR = FoE_RR, FoE_RR_confintMin = FoE_RR_confintMin, FoE_RR_confintMax = FoE_RR_confintMax,
                    fisher_pvalue = res$p.value, fisher_OR = as.numeric(res$estimate), fisher_OR_confintMin = res$conf.int[1], fisher_OR_confintMax = res$conf.int[2],
                    sig_in_gwas_num = sig_in_gwas_num, sig_num = sig_num, all_in_gwas_num = all_in_gwas_num, all_num = all_num))
            }))
        return(res_P)
    }))
    rm(ldproxy_sub)
    gc()
    return(res_pop)
})

res_OR = as.data.frame(do.call(rbind, res_OR), stringsAsFactors=F)
save(res_OR, file="res_OR.RData")

res_OR$pop = factor(res_OR$pop, levels = c("CHB","EAS","AS","ALL"), labels = c("CHB GWAS","EAS GWAS","AS GWAS","Total GWAS"))
res_OR$FoE_RR = as.numeric(res_OR$FoE_RR)
res_OR$FoE_RR_confintMin = as.numeric(res_OR$FoE_RR_confintMin)
res_OR$FoE_RR_confintMax = as.numeric(res_OR$fisher_OR_confintMax)
res_OR$fisher_pvalue = as.numeric(res_OR$fisher_pvalue)
res_OR$fisher_OR = as.numeric(res_OR$fisher_OR)
res_OR$fisher_OR_confintMin = as.numeric(res_OR$fisher_OR_confintMin)
res_OR$fisher_OR_confintMax = as.numeric(res_OR$fisher_OR_confintMax)
res_OR$sig_in_gwas_num = as.numeric(res_OR$sig_in_gwas_num)
res_OR$sig_num = as.numeric(res_OR$sig_num)
res_OR$all_in_gwas_num = as.numeric(res_OR$all_in_gwas_num)
res_OR$all_num = as.numeric(res_OR$all_num)
res_OR$maxP = factor(as.numeric(res_OR$maxP), levels=sort(as.numeric(unique(res_OR$maxP)), decreasing=T), labels = sort(signif(-log(as.numeric(unique(res_OR$maxP)), base=10),3)))
res_OR$minRsq = paste("minRsq", res_OR$minRsq)
res_OR$minBetaOR = paste("minBetaOR", res_OR$minBetaOR)

pdf("Figure2C.pdf")
ggplot(res_OR[which(res_OR$minRsq == "minRsq 0.8"), ],
        aes(maxP, FoE_RR, color=pop, group = pop)) + geom_point() + geom_line() + facet_grid(~minRsq, scales="free") +
    geom_ribbon(aes(ymin = FoE_RR_confintMin, ymax = FoE_RR_confintMax), color=NA, fill = "grey60", alpha = 0.3) +
    geom_vline(aes(xintercept=5.3), linetype = "dashed") +
    scale_y_sqrt(n.breaks = 10) +
    scale_color_d3() +
    labs(x="-log10(Pvalue)", y="Fold of enrichment", fill="") +
    theme_bw() +
    theme(legend.position="right",
        legend.title = element_blank(),
        axis.text.x = element_text(angle=0, vjust=0.5),
        aspect.ratio = 1)
dev.off()
```

## Figure2D

```R
load("pheWAS_sig_SAIGE_5e8_power0.8.RData")
library(UpSetR)
dat = unique(res_sig[, c("LDblock","gwas", "AS", "EAS","SAS","CHB","eQTL", "geneatlas", "phewas", "labwas")]) %>% group_by(LDblock) %>% summarise_all(sum) %>% as.data.frame
rownames(dat) = dat$LDblock
dat$LDblock = NULL
dat$ICD = 1
dat[dat>1] = 1
colnames(dat) = c("Total GWAS", "AS GWAS", "EAS GWAS", "SAS GWAS", "Chinese GWAS", "eQTL calalog", "EUR GeneATLAS","EUR PheWAS","EUR LabWAS", "Our PheWAS")

pdf("Figure2D.pdf", width= 6, heigh = 6)
upset(dat, sets = c("EUR LabWAS", "EUR PheWAS", "EUR GeneATLAS", "Chinese GWAS", "SAS GWAS", "EAS GWAS", "AS GWAS", "Total GWAS", "Our PheWAS"),
    nintersects = 20,
    order.by = "freq", keep.order = T, point.size = 2, line.size = 0.8, mb.ratio = c(0.6, 0.4),
    mainbar.y.label = "Intersection Size", sets.x.label = "Set size")
dev.off()
```
