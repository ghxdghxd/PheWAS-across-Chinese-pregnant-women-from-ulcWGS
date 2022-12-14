# 8. PheWAS for fetal variants

## 8.1 input

```R
files = list.files("SNPs", pattern="*.gz",full.names=T)
sample_ID <- read.table("all_sample_ID.txt", header=F,sep="\t",stringsAsFactors=F)

sample_ID = c(sample_ID1$V1, sample_ID2$V1)

clin_merge="ICD_child_mat.RData"
load(clin_merge)

library(dplyr)
ICD_mat = ICD_mat[,-c(1:4,6)] %>% group_by(ID) %>% summarise_all(sum)

FFY_RData = "cnv_addPREFACE.RData"
load(FFY_RData)
cnv_mat$PREFACE = cnv_mat$PREFACE * 100

int_sample = intersect(ICD_mat$ID, sample_ID)
ICD_mat = ICD_mat[which(ICD_mat$ID %in% int_sample),]
FF = cnv_mat$PREFACE[match(ICD_mat$ID, cnv_mat$ID)]
rm(cnv_mat)
gc()
index = match(ICD_mat$ID, sample_ID)

save(FF, ICD_mat, file="SNPs_RD_AD/ICD_FF.RData")

#########
merge_AD_RD_for_ICD = function(f, index){
    SNPs <- fread(f, header=T,sep="\t",stringsAsFactors=F)
    AD <- fread(paste0("AD/", gsub("SNPs","AD", basename(f))), header=F, sep="\t", select=SNPs$ind1, stringsAsFactors=F)
    RD <- fread(paste0("RD/", gsub("SNPs","RD", basename(f))), header=F, sep="\t", select=SNPs$ind1, stringsAsFactors=F)
    AD = AD[index]
    RD = RD[index]
    save(SNPs, AD, RD, file=paste0("SNPs_RD_AD/", gsub(".txt.gz", "_RD_AD.RData", basename(f))))
}

library(data.talbe)
for(f in files){
    merge_AD_RD_for_ICD(f, index)
}

#########
library(parallel)

get_AF_callNum = function(f){
    load(f)
    GT = AD
    GT[RD>0 & AD==0] = 0
    GT[RD>0 & AD>0] = 1
    GT[RD==0 & AD>0] = 2
    AF_callNum = lapply(1:ncol(GT),function(n){
        res = summary(snp(unlist(GT[,..n]), name.genotypes=c(0,1,2)))
        return(c(obAF_icd = res$allele.freq["B","percentage"]))
    })
    AF_callNum = as.data.frame(do.call(rbind, AF_callNum), stringsAsFactors = F)
    AF_callNum[, "obAF_icd"] = AF_callNum[,"obAF_icd"]/100
    AF_callNum$callNum_icd = colSums(!is.na(GT))
    AF_callNum$AF_sd_icd = sqrt((SNPs$AF*(1-SNPs$AF))/AF_callNum$callNum_icd)

    AF_callNum$AF_ref_min_icd = SNPs$AF_ref - 3 * AF_callNum$AF_sd_icd
    AF_callNum$AF_ref_min_icd[AF_callNum$AF_ref_min_icd < 0] = 0 #min(AF_min[which(AF_min>0)])
    AF_callNum$AF_ref_max_icd = SNPs$AF_ref + 3 * AF_callNum$AF_sd_icd
    AF_callNum$AF_ref_max_icd[AF_callNum$AF_ref_max_icd > 1] = 1

    SNPs = cbind(SNPs, AF_callNum)
    save(AD, RD, SNPs, file=gsub("SNPs_RD_AD", "SNPs_RD_AD_icd", f))
}

files = list.files("SNPs_RD_AD", pattern="*.RData", full.names=T)

library(SNPassoc)
library(foreach)
library(doParallel)

registerDoParallel(10)
foreach(a = 1:length(files)) %dopar% {
    get_AF_callNum_HWE(files[a])
}
```

## 8.2 run.R

```R
library(parallel)
library(emdbook)
library(invgamma)
library(data.table)
Posterior <- function(par, nA_1, nB_1, nA_2, nB_2, p_eas, Pi1, Pi2) {
    # cat(par, "\n")
    n1 <- length(nA_1)
    n2 <- length(nA_2)
    delta <- par[1]
    phi   <- par[2]
    # pi    <- par[3]
    mju <- par[3]
    theta <- par[4]
    # theta = 1e-5
    p <- par[5]
    alpha <- par[6]
    sigma <- par[7]
    Omega <- matrix(c(1-delta, delta, delta, 1-delta), 2, 2) %*% diag(c(2*(1-phi), 2*phi))
    N1 <- cbind(nA_1, nB_1) %*% Omega
    N2 <- cbind(nA_2, nB_2) %*% Omega
    N_A1 <- round(sum(N1[, 1]))
    N_B1 <- round(sum(N1[, 2]))
    N_A2 <- round(sum(N2[, 1]))
    N_B2 <- round(sum(N2[, 2]))
    loglik <- dbetabinom(x = N_A1, size = N_A1 + N_B1, prob = p + Pi1*mju, theta, log = TRUE) +
                dbetabinom(x = N_A2, size = N_A2 + N_B2, prob = alpha*p + Pi2*mju, theta, log = TRUE)
    PP <- sum(loglik, na.rm = T) +
        #   dbeta(pi, shape1 = 0.035, shape2 = 1, log = TRUE) +
        #   dnorm(pi, mean = 0.035, log = TRUE) +
          dbeta(delta, shape1 = 1.01, shape2 = 1.99, log = TRUE) +
          dgamma(theta, shape = 1, rate = 0.1, log = TRUE) +
          dbeta(phi, shape1 = 10, shape2 = 10, log = TRUE) +
          dbeta(p, shape1 = 1, shape2 = 1, log = TRUE) +
        #   dnorm(p, mean = p_eas, log = TRUE) +
          dnorm(log(alpha), mean = 0, sd = sigma, log = TRUE) +
          dinvgamma(sigma^2, shape = 1, scale = 1, log = TRUE) +
          dbeta(mju, shape1 = 1E-8, shape2 = 1, log = TRUE)
    # print(PP)
    return(-PP)
}

NIPT_ICD = function(ICD, RD_snp, AD_snp, AF, AF_min, AF_max, FF){
    dat = cbind(RD_snp, AD_snp, ICD, FF)
    colnames(dat) = c("RD","AD","ICD", "FF")
    dat = na.omit(dat)
    if(nrow(dat) == 0 || sum(dat$ICD, na.rm=T) == 0){
        return(rep(NA, 13))
    }
    nA_1 = dat$RD[dat$ICD == 0]        # ref alleles count ctrl
    nB_1 = dat$AD[dat$ICD == 0]        # alt alleles count ctrl
    nA_2 = dat$RD[dat$ICD > 0]        # ref alleles count case
    nB_2 = dat$AD[dat$ICD > 0]        # alt alleles count case

    Pi1 = mean(dat$FF[dat$ICD == 0], na.rm = T)
    Pi2 = mean(dat$FF[dat$ICD > 0], na.rm = T)

    delta = 1E-3        # MQ > 30 = -10 * log10(1e-3)
    delta_max = 1e-4
    phi = 0.5           # [0.4, 0.6]
    # Pi = 0.035           # [0.01, 0.05]
    mju = 1e-8        # [1e-8, 1e-4],
    mju_max = 1e-3
    theta = 1e-5          # [1e-6, 1e-5]
    alpha = 1           # [1E-6, Inf], ratio of AF betweed case and ctrl
    alpha_min = 0.05 # 1e-6,,,,0.05/0.95 <alpha < 0.95/0.05
    alpha_max = min(c((1 - Pi2 * mju_max)/AF_max, 19))     # p_case = alpha*p + Pi2*mju, p_case <= 1;
    # alpha_max = Inf
    p = sum(c(nA_1, nA_2))/sum(c(nA_1, nA_2)+c(nB_1, nB_2)) # [0.01, 0.99], reference allele frequency
    sigma = 900
    a = try(optim(par = c(delta, phi, mju, theta, p, alpha, sigma),
                        fn = Posterior, nA_1 = nA_1, nB_1 = nB_1, nA_2 = nA_2, nB_2 = nB_2, p_eas = AF, Pi1 = Pi1, Pi2 = Pi2,
                        method = "L-BFGS-B",
                        lower = c(1E-32, 0.4, 1E-8, 1e-6, AF_min, alpha_min,      500),
                        upper = c(1e-3,  0.6, 1E-3, 1e-5, AF_max, alpha_max, Inf)))
    if(class(a) == "try-error"){
        return(rep(NA, 13))
    }
    a = c(a$par, length(nA_1), length(nA_2), p, AF, Pi1, Pi2)
    names(a) = c("delta", "phi", "mju", "theta", "p", "alpha", "sigma", "ctrl_num", "case_num", "p_init", "p_eas", "Pi1", "Pi2")
    return(a)
}

run_NIPT_ICD = function(f, ICD_mat, np = 1, nsamples = NULL, FF){
    load(f)
    load(paste0("filter_SNP_with_caseNum/", gsub("SNPs_RD_AD", "snps.caseNum30", basename(f))))
    if(is.null(nsamples)){
        index = 1:nrow(snps)
    }else{
        index = sample(1:nrow(snps), nrow(snps) * nsamples)
    }
    ICD_res = mclapply(index, function(x){
        ind = which(SNPs$id == snps$snp[x])
        icd_list = strsplit(snps$icd[x], split=",")[[1]]
        icd_res = lapply(icd_list, function(icd){
            print(paste(snps$snp[x], icd))
            res = NIPT_ICD(ICD_mat[, icd], RD[, ..ind], AD[, ..ind], SNPs$AF_ref[ind], SNPs$AF_ref_min_icd[ind], SNPs$AF_ref_max_icd[ind], FF)
            return(res)
        })
        names(icd_res) = icd_list
        icd_res = as.data.frame(do.call(rbind, icd_res), stringsAsFactors=F)
        icd_res$icd = rownames(icd_res)
        icd_res$snp = snps$snp[x]
        rownames(icd_res) = NULL
        return(icd_res)
    }, mc.cores = np)
    ICD_res = rbindlist(l=ICD_res)
    ICD_res = ICD_res[!is.na(ICD_res$p),]
    return(ICD_res)
}

clin="ICD_FF.RData"
load(clin_merge)
ICD_mat = ICD_mat[, which(colSums(ICD_mat>0) > 30)]

files = list.files("SNPs_RD_AD_icd", pattern="*.SNPs_RD_AD.RData",full.names=T)

for(f in files){
    out_file = gsub("SNPs_RD_AD","ICD_res", basename(f))
    if(!file.exists(out_file)){
        ICD_res = run_NIPT_ICD(f, ICD_mat, np = 30, nsamples = NULL, FF)
        save(ICD_res, file=out_file)
    }
}
```

## 8.3 母亲的显著位点

> /data/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge_Model/cohort1_cohort2_alpla0.05-19/ICD_res/for_SAIGE

```R
P_value <- function(cdf_function, x, paramet, side=0){
  p <- cdf_function(x,paramet[1],paramet[2])
  if(side < 0) p
  else if(side > 0) 1-p
  else
    if(p < 1/2) 2*p
    else 2*(1-p)
}

get_ICDres = function(dir){
    ICD_res = lapply(list.files(dir, pattern="*.RData", full.names = T), function(x){
        load(x)
        return(ICD_res)
    })
    ICD_res = na.omit(do.call(rbind, ICD_res))
    if(nrow(ICD_res)>0){
        ICD_res$pvalue = sapply(log(ICD_res$alpha), function(x){P_value(pnorm, x, paramet=c(0, 1))})
        return(ICD_res)
    }else{
        return(rep(NA, 16))
    }
}

library(ggplot2)
ggQQplot <- function(pvalues, title = NULL, random = NULL){
    pvalues = sort(na.omit(as.numeric(pvalues)))
    pvalues = cbind(op = pvalues, ep = ppoints(length(pvalues)))
    if(!is.null(random)){
        pvalues = pvalues[sample(1:nrow(pvalues), random), ]
    }
    lab <- quantile(pvalues[,"op"], 0.475)/quantile(pvalues[,"ep"], 0.475)
    pvalues <- -log10(pvalues)
    maxP = ceiling(max(pvalues))
    if(is.null(title)){
        title = ""
    }
    p = ggplot(as.data.frame(pvalues)) +
        geom_point(aes(ep, op), size = 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, maxP)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, maxP)) +
        geom_abline(intercept = 0, slope = 1, alpha = 1) +
        labs(x = "Expected (-log10(P))", y = "Observed (-log10(P))", title = title) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 15),
            axis.text = element_text(hjust = 0, size = 10),
             plot.margin=unit(c(1,2,1,1),"lines")) +
        annotate("text", x = 0.7 * maxP, y = 0.2 * maxP, size = 5,
            label = paste0("lambda = ", signif(lab, 5)))
    return(p)
}

plot_alpha_qqplot_AF = function(dat, icd, title=NULL, filter = F){
    p0 = ggplot(dat) + geom_histogram(aes(log(alpha)), bins = 100)
    if(filter){
        dat = dat[which(dat$snp %in% SNPs$id[which(SNPs$obAF >= SNPs$AF_min & SNPs$obAF <= SNPs$AF_max)] & dat$p_init > 0.01 & dat$case_num > 50 & dat$ctrl_num > 50 & dat$p_eas > 0.05 & dat$p_eas < 0.95), ]
        title = paste(icd, title, "filtered")
    }else{
        title = paste(icd, title)
    }
    if(nrow(dat)==0){
        return(cowplot::plot_grid(p0, NULL, NULL, NULL, nrow=1))
    }
    p1 = ggplot(dat) + geom_histogram(aes(log(alpha)), bins = 100)
    if(length(dat$pvalue) > 20000){
        random_sample=10000
        index = sample(1:nrow(dat), random_sample)
    }else{
        random_sample = NULL
        index = 1:nrow(dat)
    }
    p2 = ggQQplot(dat$pvalue, title = "", random = NULL)
    dat$logAlpha = log(dat$alpha)
    p3 = ggplot(dat[index,]) + geom_point(aes(p_init, p, color=alpha)) + scale_color_gradient2() + labs(title = title)
    # p4 = ggplot(dat[index,]) + geom_point(aes(p_init, p, color=logAlpha)) + scale_color_gradient2() + xlim(0.1,1)
    return(cowplot::plot_grid(p0, p1, p2, p3, nrow=1))
}

files = list.files("ICD_res", pattern="*.RData", full.names=T))

ICD_res = lapply(files, function(x){load(x);return(ICD_res)})

ICD_res_filter = lapply(ICD_res, function(m){
    return(m[which(m$snp %in% SNPs$id[which(SNPs$obAF >= SNPs$AF_min & SNPs$obAF <= SNPs$AF_max)] & m$p_init > 0.01 & m$case_num > 50 & m$ctrl_num > 50 & m$p_eas > 0.05 & m$p_eas < 0.95 & m$alpha > 0.05 & m$alpha < 19),])
})
ICD_res_filter = do.call(rbind, ICD_res_filter)
ICD_res_filter$pvalue_norm = sapply(log(ICD_res_filter$alpha), function(x){P_value(pnorm, x, paramet=c(0, 1))})
ICD_res_filter$pvalue_test = sapply(log(ICD_res_filter$alpha), function(x){P_value(pnorm, x, paramet=c(-0.1654, 1.179))})
save(ICD_res_filter, file="all_ICD_res_filter.RData")

ICD_res_filter_sig = ICD_res_filter[which(ICD_res_filter$pvalue_norm<0.05), ]
write.table(ICD_res_filter_sig, file="ICD_res_filter_sig.txt", r=F, c=T, sep="\t", quote=F)


library(data.table)
library(stringr)
load("/share/data7/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge/avsnp150_ChainMAP_genos/all/raw/GT_for_Eagle/pheWAS_sig_merge/for_SAIGE/phewas_power_P0.001.RData")
res_sig_raw = snp_power[which(snp_power$power > 0.8), ]

ICD_res_filter_withMother = ICD_res_filter[which(ICD_res_filter$snp %in% res_sig_raw$MarkerID),]
save(ICD_res_filter_withMother, file="all_ICD_res_filter_withMother.RData")

ICD_res_filter_withMother_sig = ICD_res_filter_withMother[ICD_res_filter_withMother$pvalue_norm < 0.05, ]
write.table(ICD_res_filter_withMother_sig, file="ICD_res_filter_withMother_sig.txt", r=F, c=T, sep="\t", quote=F)

pdf("qqplot_all_pvalue_norm.pdf")
ggQQplot(ICD_res_filter_withMother$pvalue_norm, title="all pvalue_norm", random=10000)
dev.off()

pdf("ICD_res_filter_withMother.pdf", width = 12, height = 3)
plot_alpha_qqplot_AF(ICD_res_filter_withMother[ICD_res_filter_withMother$icd=="D18",], icd = "D18")
plot_alpha_qqplot_AF(ICD_res_filter_withMother[ICD_res_filter_withMother$icd=="J06",], icd = "J06")
plot_alpha_qqplot_AF(ICD_res_filter_withMother[ICD_res_filter_withMother$icd=="B34",], icd = "B34")
plot_alpha_qqplot_AF(ICD_res_filter_withMother[ICD_res_filter_withMother$icd=="Q82",], icd = "Q82")
dev.off()

```

### 统计显著突变数目 ICD_res_filter_withMother_sig

```R
library(data.table)
library(dplyr)
library(stringr)

ICD_anno <- read.csv("/share/data7/NIPT/bam_vcf/bvcfs_var_norm_all_unique_genotype_merge/avsnp150_ChainMAP_genos/all/raw/GT_for_Eagle/pheWAS_sig_merge/1e4_BMI/ICD10.en.zh.txt", header=T,sep="\t",stringsAsFactors=F)

ICD_res_filter_withMother_sig = ICD_res_filter_withMother_sig[-grep("R", ICD_res_filter_withMother_sig$icd), ]

ICD_res_filter_withMother_sig$ICDanno = ""
ICD_res_filter_withMother_sig$ICDanno[grep("[AB]", ICD_res_filter_withMother_sig$icd)] = "Certain infectious and parasitic diseases"
# ICD_res_filter_withMother_sig$ICDanno[which(ICD_res_filter_withMother_sig$icd %in% c("C34", "C50", "C73", "D06", "D17", "D18", "D22", "D24", "D25", "D26"))] = "Neoplasms"
# ICD_res_filter_withMother_sig$ICDanno[which(ICD_res_filter_withMother_sig$icd %in% c("D50", "D56", "D72"))] = "Blood and blood-forming organs and certain disorders involving the immune mechanism"
ICD_res_filter_withMother_sig$ICDanno[grep("E", ICD_res_filter_withMother_sig$icd)] = "Endocrine, nutritional and metabolic diseases"
# ICD_res_filter_withMother_sig$ICDanno[grep("F", ICD_res_filter_withMother_sig$icd)] = "Mental and behavioural disorders"
ICD_res_filter_withMother_sig$ICDanno[which(ICD_res_filter_withMother_sig$icd %in% c("H00", "H01", "H02", "H04", "H10","H11", "H16", "H20", "H35", "H40", "H43", "H52"))] = "Eye and adnexa"
ICD_res_filter_withMother_sig$ICDanno[which(ICD_res_filter_withMother_sig$icd %in% c("H60","H61","H65", "H66", "H69", "H72", "H81", "H90", "H91", "H92","H93"))] = "Ear and mastoid process"
# ICD_res_filter_withMother_sig$ICDanno[grep("I", ICD_res_filter_withMother_sig$icd)] = "Circulatory system"
ICD_res_filter_withMother_sig$ICDanno[grep("J", ICD_res_filter_withMother_sig$icd)] = "Respiratory system"
ICD_res_filter_withMother_sig$ICDanno[grep("K", ICD_res_filter_withMother_sig$icd)] = "Digestive system"
ICD_res_filter_withMother_sig$ICDanno[grep("L", ICD_res_filter_withMother_sig$icd)] = "Skin and subcutaneous tissue"
# ICD_res_filter_withMother_sig$ICDanno[grep("M", ICD_res_filter_withMother_sig$icd)] = "Musculoskeletal system and connective tissue"
# ICD_res_filter_withMother_sig$ICDanno[grep("N", ICD_res_filter_withMother_sig$icd)] = "Genitourinary system"
# ICD_res_filter_withMother_sig$ICDanno[grep("O", ICD_res_filter_withMother_sig$icd)] = "Pregnancy, childbirth and the puerperium"
ICD_res_filter_withMother_sig$ICDanno[grep("P", ICD_res_filter_withMother_sig$icd)] = "Certain conditions originating in the perinatal period"
ICD_res_filter_withMother_sig$ICDanno[grep("Q", ICD_res_filter_withMother_sig$icd)] = "Congenital malformations, deformations and chromosomal abnormalities"

res_sig = ICD_res_filter_withMother_sig

# ICD大类的snp数目
plot_ICDanno = function(res_sig, title){
    ndat_ICDanno = as.data.frame(table(unique(res_sig[,c("snp","ICDanno")])$ICDanno))
    ndat_ICDanno$Var1 = factor(ndat_ICDanno$Var1, levels = rev(ndat_ICDanno$Var1[order(ndat_ICDanno$Freq)]))
    p_ICDanno = ggplot(ndat_ICDanno) + geom_col(aes(Var1, Freq)) +
        scale_y_sqrt() + labs(x="",title=title) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10, angle = -90, hjust = 0, vjust = 0.5),
            aspect.ratio = 1)
    return(p_ICDanno)
}
p_ICDanno = plot_ICDanno(res_sig, title = "risk loci of ICDanno")
# p_ICDanno_ld = plot_ICDanno(res_sig_ld, title = "risk LDblock of ICDanno")

pdf("Model_ICDanno_snp_count.pdf", width=3, height=8)
p_ICDanno
# p_ICDanno_ld
dev.off()

# ICD小类的snp数目
plot_icd = function(res_sig, title, minFreq){
    ndat_icd = as.data.frame(table(unique(res_sig[,c("snp","icd")])$icd))
    ndat_icd$Var1 = factor(ndat_icd$Var1, levels = rev(ndat_icd$Var1[order(ndat_icd$Freq)]))
    p_icd = ggplot(ndat_icd[which(ndat_icd$Freq > minFreq),]) + geom_col(aes(Var1, Freq)) +
        scale_y_sqrt() + labs(x="",title=title) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10, angle=90, vjust=0.5, hjust=1))
    return(p_icd)
}

p_icd = plot_icd(res_sig, title = "risk loci of icd", minFreq = 0)
# p_icd_ld = plot_icd(res_sig_ld, title = "risk LDblock of icd", minFreq = 10)

pdf("Model_icd_snp_count.pdf", width=9, height=3)
p_icd
dev.off()

# 每个snp影响小类icd的数目分布
plot_snp = function(res_sig, title){
    ndat_snp = as.data.frame(table(table(unique(res_sig[,c("snp","icd")])$snp)))
    p_snp = ggplot(ndat_snp) + geom_col(aes(Var1, Freq)) +
        scale_y_sqrt() + labs(x="",title=title) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10))
    return(p_snp)
}
p_snp = plot_snp(res_sig, title = "risk loci of icd")
# p_snp_ld = plot_snp(res_sig_ld, title = "risk LDblock of icd")
pdf("Model_snp_icd_count.pdf", width=6, height=3)
p_snp
# p_snp_ld
dev.off()

# 每个snp影响大类icd的数目分布
plot_snp_ICDanno = function(res_sig, title){
    ndat_snp_ICDanno = as.data.frame(table(table(unique(res_sig[,c("snp","ICDanno")])$snp)))
    ndat_snp_ICDanno$Var1 = factor(ndat_snp_ICDanno$Var1, levels = rev(ndat_snp_ICDanno$Var1[order(ndat_snp_ICDanno$Freq)]))
    p_snp_ICDanno = ggplot(ndat_snp_ICDanno) + geom_col(aes(Var1, Freq)) +
        scale_y_sqrt() + labs(x="",title=title) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10))
    return(p_snp_ICDanno)

    ggplot(ndat_snp_ICDanno) + geom_col(aes(Var1, Freq)) +
        scale_y_sqrt() + labs(x="",title=title) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10))

}
p_snp_ICDanno = plot_snp_ICDanno(res_sig, title = "risk loci of icd")
# p_snp_ICDanno_ld = plot_snp_ICDanno(res_sig_ld, title = "risk LDblock of icd")

pdf("Model_snp_ICDanno_count.pdf", width=3, height=3)
p_snp_ICDanno
# p_snp_ICDanno_ld
dev.off()

# 每个chr中影响小类icd的snp数目分布
plot_chr = function(res_sig, title){
    ndat_chr = as.data.frame(table(str_split(unique(res_sig[,c("snp","icd")])$snp, pattern=":", simplify=T)[,1]))
    ndat_chr$Var1 = factor(ndat_chr$Var1, levels = rev(ndat_chr$Var1[order(ndat_chr$Freq)]))
    p_chr = ggplot(ndat_chr) + geom_col(aes(Var1, Freq)) +
        scale_y_sqrt() + labs(x="Chromosome", title = title) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10))
    return(p_chr)
    return(p_chr)
}

p_chr = plot_chr(res_sig, title = "risk loci of icd")
# p_chr_ld = plot_chr(res_sig_ld, title = "risk LDblock of icd")

pdf("Model_chr_icd_count.pdf", width=6, height=3)
p_chr
# p_chr_ld
dev.off()

library(stringr)
Phenogram_input = unique(res_sig[,c("snp","icd","ICDgroup")])
Phenogram_input = cbind(Phenogram_input, str_split(Phenogram_input$snp, pattern=":", simplify=T)[,1:2])
colnames(Phenogram_input) = c("snp","phenotype","colorgroup","chr","pos")
write.table(Phenogram_input, file="Phenogram_plot_input.txt", r=F,c=T,sep="\t",quote=F)
```

```sh
ruby ../05Eagle_miniMac4_PheWAS/pheWAS/PhenoGram_script/pheno_gram.rb -g ../05Eagle_miniMac4_PheWAS/pheWAS/genome_file_hg38.txt \
    -i Phenogram_plot_input.txt -Y ../05Eagle_miniMac4_PheWAS/pheWAS/cytoband_hg38.txt --shade-chromatin -f pdf --restrict-chroms -S medium -p standard -o pheno_gram -F -N -T -z -c group

for i in `sed 1d Phenogram_plot_input.txt|cut -f 3|sort -u`;
do
    awk 'NR==1||$3=="'$i'"{print}' Phenogram_plot_input.txt > test.txt
    ruby ../05Eagle_miniMac4_PheWAS/pheWAS/PhenoGram_script/pheno_gram.rb -g ../05Eagle_miniMac4_PheWAS/pheWAS/genome_file_hg38.txt -i test.txt -Y ../05Eagle_miniMac4_PheWAS/pheWAS/cytoband_hg38.txt --shade-chromatin -f pdf --restrict-chroms -S medium -p standard -o pheno_gram.$i -F -N -T -z
    rm test.txt
done
```

```R
library(ggbreak)
plot_snp_levels = function(res_sig){
    ndat_snp_levels = reshape2::melt(table(unique(res_sig[,c("MarkerID","icd","levels")])[,2:3]))
    icds = ndat_snp_levels %>% group_by(icd) %>% summarise(n = sum(value)) %>% arrange(desc(n)) %>% select(icd) %>% unlist %>% as.character
    ndat_snp_levels = ndat_snp_levels[ndat_snp_levels$icd %in% icds, ]
    ndat_snp_levels$icd = factor(ndat_snp_levels$icd, levels = icds)
    # ndat_snp_levels$levels = factor(ndat_snp_levels$levels, levels =c("L4", "L3", "L2", "L1"))
    p_snp_levels = ggplot(ndat_snp_levels) + geom_col(aes(icd, value, fill=levels)) +
        theme_bw() + labs(x="",y="Associations pairs' count") +
        scale_fill_npg() +
        scale_y_sqrt() +
        scale_y_break(c(100,150), scale='free') +
        theme(axis.text.x = element_text(size=10, angle=-90, vjust=0.5, hjust=1))
    return(p_snp_levels)
}
p_snp_ld_levels = plot_snp_levels(res_sig_ld)
pdf("pheWAS_snp_icd_count.levels.1.pdf", width=8, height=4)
p_snp_ld_levels
dev.off()


# 总ICD位点图
colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
    "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
    "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

library(ggbreak)
# Cairo::CairoPNG("all_man_ld.png", width=1500, height=300)
pdf("all_man_ld.pdf", width = 10, heigh = 10)
res_sig[,c("icd","ICDanno","LDblock", "p.value", "OR")] %>% group_by(icd, ICDanno, LDblock) %>% summarise(P = min(p.value), OR = OR[which.min(p.value)]) %>%
    ggplot() + geom_point(aes(icd, -log10(P), color = ICDanno, size = OR), alpha = 0.8) +
    facet_grid(~ICDanno, scale="free", space="free", switch="x") +
    scale_size(range = c(0, 4)) +
    scale_color_manual(values = colors37[1:17], guide = "none") + #guide_legend(ncol = 3, override.aes = list(size = 2))) +
    scale_y_sqrt() +
    # scale_y_break(c(21,30), scale='free') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 15),
        strip.background = element_blank(),
        panel.spacing.x = unit(0,"line"))
dev.off()
#################
library(UpSetR)
dat = as.data.frame.matrix(table(unique(res_sig[,c("LDblock", "ICDanno")])))
colnames(dat) = str_split(colnames(dat), pattern=" ", simplify=T)[,1]

pdf("res_sig_manhantan_ld_upsetR.pdf", width= 6, heigh = 6)
# upset(dat, sets = colnames(dat)[-1], order.by = "freq")
upset(dat, sets = names(sort(colSums(dat))), nintersects = 20,
    order.by = "freq", keep.order = T, point.size = 2, line.size = 0.8, mb.ratio = c(0.6, 0.4),
    mainbar.y.label = "Intersection Size", sets.x.label = "Set size")
aplot::plot_list(
    reshape2::melt(table(apply(dat, 1, function(x){paste(colnames(dat)[x>0], collapse=",")}))) %>% arrange(value) %>%
            mutate(Var1 = factor(Var1, levels = rev(Var1))) %>%
            ggplot() + geom_col(aes(Var1, value)) +
            geom_text(aes(Var1, value, label = value), hjust = 0, nudge_y = -0.2, size = 2.5) +
            labs(y = "Intersection Size", x = "") +
            scale_y_break(c(100, 260), scales = 0.3, ticklabels = c(270,320,350), expand = TRUE) +
            scale_y_sqrt(expand = c(0,0)) +
            theme_classic() +
            theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()),
    reshape2::melt(colSums(dat)) %>%
        tibble::rownames_to_column() %>% arrange(value) %>%
        mutate(rowname = factor(rowname, levels = rowname)) %>%
        ggplot() + geom_col(aes(value, rowname)) +
        geom_text(aes(value, rowname, label = value), hjust = 1, nudge_x = 0, size = 2.5) +
        labs(x = "Set Size", y = "") +
        scale_x_break(c(20, 90), scales = 0.2, ticklabels = 95, expand = TRUE) +
        scale_x_break(c(100, 450), scales = 0.3, ticklabels = c(460,800,1000), expand = TRUE) +
        scale_x_sqrt(expand = c(0,0)) +
        theme_classic() +
        theme(axis.ticks.y = element_blank()),
    nrow = 2, heights = c(0.6, 0.4))
dev.off()

############# manhattan
library(dplyr)

plot_all_icd_snp(res_sig, title = "res_sig_manhantan", ldblock = F)
plot_all_icd_snp(res_sig, title = "res_sig_manhantan_ld", ldblock = T)
plot_all_icd_snp(res_sig_raw, title = "res_sig_raw_manhantan", ldblock = F)
plot_all_icd_snp(res_sig_raw, title = "res_sig_raw_manhantan_ld", ldblock = T)

pdf("res_sig_raw_manhantan_ld_1e5.pdf", width=15, height=6)
plot_all_icd_snp(res_sig_raw[res_sig_raw$p.value<1e-5, ], title = "res_sig_raw_manhantan_ld_1e5", ldblock = T)
dev.off()

#test
plot_all_icd_snp = function(res_sig, title, ldblock = F){
    if(ldblock){
        res_sig = res_sig %>% group_by(ICDanno, LDblock) %>%
            summarise(CHR = CHR[which.min(p.value)], BP = BP[which.min(p.value)], OR = OR[which.min(p.value)], BETA = BETA[which.min(p.value)], P = min(p.value), cytoband = cytoband[which.min(p.value)]) %>% as.data.frame(stringsAsFactors = F)
    }else{
        res_sig = res_sig %>% mutate(P = p.value, OR_mean = OR)
    }
    res_sig$OR1 = res_sig$OR
    res_sig$OR1[res_sig$OR >= 15] = ">=15"
    res_sig$OR1[res_sig$OR >= 10 & res_sig$OR < 15] = "[10,15)"
    res_sig$OR1[res_sig$OR >= 5 & res_sig$OR < 10] = "[5,10)"
    res_sig$OR1[res_sig$OR >= 1 & res_sig$OR < 5] = "[1,5)"
    res_sig$OR1[res_sig$OR > 0 & res_sig$OR < 1] = "(0,1)"
    don <- res_sig %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len=max(as.numeric(BP))) %>% 
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        # Add this info to the initial dataset
        left_join(res_sig, ., by=c("CHR"="CHR")) %>%
        # Add a cumulative position of each SNP
        arrange(CHR, BP) %>%
        mutate(BPcum=BP+tot)
    axisdf = don %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2, minBPcum = min(BPcum))
    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
        "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
    don_sub = don # %>% group_by(icd) %>% mutate(ind = 1:length(LDblock)) #%>% filter(ind < 100)
    don_sub$ICDanno = factor(don_sub$ICDanno, levels = rev(names(sort(table(don_sub$ICDanno)))))
    don_sub1 = reshape2::melt(table(unique(don_sub[,c("ICDanno","LDblock")])[,1])) %>% mutate(percent = value/length(unique(don_sub$LDblock))) %>% arrange(desc(value))
    don_sub1$percent = paste0(signif(don_sub1$percent,3)*100, "%")
    don_sub1$Var1 = factor(don_sub1$Var1, levels = don_sub1$Var1)
    don_sub1$percent1 = don_sub1$percent
    p1 = ggplot(don_sub, aes(x=BPcum, y=-log10(P), color=ICDanno, size=OR1)) +
        geom_point(alpha=0.8, position=position_jitter(h = 0, w = 0.5)) +
        geom_hline(yintercept=7.3, linetype = "dashed") +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 5))) +
        scale_x_continuous(expand = c(0.01, 0.2), label = axisdf$CHR, breaks= axisdf$center, minor_breaks = axisdf$minBPcum) +
        scale_size_manual(breaks = c("(0,1)", "[1,5)", "[5,10)", "[10,15)", ">=15"), values = c(1, 1.5, 2, 2.5, 3)) +
        # scale_y_break(c(10, 10), scale='free', space = 0) +
        scale_y_cut(10) +
        theme_bw() +
        theme(
            axis.text = element_text(size=10),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            legend.position= 'bottom',
            legend.direction = "vertical",
            legend.box = "horizontal",
            panel.grid.major.x = element_blank())
    p2 = ggplot(don_sub1 %>% mutate(Var1 = factor(Var1, levels = rev(Var1))), aes(value, Var1, fill=Var1, label = percent)) + geom_col() +
        scale_fill_manual(values=rev(colors37[1:nrow(don_sub1)])) +
        # coord_fixed(ratio = 1/1000) +
        scale_x_sqrt(expand=c(0,0)) +
        geom_text(aes(x=max(don_sub1$value) * 0.75)) +
        theme_classic() +
        theme(axis.title=element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            legend.position= 'none') #axis.text.x = element_text(angle = 60, vjust=1,hjust=1),
    return(list(p1, p2))
    # pdf(paste0(title, ".pdf"), width = 15, heigh = 6)
    # plot(p1)
    # dev.off()
    # pdf(paste0(title, ".count.pdf"), width = 8, heigh = 3)
    # plot(p2)
    # dev.off()
}

plot_all_icd_snp = function(res_sig, title, ldblock = F){
    if(ldblock){
        res_sig = res_sig %>% group_by(ICDanno, LDblock) %>%
            summarise(CHR = CHR[which.min(p.value)], BP = BP[which.min(p.value)], OR = OR[which.min(p.value)], P = min(p.value), cytoband = cytoband[which.min(p.value)]) %>% as.data.frame(stringsAsFactors = F)
    }else{
        res_sig = res_sig %>% mutate(P = p.value)
    }
    don <- res_sig %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len=max(as.numeric(BP))) %>% 
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        # Add this info to the initial dataset
        left_join(res_sig, ., by=c("CHR"="CHR")) %>%
        # Add a cumulative position of each SNP
        arrange(CHR, BP) %>%
        mutate(BPcum=BP+tot)

    axisdf = don %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum))/2)

    colors37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda",
        "#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977",
        "#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

    don_sub = don # %>% group_by(icd) %>% mutate(ind = 1:length(LDblock)) #%>% filter(ind < 100)
    don_sub$ICDanno = factor(don_sub$ICDanno, levels = rev(names(sort(table(don_sub$ICDanno)))))
    don_sub1 = reshape2::melt(table(unique(don_sub[,c("ICDanno","LDblock")])[,1])) %>% mutate(percent = value/length(unique(don_sub$LDblock))) %>% arrange(desc(value))
    don_sub1$percent = paste0(signif(don_sub1$percent,3)*100, "%")
    don_sub1$Var1 = factor(don_sub1$Var1, levels = don_sub1$Var1)
    don_sub1$percent1 = don_sub1$percent
    # don_sub1$percent1[13:16] = ""
    # don_sub$icd = factor(don_sub$icd, levels = rev(names(sort(table(don_sub$icd)))))
    # don_sub$icd_group = substr(don_sub$icd, 1, 1)
    library(plotly)
    library(htmlwidgets)
    p1 = ggplot(don_sub, aes(x=BPcum, y=-log10(P), color=ICDanno, size=OR)) +
        # Show all points
        geom_point(alpha=0.8, position=position_jitter(h = 0, w = 0.5)) +
        # geom_point(alpha=0.8) +
        geom_hline(yintercept=7.3) +
        scale_color_manual(values=colors37, guide = guide_legend(ncol = 2, override.aes = list(size = 5))) +
        # custom X axis:
        scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
        scale_y_sqrt() + # expand = c(0,0)  remove space between plot area and x axis
        # Custom the theme:
        theme_bw() +
        theme(
            axis.text = element_text(size=10),
            axis.title.y = element_text(size=12),
            axis.title.x=element_blank(),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            legend.key.size = unit(1,"line"),
            legend.position= 'bottom',
            legend.direction = "vertical",
            legend.box = "horizontal")
    p2 = ggplot(don_sub1 %>% mutate(Var1 = factor(Var1, levels = rev(Var1))), aes(value, Var1, fill=Var1, label = percent)) + geom_col() +
        scale_fill_manual(values=rev(colors37[1:nrow(don_sub1)])) +
        # coord_fixed(ratio = 1/1000) +
        scale_x_sqrt(expand=c(0,0)) +
        geom_text(aes(x=max(don_sub1$value) * 0.75)) +
        theme_classic() +
        theme(axis.title=element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            legend.position= 'none') #axis.text.x = element_text(angle = 60, vjust=1,hjust=1),
    pdf(paste0(title, ".pdf"), width = 15, heigh = 4)
    plot(p1)
    dev.off()
    pdf(paste0(title, ".count.pdf"), width = 8, heigh = 3)
    plot(p2)
    dev.off()
}
```
