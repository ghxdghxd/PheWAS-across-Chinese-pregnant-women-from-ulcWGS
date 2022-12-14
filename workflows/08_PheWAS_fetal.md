# 8. PheWAS for newborns

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

## 8.2 PheWAS for fetal variants

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

## 8.3 Maternal Pleiotropic Loci associate with neonatal phenotypes

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

library(data.table)
library(stringr)
load("phewas_power_P0.001.RData")
res_sig_raw = snp_power[which(snp_power$power > 0.8), ]

ICD_res = ICD_res_filter[which(ICD_res_filter$snp %in% res_sig_raw$MarkerID),]
save(ICD_res, file="ICD_res.RData")

ICD_res_sig = ICD_res[ICD_res$pvalue_norm < 0.05, ]
write.table(ICD_res_sig, file="ICD_res_sig.txt", r=F, c=T, sep="\t", quote=F)
```
