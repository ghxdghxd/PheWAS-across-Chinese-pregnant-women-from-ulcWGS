# Figure S5

```sh
################################ LD-R2 of imputed variants
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

index = sample(1:nrow(dat), 10000)
index_1 = sample(which(abs(dat$R2_diff) < 0.01), 10000)

pdf("FigureS5.pdf")
ggplot(dat[index_1, ], aes(R2, R2_1kg)) +
        geom_density_2d_filled(contour_var = "ndensity") +
        annotate("text", x = 0.3, y = 0.9, color = "white", 
            label = paste("rho =", signif(cor(dat$R2[which(abs(dat$R2_diff) < 0.01)], dat$R2_1kg[which(abs(dat$R2_diff) < 0.01)]), 3))) +
        labs(x="LD R2 minimac4", y="LD R2 in 1KG") +
        scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
        theme_bw() +
        theme(aspect.ratio = 1,
            legend.position = c(0.8, 0.3),
            legend.position = "none")
dev.off()
```
