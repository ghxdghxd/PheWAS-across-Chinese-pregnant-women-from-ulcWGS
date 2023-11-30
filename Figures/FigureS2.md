# Figure S2

[Ancestry analysis](../workflows/04_ancestryAnalysis.md)

```R
load("Pop_pca_1kg_com0.1_callRate0.99_conf95_SYNCSA.RData")
a = as.data.frame(res$individuals[,1:10], stringsAsFactors = F)
a$ID = rownames(a)
a$id = ID_1000G$Superpopulation.code[match(a$ID, ID_1000G$Sample.name)]
a$id1 = ID_1000G$Population.code[match(a$ID, ID_1000G$Sample.name)]
a$id2 = a$id1
a$id2[grep("CH",a$id2)] = "HAN"

pb = ggplot() + geom_point(data=a, aes(Axis.1, Axis.2, color = id)) +
    theme_bw() + theme(aspect.ratio=1)


eig = res$eigenvalues
eig$PC = 1:nrow(eig)
eig$value = eig$Perc.Inertia * 100

pa = ggplot(head(eig, 50), aes(PC, value)) + geom_point() + geom_line() + theme_bw() +
    labs(x="Dimensions", y="Percentage of explained variances")


load("Pop_pca_NIPTadd1kgCHB_JPT_com0.1_callRate0.99_conf95_SYNCSA.RData")

a = as.data.frame(res$individuals[,1:10], stringsAsFactors = F)
a$ID = c(sample_ID[!duplicated(sample_ID)], ID_1000G$Sample.name[ID_1000G$Population.code %in% c("CHB","CHS","CDX", "JPT")])
a$id = c(rep("Case", 48734), ID_1000G$Superpopulation.code[match(ID_1000G$Sample.name[ID_1000G$Population.code %in% c("CHB","CHS","CDX", "JPT")], ID_1000G$Sample.name)])
a$id1 = c(rep("Case", 48734), ID_1000G$Population.code[match(ID_1000G$Sample.name[ID_1000G$Population.code %in% c("CHB","CHS","CDX", "JPT")], ID_1000G$Sample.name)])
a$id2 = a$id1
a$id2[grep("CH",a$id2)] = "HAN"

eig = res$eigenvalues
eig$PC = 1:nrow(eig)
eig$value = eig$Perc.Inertia * 100


pd = ggplot() + geom_point(data=a, aes(Axis.1, Axis.2, color = id)) +
    theme_bw() + theme(aspect.ratio=1)


library(factoextra)
eig = res$eigenvalues
eig$PC = 1:nrow(eig)
eig$value = eig$Perc.Inertia * 100

pc = ggplot(head(eig, 50), aes(PC, value)) + geom_point() + geom_line() + theme_bw() +
    labs(x="Dimensions", y="Percentage of explained variances")


pdf("FigureS2.pdf", widht = 4, height = 4)
pa
pb
pc
pd
dev.off()
```
