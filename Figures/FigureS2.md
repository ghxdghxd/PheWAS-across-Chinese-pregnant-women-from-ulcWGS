# Figure S2

[Ancestry analysis](../workflows/04_ancestryAnalysis.md)

```R
load("Pop_pca_NIPTadd1kgCHB_JPT_com0.1_callRate0.99_conf95_SYNCSA.RData")

pdf("FigureS2.pdf", width = 5, height = 5)
plot(res, show = "variables", arrows = TRUE)
plot(res, show = "individuals", axis = c(1, 2), text = TRUE)
plot(res, show = "individuals", text = FALSE, points = TRUE)
ggplot() + geom_point(data=a, aes(Axis.1, Axis.2, color = id)) +
    theme_bw() + theme(aspect.ratio=1)
ggplot() + geom_point(data=a[a$id!="Case",], aes(Axis.1, Axis.2, color = id1)) +
    geom_point(data=a[a$id=="Case",], aes(Axis.1, Axis.2), alpha = 50, shape = 1) +
    theme_bw() + theme(aspect.ratio=1)
ggplot() + geom_point(data=a[a$id=="Case",], aes(Axis.1, Axis.2), alpha = 50, shape = 1) +
    geom_point(data=a[a$id!="Case",], aes(Axis.1, Axis.2, color = id1)) +
    theme_bw() + theme(aspect.ratio=1)
ggplot() + geom_point(data=a[a$id=="Case",], aes(Axis.1, Axis.2), alpha = 50, shape = 1) +
    geom_point(data=a[a$id!="Case",], aes(Axis.1, Axis.2, color = id2)) +
    theme_bw() + theme(aspect.ratio=1)
ggplot() + geom_point(data=a[a$id=="Case",], aes(Axis.1, Axis.2)) +
    theme_bw() + theme(aspect.ratio=1)
dev.off()
```
