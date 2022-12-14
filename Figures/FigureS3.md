# Figure S3

```R
load("SNPs_AF_callNum_bed.RData")

library(ggplot2)
library(RColorBrewer)

p1 = ggplot(bed[bed$type!="E", ], aes(obAFmean)) + geom_histogram(binwidth=0.01) +
    scale_fill_distiller(palette=4, direction=1) +
    scale_x_continuous(expand = c(0.01, 0), limits = c(0,1)) +
    theme_bw() +
    theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
p2 = ggplot(bed[bed$type!="E", ], aes(AF)) + geom_histogram(binwidth=0.01) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    scale_x_continuous(expand = c(0.01, 0), limits = c(0,1)) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
p3 = ggplot(bed[bed$type!="E", ], aes(obAFmean, AF)) +
    stat_density_2d(aes(fill = sqrt(..density..)), geom = "raster", contour = F) +
    annotate("text", x = 0.3, y = 0.9, label = paste("rho =", signif(as.numeric(cor.test(bed$AF[bed$type!="E"], bed$obAFmean[bed$type!="E"])$estimate),3))) +
    labs(x="NIPT AF", y = "AF in ChinaMap") +
    scale_fill_gradient(low = "white", high = "cornflowerblue") +
    scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    theme_bw() +
    theme(aspect.ratio = 1,
        legend.position = "none")

p31 = ggplot(bed[bed$type!="E", ], aes(obAFmean, AF)) +
    stat_density_2d(aes(fill = sqrt(..density..)), geom = "raster", contour = F) +
    annotate("text", x = 0.3, y = 0.9, label = paste("rho =", signif(as.numeric(cor.test(bed$AF[bed$type!="E"], bed$obAFmean[bed$type!="E"])$estimate),3))) +
    labs(x="NIPT AF", y = "AF in ChinaMap") +
    # scale_fill_discrete() +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    theme_bw() +
    theme(aspect.ratio = 1,
        legend.position = "none")

pdf("FigureS3.pdf", width=6, height=6)
plot_grid(p1, NULL, p3, p2, nrow = 2, align = 'hv', axis = "tblr", rel_widths = c(1, 0.5), rel_heights = c(0.5,1))
dev.off()
```
