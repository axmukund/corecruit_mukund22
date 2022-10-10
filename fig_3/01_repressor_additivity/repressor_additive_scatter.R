# imports -----

colorblind_palette = c(
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
)
library(tidyverse)
library(ggplot2)
library(stringr)
library(hrbrthemes)
library(tidyverse)
library(cowplot)
library(deSolve)
library(data.table)
library(patchwork)
library(comprehenr)
library(tidyr)
library(ggthemes)
library(ggeasy)
library(dplyr)
library(viridis)
library(gtable)
library(corrr)
library(ggsignif)
library(ggforce)
library(rlist)
library(ggridges)
library(ggHoriPlot)
library(latex2exp)
library(ggrepel)
library(interp)
library(ggpubr)
library(scico)
library(fitdistrplus)
library(colorspace)
library(biclust)
library(heatmaply)
library(webshot)
library(dendextend)
library(gplots)
library(ComplexHeatmap)
library(circlize)

d2_threshold  =  -0.0236
d5_threshold  = 0.771

# d2_baseline_threshold = -0.938
# d5_baseline_threshold = -0.300
d5_color = "#28bedb"
d2_color = "#dbc60d"

df = read_csv('../../fig_1/04_scatter/pairs_baselined.csv')
df %>%
    dplyr::mutate(baseline_sum_d2 = d1_med_d2 + d2_med_d2,
                  baseline_sum_d5 = d1_med_d5 + d2_med_d5) -> df
write_csv(df, "./pairs_baselinesums.csv")

# plots -----

df %>%
    dplyr::mutate(
        makeup = case_when(
            composition == "C-C" ~ "Ctrl + Ctrl",
            # (composition == "C-R" | composition == "R-C") ~ "Ctrl + Rep",
            composition == "R-R" ~ "Rep + Rep",
            composition %in% c("D-R", "R-D") ~ "Rep + Dual",
            composition == "D-D" ~ "Dual + Dual",
            TRUE ~ "N/A"
        )
    ) %>%
    dplyr::filter(composition %in% c("C-C", "R-R", "R-D", "D-R", "D-D")) %>%
    dplyr::mutate(pname = paste(d1_Gene, d2_Gene, sep="-")) -> rdf
rdf$makeup = factor(rdf$makeup,
                    levels = c("Ctrl + Ctrl", "Dual + Dual", "Rep + Dual", "Rep + Rep"))

p1 = ggplot(data = dplyr::arrange(rdf, desc(makeup)))
p1 = p1 + geom_hline(yintercept = d5_threshold,
                     color = "#bababa",
                     linetype = "solid")
p1 = p1 + geom_vline(xintercept = d5_threshold,
                     color = "#bababa",
                     linetype = "solid")
p1 = p1  + geom_point(aes(x = baseline_sum_d5, y = avg_enrichment_d5, fill = makeup),
                      size = 1.0, stroke = 0.15, shape=21, color = "#bababa")
p1 = p1 + geom_abline(
    slope = 1,
    intercept = 0,
    color = "#777777",
    linetype = "dashed"
)
p1 = p1 + scale_fill_manual(
    breaks = c("Ctrl + Ctrl", "Dual + Dual", "Rep + Dual", "Rep + Rep"),
    values = c("#777777", "#be6eff", "#7196ed", "#23bedb"),
    guide = guide_legend(override.aes = list(size = 3)),
    name = ""
)
p1 = p1 + geom_label_repel(
    data = dplyr::filter(rdf, pname %in% c("FOXO3-ZNF10", "ZNF10-CBX1", "BIN1-FOXO3")),
    aes(
        x = baseline_sum_d5,
        y = avg_enrichment_d5,
        label = pname,
    ),
    size = 3,
    ylim = c(-8, -5),
    force = 10,
    box.padding = 0.1,
    max.time = 5,
    max.iter = 100000
)
p1 = p1 + coord_fixed(ratio = 10 / 12,
                      xlim = c(-7.5, 5),
                      ylim = c(-7.5, 5))
p1 = p1 + labs(x = "Sum of Control-Paired log2(ON:OFF)", y = "Concatenation log2(ON:OFF)")
p1 = p1 + theme_bw()
p1 = p1 + theme(
    legend.title = element_blank(),
    legend.position = c(0.18, 0.86),
    legend.margin = margin(-10, 0, 0, 0),
    panel.grid = element_blank(),
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
)
ggsave("./rep_scatter.pdf", p1, height = 2.75, width = 3.1)
p1
