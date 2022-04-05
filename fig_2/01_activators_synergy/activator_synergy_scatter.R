colorblind_palette = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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

d2_threshold  = 0.639
d5_threshold  = 0.103

d2_baseline_threshold = -0.938
d5_baseline_threshold = -0.300
d5_color = "#28bedb"
d2_color = "#dbc60d"

df = read_csv('../../fig_1/04_scatter/pairs_baselined.csv')
df %>%
    dplyr::mutate(
        baseline_sum_d2 = d1_med_d2 + d2_med_d2,
        baseline_sum_d5 = d1_med_d5 + d2_med_d5) -> df
write_csv(df, "./pairs_baselinesums.csv")

df %>%
    dplyr::mutate(
        makeup = case_when(
            composition == "C-C" ~ "Ctrl + Ctrl",
            (composition == "C-A" | composition == "A-C") ~ "Ctrl + Act",
            composition == "A-A" ~ "Act + Act",
            TRUE ~ "N/A"
        )
    ) %>%
    dplyr::filter(composition %in% c("C-C", "A-A")) -> adf


p1 = ggplot(data = adf)
p1 = p1 + geom_hline(yintercept = d2_threshold,
                     color = "#bababa",
                     linetype = "solid")
p1 = p1 + geom_vline(xintercept = d2_threshold,
                     color = "#bababa",
                     linetype = "solid")
p1 = p1  + geom_point(aes(x = baseline_sum_d2, y = avg_enrichment_d2, color = makeup),
                      size = 0.75)
p1 = p1 + geom_abline(
    slope = 1,
    intercept = 0,
    color = "#777777",
    linetype = "dashed"
)
p1 = p1 + geom_text_repel(
    data = dplyr::filter(adf,
                         paste(d1_Gene, d2_Gene, sep =
                                   '-') %in% c("ANM2-KIBRA",
                                               "NOTC2-ANM2",
                                               "NOTC2-KIBRA")),
    aes(
        label = paste(d1_Gene, d2_Gene, sep = '-'),
        x = baseline_sum_d2,
        y = avg_enrichment_d2
    ),
    min.segment.length = 0,
    nudge_x = c(-5, 3, 3),
    nudge_y = c(0.5, 5, -2)
)
p1 = p1 + scale_color_manual(
    breaks = c("Ctrl + Ctrl", "Ctrl + Act", "Act + Act"),
    values = c("#777777", "#dbc60d", "#ff8d6e"),
    guide = guide_legend(override.aes = list(size = 3)),
    name = ""
)
p1 = p1 + coord_fixed(ratio = 10 / 12,
                      xlim = c(-4, 7),
                      ylim = c(-4.5, 6.5))
p1 = p1 + labs(x = "Sum of Control-Paired log2(ON:OFF)", y = "Concatenation log2(ON:OFF)")
p1 = p1 + theme_bw()
p1 = p1 + theme(
    legend.position = c(0.8, 0.2),
    legend.margin = margin(-10, 0, 0, 0),
    panel.grid = element_blank()
)
ggsave("./act_scatter.pdf", p1, height = 2.75, width = 3.1)