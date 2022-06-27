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

d2_threshold  = 0.639
d5_threshold  = 0.103

d2_baseline_threshold = -0.938
d5_baseline_threshold = -0.300

d5_color = "#28bedb"
d2_color = "#dbc60d"

concats = readr::read_csv("../../fig_2/01_activators_synergy/pairs_baselinesums.csv")

concats %>%
    dplyr::filter(composition %in% c("A-R", "R-A", "A-N", "N-A", "R-N", "N-A", "C-C")) %>%
    dplyr::mutate(ctype = case_when(
        composition %in% c("A-R", "R-A") ~ "Both",
        composition %in% c("A-N", "N-A") ~ "Act",
        composition %in% c("R-N", "N-A") ~ "Rep",
        composition == "C-C" ~ "Control"
    )) -> cdf

p = ggplot()
p = p + geom_point(
    data = cdf,
    aes(
        x = avg_enrichment_d2,
        y = avg_enrichment_d5,
        color = ctype,
    ),
    size = 0.15
)
p = p + geom_hline(yintercept = d5_threshold, color = 'black', linetype = 'dashed') 
p = p + geom_vline(xintercept = d2_threshold, color = 'black', linetype = 'dashed')
p = p + scale_color_manual(
    name = "",
    breaks = c("Both", "Act", "Rep", "Control"),
    labels = c("A-R",
               "A-N/A",
               "R-N/A",
               "Control"),
    values = c("#0fa616", d2_color, d5_color, '#777777'))
p = p + coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-7, 4.5))
p = p + guides(color = guide_legend(override.aes = list(size = 3)))
p = p + labs(x = bquote("Activation" ~ log[2] ~ "(ON:OFF)"),
             y = bquote("Repression" ~ log[2] ~ "(ON:OFF)"))
p = p + theme_bw() 
p = p + theme(legend.position = c(0.8, 0.2), legend.title = element_blank(), 
              panel.grid = element_blank(), legend.margin = margin(unit(c(0, 0, 0, 0), units = "mm"),),
              legend.key.size = unit(c(4), units = "mm"))
ggsave("./d2_vs_d5_scatter.pdf", p, height = 2.5, width = 2.5)
