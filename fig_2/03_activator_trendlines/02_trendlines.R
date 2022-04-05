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

# trendline... trends
df = read_csv('./activators_fitted.csv')

# activation vs Kd
p = ggplot()
p = p + geom_smooth(
    data = df,
    aes(
        x = med_d2,
        y = Kd
    ),
    size = 0.75,
    se = FALSE,
    n = 100,
    color = "#e15759"
)
p = p + geom_point(
    data = df,
    aes(
        x = med_d2,
        y = Kd
    ),
    size = 1,
    color = "#4e79a7"
)
p = p + geom_point(
    data = df,
    aes(
        x = med_d2,
        y = Kd
    ),
    size = 0.005,
    color = "white"
)
p = p + geom_text(data = data.frame(x = 0.6, y = -0.05, label = "Pearson R = -0.96"),
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-1, 3.3), ylim=c(-2, 0), ratio = 4.3/2)
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Median Act." ~ log[2](ON:OFF)))
p = p + ylab(bquote("Partner Strength at\nHalf-Maximal Pair Activation"))
ggsave("./kd_vs_med_d2.pdf", p, height = 2.5, width = 2.6)

# activation vs. hill coeff.
p = ggplot(data = df, aes(x = med_d2, y = n))
p = p + geom_smooth(
    method = "lm",
    size = 0.75,
    se = FALSE,
    n = 100,
    color = "#e15759"
)
p = p + geom_point(
    size = 1,
    color = "#4e79a7"
)
p = p + geom_point(
    size = 0.005,
    color = "white"
)
p = p + geom_text(data = data.frame(x = 0.35, y = 24, label = "Pearson R = 0.72"),
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-1, 3.3), ylim=c(0, 24.5), ratio = 4.3/25)
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Median Act." ~ log[2](ON:OFF)))
p = p + ylab("Hill Coefficient")
ggsave("./hn_vs_med_d2.pdf", p, height = 2.5, width = 2.5)

# activation vs synergy
df = read_csv("../../fig_2/01_activators_synergy/pairs_baselinesums.csv")
df %>%
    dplyr::filter(composition == "A-A") %>%
    dplyr::group_by(domain1) %>%
    dplyr::summarise(
        mean_syn = mean((avg_enrichment_d2 - d1_med_d2 - d2_med_d2), na.rm = TRUE),
        d1_med_d2 = mean(d1_med_d2, na.rm = TRUE),
        mdiv_syn = mean_syn/d1_med_d2,
    ) -> syndf
p = ggplot(syndf, mapping = aes(x = d1_med_d2, y = mean_syn))
p = p + geom_smooth(method = "lm", size = 0.75, se = FALSE, n = 100, color = "#e15759")
p = p + geom_point(size = 1, color = "#4e79a7")
p = p + geom_point(size = 0.005, color = "white")
p = p + geom_text(data = data.frame(x = 1.9, y = 2.9, label = bquote("Pearson R=-0.89")), 
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE) 
p = p + coord_fixed(xlim = c(-1, 3.3), ylim = c(-2, 3), ratio = 4.3/5)
p = p + xlab("Median Act." ~ log[2](ON:OFF)) + ylab("Avg. Synergy with Activators")
p = p + theme_linedraw() 
p = p + theme(panel.grid = element_blank())
ggsave("./syn_vs_med_d2.pdf", p, height = 2.5, width = 2.5)
