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

d2_threshold = -0.0236
d5_threshold = 0.771

d5_color = "#28bedb"
d2_color = "#dbc60d"

concats = readr::read_csv("../../fig_2/01_activators_synergy/pairs_baselinesums.csv")

# avr scatter -----
concats %>%
    #    dplyr::filter(composition %in% c("A-R", "R-A", "A-N", "N-A", "R-N", "N-A", "C-C")) %>%
    dplyr::filter(composition %in% c(
        "A-R",
        "R-A",
        "A-A",
        "R-R",
        "C-C",
        "D-D",
        "A-D",
        "D-A",
        "A-R",
        "R-A"
    )) %>%
    dplyr::mutate(
        ctype = case_when(
            composition %in% c("A-R", "R-A", "D-A", "D-R", "A-D", "R-D", "D-D") ~ "Both",
            composition %in% c("A-N", "N-A", "A-A") ~ "Act",
            composition %in% c("R-N", "N-A", "R-R") ~ "Rep",
            composition == "C-C" ~ "Control"
        )
    ) -> cdf
cdf$ctype = factor(cdf$ctype, levels = c("Rep", "Both", "Control", "Act"))
cdf = dplyr::arrange(cdf, ctype)

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
    breaks = c("Act", "Both", "Rep", "Control"),
    labels = c("Act", "Both", "Rep", "Control"),
    values = c(d2_color, "#0fa616", d5_color, "#777777"))
    # values = c("#0fa616", d2_color, d5_color, '#777777'))
p = p + coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-7, 4.5))
p = p + guides(color = guide_legend(override.aes = list(size = 3)))
p = p + labs(x = bquote("Activation" ~ log[2] ~ "(ON:OFF)"),
             y = bquote("Repression" ~ log[2] ~ "(ON:OFF)"))
p = p + theme_bw()
p = p + theme(legend.position = c(0.8, 0.2), legend.title = element_blank(),
              panel.grid = element_blank(), legend.margin = margin(unit(c(0, 0, 0, 0), units = "mm"),),
              legend.key.size = unit(c(4), units = "mm"))
ggsave("./d2_vs_d5_scatter.pdf", p, height = 2.5, width = 2.5)
p

# scatters ----------

cdf %>%
    dplyr::filter(ctype == "Both") %>%
    # dplyr::mutate(ar_name = case_when(
    #     composition %in% c("D-D", "A-R", "A-D") ~ paste(d1_Gene, d2_Gene),
    #     composition %in% c("R-A", "D-A") ~ paste(d2_Gene, d1_Gene)
    # )) %>%
    dplyr::mutate(ar_name = if_else(
        composition %in% c("D-D", "A-R", "A-D"),
        paste(d1_Gene, d2_Gene),
        paste(d2_Gene, d1_Gene)
    )) %>%
    tidyr::pivot_wider(
        id_cols = c(ar_name, ctype),
        names_from = composition,
        values_from = c(avg_enrichment_d2, avg_enrichment_d5),
    ) %>%
    dplyr::rename(
        d2_ar = `avg_enrichment_d2_A-R`,
        d2_ra = `avg_enrichment_d2_R-A`,
        d5_ar = `avg_enrichment_d5_A-R`,
        d5_ra = `avg_enrichment_d5_R-A`,
    ) -> cross_df


aplot = ggplot(data = cross_df, aes(x = d2_ar, y = d2_ra)) +
    geom_smooth(
        formula = y ~ x,
        method = "lm",
        size = 0.75,
        se = FALSE,
        n = 100,
        color = "#e15759"
    ) +
    # geom_point(size = 0.002, color = "#4e79a7") +
    geom_point(size = 0.005, color = "#4e79a7") +
    geom_text(
        data = data.frame(
            x = -3.5,
            y = 4.25,
            label = bquote("Pearson\nR=0.82")
        ),
        mapping = aes(x = x, y = y, label = label),
        label.size = NA,
        inherit.aes = FALSE
    ) +
    coord_fixed(xlim = c(-5, 5), ylim = c(-5, 5)) +
    geom_vline(aes(xintercept = d2_threshold), linetype = 'dashed') +
    geom_hline(aes(yintercept = d2_threshold), linetype = 'dashed') +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(title = "Activation", x = 'Activator-Repressor', y = 'Repressor-Activator')
ggsave("./d2_ar_vs_ra_scatter.pdf", aplot, height = 2.5, width=2.5)

rplot = ggplot(data = cross_df, aes(x = d5_ar, y = d5_ra)) +
    geom_smooth(
        formula = y ~ x,
        method = "lm",
        size = 0.75,
        se = FALSE,
        n = 100,
        color = "#e15759"
    ) +
    # geom_point(size = 0.002, color = "#4e79a7") +
    geom_point(size = 0.005, color = "#4e79a7") +
    coord_fixed(ratio = 8/9.5, xlim = c(-4, 4), ylim = c(-5.5, 4)) +
    geom_vline(aes(xintercept = d5_threshold), linetype = 'dashed') +
    geom_hline(aes(yintercept = d5_threshold), linetype = 'dashed') +
    geom_label(
        data = data.frame(
            x = 1.5,
            y = -5.25,
            label = bquote("Pearson R=0.80")
        ),
        mapping = aes(x = x, y = y, label = label),
        label.size = NA,
        inherit.aes = FALSE
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(title = "Repression", x = 'Activator-Repressor', y = 'Repressor-Activator')
ggsave("./d5_ar_vs_ra_scatter.pdf", rplot, height = 2.5, width=2.5)


aplot + rplot
