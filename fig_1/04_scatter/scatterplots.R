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
library(ggpubr)

d2_color = "#28bedb"
d5_color = "#dbc60d"
d10_color = "#db9d0d"

df = readr::read_csv('../02_filtering_controls/02_pairs_priored_controlfiltered.csv')
baseline_df = readr::read_csv('../03_computing_baselines/baseline_scores.csv')

b1df = copy(baseline_df)
b2df = copy(baseline_df)
names(b1df) = paste0("d1_", names(b1df))
names(b2df) = paste0("d2_", names(b2df))

df %>% 
    dplyr::select(names(df)[1:36]) -> df

df %>%
    dplyr::left_join(
        b1df,
        by = c("domain1" = "d1_domain")
    ) %>%
    dplyr::left_join(
        b2df,
        by = c("domain2" = "d2_domain")
    ) -> df

abbrvchar = function(d1, d2) { return(paste(substr(d1, 1, 1), substr(d2, 1, 1), sep = "-")) }

df %>% # add compositions
    dplyr::mutate(
        composition = abbrvchar(d1_baseline_type, d2_baseline_type),
    )  %>%
    dplyr::mutate(
        character = case_when(
            str_detect(composition, "C-C") ~ "Control",
            str_detect(composition, "C-R") ~ "Repressor",
            str_detect(composition, "C-A") ~ "Activator",
            str_detect(composition, "R-C") ~ "Repressor",
            str_detect(composition, "R-R") ~ "Repressor",
            str_detect(composition, "R-A") ~ "Both",
            str_detect(composition, "A-C") ~ "Activator",
            str_detect(composition, "A-R") ~ "Both",
            str_detect(composition, "A-A") ~ "Activator",
            TRUE ~ "Other"
        )
    ) -> df
df$character = factor(df$character, levels = c("Control", "Activator",
                                               "Repressor", "Both", "Other"))


df %>%
    dplyr::filter(character == "Control") %>%
    tidyr::pivot_longer(
        cols = c(avg_enrichment_d2, avg_enrichment_d5), 
        names_to = "day",
        names_prefix = "avg_enrichment_d",
        values_to = "avg_enrichment"
    ) %>%
    dplyr::select(pair, day, avg_enrichment) %>%
    tidyr::drop_na() %>%
    dplyr::group_by(day) %>%
    summarise(
        mean_enrichment = mean(avg_enrichment),
        sd_enrichment   = sd(avg_enrichment)
    ) %>%
    dplyr::mutate(
        threshold = if_else(day==2, mean_enrichment + 3 * sd_enrichment,
                            mean_enrichment - 3 * sd_enrichment)) -> control_df
df %>%
    dplyr::mutate(
        act_pair_hit = avg_enrichment_d2 >= 0.639,
        rep_pair_hit = avg_enrichment_d5 <= 0.103
    ) -> df

# Day 2 Scatterplot
corall = cor(df$enrichment_ratio_r1_d2, df$enrichment_ratio_r2_d2, use = "complete.obs")
p2 = ggplot(data = df)
p2 = p2 + geom_point(aes(x = enrichment_ratio_r1_d2, y = enrichment_ratio_r2_d2, color = character), size = 0.15)
p2 = p2 + geom_segment(
    x = 3.5, y = dplyr::filter(control_df, day==2)$threshold - 3.5,
    xend = dplyr::filter(control_df, day==2)$threshold - 3.5, yend = 3.5,
    linetype = "dashed", size = 0.5
)
p2 = p2 + annotate("text", x = -4.5, y = 6.5, size = 9 / .pt, family = "Helvetica",
                   label = bquote("Pearson" ~ rho ~ "=" ~ .(round(corall, digits= 2))))# ~ "for all domains"))
p2 = p2 + coord_fixed(ratio = 1, xlim = c(-7.5, 6.5), ylim=c(-7.5, 6.5))
p2 = p2 + theme_linedraw() + theme(panel.grid = element_blank())
p2 = p2 + guides(color = guide_legend(override.aes = list(size = 3)))
p2 = p2 + labs(x = bquote("Activation" ~ log[2] ~ "(ON:OFF)" ~ "R1"),
               y = bquote("Activation" ~ log[2] ~ "(ON:OFF)" ~ "R2"))
p2 = p2 + scale_color_manual(
    name = "", 
    breaks = c("Control", "Repressor", "Activator", "Both", "Other"),
    labels = c("C + C",
               "R + C/R", 
               "A + C/R",
               "R + A",
               "Other"),
    values = c("#777777", '#23bedb', '#dbc60d', "#0fa616", '#dedede'))

# Day 5 Scatterplot
corall = cor(df$enrichment_ratio_r1_d5, df$enrichment_ratio_r2_d5, use = "complete.obs")
p5 = ggplot(data = df)
p5 = p5 + geom_point(aes(x = enrichment_ratio_r1_d5, y = enrichment_ratio_r2_d5, color = character), size = 0.15)
p5 = p5 + geom_segment(
    x = 4, y = dplyr::filter(control_df, day==5)$threshold - 4,
    xend = dplyr::filter(control_df, day==5)$threshold - 4, yend = 4,
    linetype = "dashed", size = 0.5
)
# p5 = p5 + geom_abline(slope = 1, intercept = 0, size = 0.5, color = "#444444")
p5 = p5 + annotate("text", x = -4, y = 7, size = 9 / .pt, family = "Helvetica",
                   label = bquote("Pearson" ~ rho ~ "=" ~ .(round(corall, digits= 2))))# ~ "for all domains"))
p5 = p5 + coord_fixed(ratio = 1, xlim = c(-7, 7), ylim = c(-7, 7))
p5 = p5 + theme_linedraw() + theme(panel.grid = element_blank())
p5 = p5 + scale_color_manual(
    name = "", 
    breaks = c("Control", "Repressor", "Activator", "Both", "Other"),
    labels = c("C + C",
               "R + C/R", 
               "A + C/A",
               "R + A",
               "Other"),
    values = c("#777777", '#23bedb', '#dbc60d', "#0fa616", '#dedede'))
p5 = p5 + guides(color = guide_legend(override.aes = list(size = 3)))
p5 = p5 + labs(x = bquote("Repression" ~ log[2] ~ "(ON:OFF)" ~ "R1"),
               y = bquote("Repression" ~ log[2] ~ "(ON:OFF)" ~ "R2"))

# Put the plots together
p2 = p2 + theme(legend.text = element_text(size = 8))
p5 = p5 + theme(legend.text = element_text(size = 8))
p = ggarrange(p2, p5, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggsave("./d2_d5_scatter.pdf", p, height = 3, width = 6, useDingbats=FALSE)

write_csv(df, './pairs_baselined.csv')