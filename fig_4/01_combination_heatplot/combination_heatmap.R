colorblind_palette = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(tidyverse)
library(ggplot2)
library(stringr)
library(hrbrthemes)
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
library(MetBrewer)

d2_threshold  = 0.639
d5_threshold  = 0.103

d2_baseline_threshold = -0.938
d5_baseline_threshold = -0.300

d5_color = "#28bedb"
d2_color = "#dbc60d"

df = read_csv("../../fig_2/01_activators_synergy/pairs_baselinesums.csv")
bases = read_csv("../../fig_1/03_computing_baselines/baseline_scores.csv")
oligos = read_csv("../../fig_1/01_raw_counts/csvs/base_oligo_library.csv")

oligos %>%
    dplyr::left_join(
        bases,
        by = c("label" = "domain"),
        suffix = c("_baseline", "_prior")
    ) -> singlets

df %>%
    dplyr::filter(df$composition %in% c("A-R", "R-A")) %>%
    dplyr::mutate(
        act_dom = if_else(composition == "A-R", d1_Gene, d2_Gene),
        oth_dom = if_else(composition == "A-R", d2_Gene, d1_Gene),
        ahit = avg_enrichment_d2 >= d2_threshold
    ) %>%
    dplyr::group_by(act_dom) %>%
    dplyr::summarise(num_hits = sum(ahit, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
        singlets,
        by = c("act_dom" = "Gene")
    ) -> act_hitted

p = ggplot(data = act_hitted, aes(x = med_d2, y = num_hits))
p = p + geom_smooth(method = "lm", color = "#e15759", size = 0.75, se = FALSE)
p = p + geom_point(size = 1, color = "#4e79a7")
p = p + geom_point(size = 0.005, color = "white")
p = p + geom_text(data = data.frame(x = 0.5, y = 7, label = bquote("Pearson R=0.81")),
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-1, 4), ylim = c(0, 7), ratio = 5/7)
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Act. Ctrl-Paired" ~ log[2] ~ "(ON:OFF)"))
p = p + ylab("Num. Act. Pairs w/ Reprs.")
ggsave("./activators_beating_repressors.pdf", p, height = 2.5, width = 2.5)


df %>%
    dplyr::filter(df$composition %in% c("A-R", "R-A")) %>%
    dplyr::mutate(
        rep_dom = if_else(composition == "A-R", d2_Gene, d1_Gene),
        oth_dom = if_else(composition == "A-R", d1_Gene, d2_Gene),
        rhit = avg_enrichment_d5 <= d5_threshold
    ) %>%
    dplyr::group_by(rep_dom) %>%
    dplyr::summarise(num_hits = sum(rhit, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(
        singlets,
        by = c("rep_dom" = "Gene")
    ) -> rep_hitted

p = ggplot(data = rep_hitted, aes(x = med_d5, y = num_hits))
p = p + geom_smooth(method = "lm", color = "#e15759", size = 0.75, se = FALSE)
p = p + geom_point(size = 1, color = "#4e79a7")
p = p + geom_point(size = 0.005, color = "white")
p = p + geom_text(data = data.frame(x = -2.25, y = 0, label = bquote("Pearson R=-0.45")),
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-3.5, 0), ylim = c(0, 40), ratio = 3.5/40)
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Rep. Ctrl-Paired" ~ log[2] ~ "(ON:OFF)"))
p = p + ylab("Num. Repr. Pairs w/ Acts.")
ggsave("./repressors_beating_activators.pdf", p, height = 2.5, width = 2.5)

df %>%
    dplyr::mutate(
        d1_combo_type = if_else(
            d1_baseline_type == "Non-hit",
            if_else(str_detect(d1_Description, "Act"), "Act", "Rep"),
            stringr::str_sub(d1_baseline_type, 1, 3)
        ),
        d2_combo_type = if_else(
            d2_baseline_type == "Non-hit",
            if_else(str_detect(d2_Description, "Act"), "Act", "Rep"),
            stringr::str_sub(d2_baseline_type, 1, 3)
        ),
        combo_type = paste(
            stringr::str_sub(d1_combo_type, 1, 1), 
            stringr::str_sub(d2_combo_type, 1, 1), 
            sep = "-")
    ) %>% 
    dplyr::filter(combo_type %in% c("A-R", "R-A")) %>%
    # dplyr::filter(df$composition %in% c("A-R", "R-A")) %>%
    dplyr::mutate(
        act_dom = if_else(d1_combo_type == "Act", d1_Gene, d2_Gene),
        act_dom = unlist(map(stringr::str_split(act_dom, " "), first)),
        rep_dom = if_else(d1_combo_type == "Act", d2_Gene, d1_Gene),
        rep_dom = unlist(map(stringr::str_split(rep_dom, " "), first)),
        act_val = if_else(d1_combo_type == "Act", d1_med_d2, d2_med_d2),
        rep_val = if_else(d1_combo_type == "Act", d2_med_d5, d1_med_d5),
        ar_act = avg_enrichment_d2,
        ar_rep = avg_enrichment_d5,
        pair_type = case_when(
            ar_act < d2_threshold & ar_rep > d5_threshold ~ "Neither",
            ar_act >= d2_threshold & ar_rep > d5_threshold ~ "Activator",
            ar_act < d2_threshold & ar_rep <= d5_threshold ~ "Repressor",
            ar_act >= d2_threshold & ar_rep <= d5_threshold ~ "Dual"
        )
    ) %>%
    dplyr::filter(pair_type %in% c("Neither", "Activator", "Repressor", "Dual")) -> comb

alist = unique(dplyr::arrange(comb, act_val)$act_dom)
rlist = rev(unique(dplyr::arrange(comb, rep_val)$rep_dom))

comb %>%
    dplyr::filter(
        (rep_dom %in% rlist) & (act_dom %in% alist)
    ) %>%
    dplyr::mutate(
        r_ind = match(rep_dom, rlist),
        a_ind = match(act_dom, alist)
    ) -> comb

comb$act_dom = factor(comb$act_dom, levels = alist)
comb$rep_dom = factor(comb$rep_dom, levels = rlist)

getmode = function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

comb %>%
    dplyr::group_by(act_dom, rep_dom) %>%
    dplyr::summarise(
        pair_type = getmode(pair_type),
        a_ind = a_ind[1],
        r_ind = r_ind[1]) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(a_ind), r_ind) %>%
    dplyr::select(act_dom, rep_dom, pair_type) %>%
    tidyr::pivot_wider(
        id_cols = act_dom,
        names_from = rep_dom,
        values_from = pair_type
    ) %>%
    tibble::column_to_rownames('act_dom') -> amat
amat[is.na(amat)] = NA_character_
amat = amat[, rlist]

get_activation = function(a) { 
    df %>% 
        dplyr::filter(str_detect(d1_Gene, a)) %>%
        dplyr::select(d1_med_d2) -> smola
    return(mean(smola$d1_med_d2, na.rm = TRUE))
}

get_repression = function(r) {
    df %>% 
        dplyr::filter(str_detect(d1_Gene, r)) %>%
        dplyr::select(d1_med_d5) -> smolr
    return(mean(smolr$d1_med_d5, na.rm = TRUE))
}

colors = structure(
    c('#23bedb', '#dbc60d', "#0fa616", '#999999'),
    names = c("Repressor", "Activator", "Dual", "Neither")
)

rep_col_fun = circlize::colorRamp2(c(-5, 1, 5),
                                   c(d5_color, "#dddddd", "#dbd7af"))
act_col_fun = circlize::colorRamp2(c(-5, -1, 5),
                                   c("#c2d3d6", "#dddddd", d2_color)) 

rep_annot = purrr::map_dbl(rlist, get_repression)
act_annot = purrr::map_dbl(alist, get_activation)

m_leg = Legend(
    title = "Combined Effect",
    labels = c("Repressor", "Activator", "Dual", "Neither"),
    legend_gp = gpar(fill = c(
        '#23bedb', '#dbc60d', "#0fa616", '#999999'
    )),
    labels_gp = gpar(fontsize = 14),
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    nr = 1
)
a_leg = Legend(
    title = "Med. Activation",
    col_fun = act_col_fun,
    direction = "horizontal",
    labels_gp = gpar(fontsize = 14),
    title_gp = gpar(fontsize = 14, fontface = "bold")
)
r_leg = Legend(
    title = "Med. Repression",
    col_fun = rep_col_fun,
    direction = "horizontal",
    labels_gp = gpar(fontsize = 14),
    title_gp = gpar(fontsize = 14, fontface = "bold")
)
pd = packLegend(a_leg,
                r_leg,
                m_leg,
                direction = "horizontal",
                column_gap = unit(1, "cm"))

padding = unit(3, "mm")
ht_opt$COLUMN_ANNO_PADDING = padding
ht_opt$ROW_ANNO_PADDING = padding

ht = Heatmap(
    amat,
    col = colors,
    rect_gp = gpar(col = "white", lwd = 0.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    na_col = "#FFFFFF",
    row_names_side = 'left',
    show_heatmap_legend = FALSE,
    row_title = "Activator %s",
    column_title = "Repressor %s",
    column_title_side = "bottom",
    row_split = factor(rev(c(rep("Non-Hit", 12), rep("Hit", 19))), levels = c("Hit", "Non-Hit")),
    row_order = rev(alist),
    row_gap = unit(2, "mm"),
    column_split = factor(c(rep("Non-Hit", 13), rep("Hit", 28)), levels = c("Non-Hit", "Hit")),
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    column_order = rlist,
    column_gap = unit(2, "mm"),
    bottom_annotation = HeatmapAnnotation(
        avg_rep = rep_annot,
        col = c(avg_rep = rep_col_fun),
        na_col = "#ffffff",
        annotation_label = "Med. Repr. →",
        annotation_name_side = "left",
        show_legend = FALSE,
        gp = gpar(col = "white", lwd = 0.5),
        gap = unit(9, "mm")
    ),
    left_annotation = rowAnnotation(
        avg_act = rev(act_annot),
        col = c(avg_act = act_col_fun),
        na_col = "#ffffff",
        annotation_label = "Med\nAct.\n↓",
        annotation_name_side = "top",
        annotation_name_rot = 0,
        show_legend = FALSE,
        gp = gpar(col = "white", lwd = 0.5),
        gap = unit(9, "mm")
    ),
)
cairo_pdf(
    "./activator_repressor_combos.pdf",
    width = 10,
    height = 8.5,
    family = "Arial")
draw(ht, annotation_legend_list = pd, annotation_legend_side = "top")
dev.off()

ht_opt(RESET = TRUE)

