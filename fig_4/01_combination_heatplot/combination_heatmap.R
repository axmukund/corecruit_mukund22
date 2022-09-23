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

# d2_threshold  = 0.639
# d5_threshold  = 0.103
d2_threshold = -0.0236
d5_threshold = 0.771

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

# activators beating repressors -----
df %>%
    dplyr::filter(composition %in% c("A-R", "R-A", "R-D", "D-R")) %>%
    dplyr::mutate(
        act_dom = if_else(composition %in% c("A-R", "D-R"), d1_Gene, d2_Gene),
        oth_dom = if_else(composition %in% c("A-R", "D-R"), d2_Gene, d1_Gene),
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
p = p + geom_text(data = data.frame(x = 0.85, y = 24, label = bquote("Pearson R=0.79")),
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-1, 4), ylim = c(0, 25), ratio = 5/25)
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Act. Ctrl-Paired" ~ log[2] ~ "(ON:OFF)"))
p = p + ylab("Num. Act. Pairs w/ Reprs.")
p
ggsave("./activators_beating_repressors.pdf", p, height = 2.5, width = 2.5)

# repressors beating activators ------
df %>%
    dplyr::filter(df$composition %in% c("A-R", "R-A", "A-D", "D-A")) %>%
    dplyr::mutate(
        rep_dom = if_else(composition %in% c("A-R", "A-D"), d2_Gene, d1_Gene),
        oth_dom = if_else(composition %in% c("A-R", "A-D"), d1_Gene, d2_Gene),
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
p = p + geom_text(data = data.frame(x = -0.2, y = 15, label = bquote("Pearson R=-0.54")),
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-3.5, 1.5), ylim = c(0, 15), ratio = 5/15)
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Rep. Ctrl-Paired" ~ log[2] ~ "(ON:OFF)"))
p = p + ylab("Num. Repr. Pairs w/ Acts.")
p
ggsave("./repressors_beating_activators.pdf", p, height = 2.5, width = 2.5)

# combinations -----

df %>%
    dplyr::mutate(
        d1_combo_type = if_else(
            d1_baseline_type == "Non-hit",
            if_else(str_detect(d1_Description, "Act"), "Act", "Rep"),
            d1_baseline_type
        ),
        d2_combo_type = if_else(
            d2_baseline_type == "Non-hit",
            if_else(str_detect(d2_Description, "Act"), "Act", "Rep"),
            d2_baseline_type
        ),
        combo_type = paste(
            stringr::str_sub(d1_combo_type, 1, 1),
            stringr::str_sub(d2_combo_type, 1, 1),
            sep = "-")
    ) %>%
    dplyr::filter(combo_type %in% c("A-R", "R-A", "A-D", "D-A", "R-D", "D-R", "D-D")) |>
    dplyr::mutate(
        act_dom = case_when(
            combo_type == "A-R" ~ d1_Gene,
            combo_type == "R-A" ~ d2_Gene,
            combo_type == "A-D" ~ d1_Gene,
            combo_type == "D-A" ~ d2_Gene,
            combo_type == "D-R" ~ d1_Gene,
            combo_type == "R-D" ~ d2_Gene,
            combo_type == "D-D" ~ d1_Gene,
        ),
        act_dom = unlist(map(stringr::str_split(act_dom, " "), first)),
        rep_dom = case_when(
            combo_type == "A-R" ~ d2_Gene,
            combo_type == "R-A" ~ d1_Gene,
            combo_type == "A-D" ~ d2_Gene,
            combo_type == "D-A" ~ d1_Gene,
            combo_type == "D-R" ~ d2_Gene,
            combo_type == "R-D" ~ d1_Gene,
            combo_type == "D-D" ~ d2_Gene,
        ),
        rep_dom = unlist(map(stringr::str_split(rep_dom, " "), first)),
        act_val = if_else(act_dom == d1_Gene, d1_med_d2, d2_med_d2),
        rep_val = if_else(rep_dom == d1_Gene, d1_med_d5, d2_med_d5),
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

get_activation = function(a) {
    comb %>%
        dplyr::filter(str_detect(d1_Gene, a)) %>%
        dplyr::select(d1_med_d2) -> smola
    return(mean(smola$d1_med_d2, na.rm = TRUE))
}

get_repression = function(r) {
    comb %>%
        dplyr::filter(str_detect(d1_Gene, r)) %>%
        dplyr::select(d1_med_d5) -> smolr
    return(mean(smolr$d1_med_d5, na.rm = TRUE))
}

rep_annot = purrr::map_dbl(rlist, get_repression)
ria_annot = purrr::map_dbl(rlist, get_activation)

act_annot = purrr::map_dbl(alist, get_activation)
air_annot = purrr::map_dbl(alist, get_repression)

tibble(dom = rlist, score=rep_annot, other=ria_annot) |>
    dplyr::arrange(desc(score)) -> rdf
rlist = rdf$dom
rep_annot = rdf$score
ria_annot = rdf$other

tibble(dom = alist, score=act_annot, other=air_annot) |>
    dplyr::arrange(score) -> adf
alist = adf$dom
act_annot = adf$score
air_annot = adf$other

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

colors = structure(
    c('#23bedb', '#dbc60d', "#0fa616", '#999999'),
    names = c("Repressor", "Activator", "Dual", "Neither")
)

rep_col_fun = circlize::colorRamp2(c(-4, 2, 4),
                                   c(d5_color, "#dddddd", "#dddddd"))
act_col_fun = circlize::colorRamp2(c(-4, -2, 4),
                                   c("#dddddd", "#dddddd", d2_color))
get_gene_color = function(g) {
    df |>
        dplyr::filter(str_detect(d1_Gene, g)) |>
        dplyr::pull(d1_baseline_type) -> type_list
    gtype = type_list[1]
    if (gtype == "Activator") {
        return("#DBC60D")
    } else if (gtype == "Repressor") {
        return("#23BEDB")
    } else if (gtype == "Dual") {
        return("#0FA616")
    } else {
        return("#999999")
    }
}
a_name_colors = purrr::map_chr(alist, get_gene_color)
r_name_colors = purrr::map_chr(rlist, get_gene_color)

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

padding = unit(2, "mm")
ht_opt$COLUMN_ANNO_PADDING = padding
ht_opt$ROW_ANNO_PADDING = padding

ht = Heatmap(
    amat,
    col = colors,
    rect_gp = gpar(col = "white", lwd = 0.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    na_col = "#FFFFFF",
    show_heatmap_legend = FALSE,
    row_title = "Activators and Dual-Functional Effectors",
    row_order = rev(alist),
    row_gap = unit(2, "mm"),
    row_names_side = 'left',
    row_names_gp = gpar(col=rev(a_name_colors)),
    column_title = "Repressors and Dual-Functional Effectors",
    column_title_side = "bottom",
    column_names_gp = gpar(col=r_name_colors),
    cluster_column_slices = FALSE,
    cluster_row_slices = FALSE,
    column_order = rlist,
    column_gap = unit(2, "mm"),
    bottom_annotation = HeatmapAnnotation(
        avg_rep = rep_annot,
        avg_act = ria_annot,
        col = c(avg_rep = rep_col_fun, avg_act = act_col_fun),
        na_col = "#ffffff",
        annotation_label = c("Rep →", "Act →"),
        annotation_name_side = "left",
        show_legend = FALSE,
        gp = gpar(col = "white", lwd = 0.5),
        gap = unit(0.25, "mm")
    ),
    left_annotation = rowAnnotation(
        avg_rep = rev(air_annot),
        avg_act = rev(act_annot),
        col = c(avg_rep = rep_col_fun, avg_act = act_col_fun),
        na_col = "#ffffff",
        annotation_label = c("← Rep", "← Act"),
        annotation_name_side = "top",
        annotation_name_rot = 90,
        show_legend = FALSE,
        gp = gpar(col = "white", lwd = 0.5),
        gap = unit(0.25, "mm")
    ),
)
cairo_pdf(
    "./activator_repressor_combos.pdf",
    width = 10.5,
    height = 7,
    family = "Arial")
draw(ht, annotation_legend_list = pd, annotation_legend_side = "top")
dev.off()

ht_opt(RESET = TRUE)

