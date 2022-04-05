# Imports
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
library(FactoMineR)

d5_color = "#28bedb"
d2_color = "#dbc60d"

# add prior information to the dataframe
df = readr::read_csv('../01_raw_counts/csvs/concatenation_filtered_scored_pairs.csv')
oligos = read_csv('../01_raw_counts/csvs/base_oligo_library.csv')
pfam = read_csv('../01_raw_counts/csvs/NucPfam_ReprActStbl_data.csv')
siltile = read_csv('../01_raw_counts/csvs//SilTile_dfAnalysis.csv')

pfam %>%
    dplyr::select(
        label, "Avg ReprD5", "Avg ReprD9", "Avg Act"
    ) %>%
    dplyr::rename(
        # new = old
        avg_d5 = "Avg ReprD5",
        avg_dm = "Avg ReprD9",
        avg_d2 = "Avg Act"
    ) -> pfam
siltile %>%
    dplyr::select(
        label, "Avg D5", "Avg D13"
    ) %>%
    dplyr::rename(
        avg_d5 = "Avg D5",
        avg_dm = "Avg D13"
    ) %>%
    tibble::add_column(
        avg_d2 = NA
    ) -> siltile
prior = dplyr::bind_rows(pfam, siltile)
p1 = copy(prior)
p2 = copy(prior)
colnames(p1) = paste("d1_prior", colnames(p1), sep="_")
colnames(p2) = paste("d2_prior", colnames(p2), sep="_")
df %>%
    dplyr::left_join(
        p1,
        c(domain1 = "d1_prior_label"),
        copy = TRUE,
    ) %>%
    dplyr::left_join(
        p2,
        c(domain2 = "d2_prior_label"),
        copy = TRUE
    ) -> dfx
write_csv(dfx, "./02_pairs_priored.csv")

# eliminate obvious outlier controls
order_df = readr::read_csv("./01_one_control_ordered_df.csv")
order_df %>%
    dplyr::mutate(
        pair_filtered_enrichment = case_when(
            (str_detect(ctrl_character,"Act")) & (day == 2) ~ avg_pair_enrichment,
            (str_detect(ctrl_character,"Rep")) & (day == 5) ~ avg_pair_enrichment,
            TRUE ~ NA_real_
        ),
    ) %>%
    tidyr::drop_na(pair_filtered_enrichment) -> order_df
order_df %>%
    dplyr::mutate(day = as.numeric(day)) %>%
    dplyr::group_by(control, day) %>%
    tidyr::nest()  %>%
    tidyr::spread(key = day, value = data) %>%
    dplyr::mutate(norm_fit_d2 = map(`2`,
                                    function(x) { 
                                        return(fitdistrplus::fitdist(x$avg_pair_enrichment, 
                                                                     "norm", 
                                                                     "mle")$estimate)}),
                  norm_fit_d5 = map(`5`,
                                    function(x) {
                                        return(fitdistrplus::fitdist(x$avg_pair_enrichment, 
                                                                     "norm",
                                                                     "mle")$estimate)})) %>%
    tidyr::unnest_wider(norm_fit_d2, names_sep="_") %>%
    tidyr::unnest_wider(norm_fit_d5, names_sep="_") %>%
    tidyr::pivot_longer(
        cols = c(norm_fit_d2_mean, norm_fit_d5_mean),
        names_to = "day",
        names_pattern = ".*d(.).*",
        values_to = "mean_score"
    ) %>%
    dplyr::select(control, day, mean_score) %>%
    dplyr::ungroup() -> pdf

pdf %>%
    dplyr::group_by(day) %>%
    tidyr::nest() %>%
    tidyr::spread(key = day, value = data) -> tdf

tdf %>% 
    dplyr::mutate(
        mean_d2 = map(`2`, function(x) {return(mean(x$mean_score, na.rm = TRUE))}),
        sd_d2 = map(`2`, function(x) {return(sd(x$mean_score, na.rm = TRUE))}),
        mean_d5 = map(`5`, function(x) {return(mean(x$mean_score, na.rm = TRUE))}),
        sd_d5 = map(`5`, function(x) {return(sd(x$mean_score, na.rm = TRUE))}),
        mean_d2 = mean_d2[[1]],
        sd_d2 = sd_d2[[1]],
        mean_d5 = mean_d5[[1]],
        sd_d5 = sd_d5[[1]]
    ) -> tdf

c_avg_d2 = -0.302
c_std_d2 = 0.511
c_avg_d5 = -0.531
c_std_d5 = 0.565

pdf %>%
    dplyr::mutate(
        c_mean = if_else(day == 2, c_avg_d2, c_avg_d5),
        c_std = if_else(day == 2, c_std_d2, c_std_d5),
        sig = (mean_score >= c_mean + 2*c_std) | (mean_score <= c_mean - 2*c_std)
    ) %>%
    dplyr::mutate(screen = if_else(day == 2, "Activation", "Repression")) -> pdf
p = ggplot(data = pdf)
p = p + geom_point(
    aes(
        x = control,
        y = mean_score,
        color = sig
    )
)
p = p + geom_hline(aes(yintercept = c_mean + 2 * c_std),
                   linetype = "dashed")
p = p + geom_hline(aes(yintercept = c_mean - 2 * c_std),
                   linetype = "dashed")
p = p + geom_hline(aes(yintercept = c_mean),
                   linetype = "dashed", color = "#999999")
p = p + facet_wrap(facets = vars(screen))
p = p + labs(x = "", y = "Average log(ON:OFF)")
p = p + coord_fixed(ylim = c(-2, 2), ratio = 1)
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
),
legend.position = "none")
ggsave("./02_control_outliers.pdf", p, height = 4, width = 8)
df = readr::read_csv("./02_pairs_priored.csv")
df %>%
    dplyr::mutate(
        avg_enrichment_d2 = if_else(
            str_detect(pair, "DMD_control_tiles;ENSG00000198947;;31;"),
            NA_real_,
            avg_enrichment_d2
        )) %>%
    dplyr::mutate(
        avg_enrichment_d5 = if_else(
            str_detect(pair, "DMD_control_tiles;ENSG00000198947;;31;"),
            NA_real_,
            avg_enrichment_d5
        )) %>%
    dplyr::mutate(
        avg_enrichment_d5 = if_else(
            str_detect(pair, "Random_control;;;89"),
            NA_real_,
            avg_enrichment_d5
        )) -> df
df %>% # drop d10 data
    dplyr::select(-one_of(Filter(function(x){str_detect(x, "_d10")}, names(df)))) -> df
readr::write_csv(df, "./02_pairs_priored.csv")

# get list of controls
df = readr::read_csv('./02_pairs_priored.csv')
df %>%
    dplyr::filter(type == "2 control Pair") -> cdf
controls = unique(cdf$domain1)

# get list of non-control domains
df %>%
    dplyr::filter(type == "0 control Pair") -> ddf
domains = unique(ddf$domain1)

# for each control, get its avg. score with each domain
get_avg_score = function(control, domain) {
    df %>%
        dplyr::filter((domain1 == control & domain2 == domain) |
                          (domain2 == control & domain1 == domain)) %>%
        dplyr::summarise(
            control = control,
            domain = domain,
            avg_act = mean(avg_enrichment_d2, na.rm = TRUE),
            avg_rep = mean(avg_enrichment_d5, na.rm = TRUE),
            pri_act = mean(if_else(domain1 == control, d2_prior_avg_d2, d1_prior_avg_d2), na.rm = TRUE),
            pri_rep = mean(if_else(domain1 == control, d2_prior_avg_d5, d1_prior_avg_d5), na.rm = TRUE)
        ) -> smoldf
    return(smoldf)
}
all_pairs = expand.grid(controls, domains)
pair_ctls = all_pairs$Var1
pair_doms = all_pairs$Var2
pair_scores = map2_dfr(pair_ctls, pair_doms, get_avg_score)

# building the activator matrix
pair_scores %>%
    dplyr::mutate(
        control = as.character(control),
        domain = as.character(domain)
    ) %>%
    dplyr::mutate(control = if_else(
        str_detect(control, "DMD"),
        paste("DMD Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
        paste("Random Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
    )) %>%
    dplyr::select(control, domain, avg_act) %>%
    tidyr::pivot_wider(
        id_cols = control, #vertical
        names_from = domain, #horizontal
        values_from = avg_act
    ) %>%
    tibble::column_to_rownames("control") -> act_mat
dend = as.dendrogram(hclust(dist(act_mat)), method = "ward.D2")
dend = color_branches(dend, k=3, col = c("#7570b3", "#d95f02", "#1b9e77"))

col_fun = circlize::colorRamp2(c(-5,-1, 5),
                               c("#c2d3d6", "#dddddd", d2_color)) 
tm = cbind(act_mat) # copy into new memory location
colnames(tm) = NULL
rownames(tm) = NULL
colmeans = colMeans(tm, na.rm = TRUE)
pair_scores %>%
    dplyr::mutate(
        control = as.character(control),
        domain = as.character(domain)
    ) %>%
    dplyr::mutate(control = if_else(
        str_detect(control, "DMD"),
        paste("DMD Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
        paste("Random Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
    )) %>%
    dplyr::select(control, domain, pri_act) %>%
    tidyr::pivot_wider(
        id_cols = control, #vertical
        names_from = domain, #horizontal
        values_from = pri_act 
    ) %>%
    tibble::column_to_rownames("control") -> act_pri_mat
colnames(act_pri_mat) = NULL
rownames(act_pri_mat) = NULL
primeans = colMeans(act_pri_mat, na.rm = TRUE)
h = Heatmap(
    act_mat,
    name = "Activation",
    cluster_rows = dend,
    row_split = 3,
    row_title = " ",
    column_title = "← More Activating Effector Domains",
    row_title_side = "right",
    row_names_side = "left",
    row_dend_side = "right",
    show_column_names = FALSE,
    show_column_dend = FALSE,
    col = col_fun,
    heatmap_legend_param = list(
        col_fun = col_fun,
        title = bquote(atop("Activation", log[2] ~ "(ON:OFF)")),
        title_position = "topcenter",
        legend_height = unit(4, "cm"),
        na_col = c("#ffffff", "#ffffff"),
        at = c(-5,-2.5, 0, 2.5, 5)
    ),
    na_col = "#ffffff",
    top_annotation = HeatmapAnnotation(
        avg = colmeans,
        pri = primeans * -1,
        col = list(avg = col_fun, pri = col_fun),
        na_col = c("#ffffff", "#ffffff"),
        annotation_label = c("Average", "Prior Screen Data"),
        annotation_name_side = "left",
        show_legend = FALSE
    ),
)
cairo_pdf(
    "./02_control_activator_clusters.pdf",
    width = 9,
    height = 4.5,
    family = "Arial")
draw(h, padding = unit(c(2, 6, 2, 10), "mm"))
dev.off()

# now for repressors
pair_scores %>%
    dplyr::mutate(
        control = as.character(control),
        domain = as.character(domain)
    ) %>%
    dplyr::mutate(control = if_else(
        str_detect(control, "DMD"),
        paste("DMD Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
        paste("Random Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
    )) %>%
    dplyr::select(control, domain, avg_rep) %>%
    tidyr::pivot_wider(
        id_cols = control, #vertical
        names_from = domain, #horizontal
        values_from = avg_rep
    ) %>%
    dplyr::filter(if_any(everything(), ~ !is.nan(.))) %>%
    tibble::column_to_rownames("control") -> rep_mat
rep_mat = subset(rep_mat, rownames(rep_mat) != "Random Control 89") # all NaNs, toss
dend = as.dendrogram(hclust(dist(rep_mat)), method = "ward.D2")
dend = color_branches(dend, k=3, col = c("#7570b3", "#d95f02", "#1b9e77"))

col_fun = circlize::colorRamp2(c(-5,1, 5),
                               c(d5_color, "#dddddd", "#dbd7af")) 

# histogram on top
tm = cbind(rep_mat) # copy into new memory location
colnames(tm) = NULL
rownames(tm) = NULL
colmeans = colMeans(tm, na.rm = TRUE)

# prior repression scores
pair_scores %>%
    dplyr::mutate(
        control = as.character(control),
        domain = as.character(domain)
    ) %>%
    dplyr::mutate(control = if_else(
        str_detect(control, "DMD"),
        paste("DMD Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
        paste("Random Control", purrr::map(control, function(x) { strsplit(x, ';')[[1]][4] })),
    )) %>%
    dplyr::select(control, domain, pri_rep) %>%
    tidyr::pivot_wider(
        id_cols = control, #vertical
        names_from = domain, #horizontal
        values_from = pri_rep
    ) %>%
    tibble::column_to_rownames("control") -> rep_pri_mat
colnames(act_pri_mat) = NULL
rownames(act_pri_mat) = NULL
primeans = colMeans(rep_pri_mat, na.rm = TRUE)
h = Heatmap(
    rep_mat,
    name = "Repression",
    cluster_rows = dend,
    row_split = 3,
    row_title = " ",
    column_title = "More Repressing Effector Domains →",
    row_title_side = "right",
    row_names_side = "left",
    row_dend_side = "right",
    show_column_names = FALSE,
    show_column_dend = FALSE,
    col = col_fun,
    heatmap_legend_param = list(
        col_fun = col_fun,
        title = bquote(atop("Repression", log[2] ~ "(ON:OFF)")),
        title_position = "topcenter",
        legend_height = unit(4, "cm"),
        at = c(-5,-2.5, 0, 2.5, 5)
    ),
    na_col = "#ffffff",
    top_annotation = HeatmapAnnotation(
        avg = colmeans,
        pri = -1 * primeans,
        col = list(avg = col_fun, pri = col_fun),
        na_col = c("#ffffff", "#ffffff"),
        annotation_label = c("Average", "Prior Screen Data"),
        annotation_name_side = "left",
        show_legend = FALSE
    ),
)
cairo_pdf(
    "./02_control_repressor_clusters.pdf",
    width = 9,
    height = 4.5,
    family = "Arial")
draw(h, padding = unit(c(2, 6, 2, 10), "mm"))
dev.off()

# want to filter out certain controls
# activators -- Rand 89, DMD 130, Rand 257, DMD 297
# repressors -- DMD 182, DMD 297, DMD 130 
act_cont_filter_list = c(
    "Random_control;;;89;",
    "Random_control;;;257;",
    "DMD_control_tiles;ENSG00000198947;;130;",
    "DMD_control_tiles;ENSG00000198947;;297;"
)
rep_cont_filter_list = c(
    "DMD_control_tiles;ENSG00000198947;;130;",
    "DMD_control_tiles;ENSG00000198947;;182;",
    "DMD_control_tiles;ENSG00000198947;;297;"
)
df = readr::read_csv('./02_pairs_priored.csv')
df %>%
    dplyr::mutate(
        avg_enrichment_d2 = if_else(
            (domain1 %in% act_cont_filter_list) | (domain2 %in% act_cont_filter_list),
            NA_real_,
            avg_enrichment_d2
        ),
        avg_enrichment_d5 = if_else(
            (domain1 %in% rep_cont_filter_list) | (domain2 %in% rep_cont_filter_list),
            NA_real_,
            avg_enrichment_d5
        )
    ) -> fdf
readr::write_csv(fdf, "./02_pairs_priored_controlfiltered.csv")
readr::write_csv(df, "./02_pais_priored_control_UNfiltered.csv")
