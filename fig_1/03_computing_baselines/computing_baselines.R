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

# first, filter out/drop domains where the orientation will break them
onectrls = dplyr::filter(df, type == "1 control Pair")
onectrls %>%
    dplyr::mutate(
        d1_Gene = if_else(
            is.na(d1_Gene),
            "RANDCTRL",
            d1_Gene
        ),
        d2_Gene = if_else(
            is.na(d2_Gene),
            "RANDCTRL",
            d2_Gene
        )) -> onectrls
genelist = sort(unique(append(unique(onectrls$d1_Gene), unique(onectrls$d2_Gene))))
hits = list.filter(genelist, . != "RANDCTRL" & . != "DMD")
ctrls = list.filter(genelist, (. == "RANDCTRL" | . == "DMD"))

lablist = sort(unique(append(unique(onectrls$domain1), unique(onectrls$domain2))))
ctr_labs = list.filter(lablist, str_detect(., "control"))
cart = expand_grid(hits, ctr_labs)

get_score_df = function(h, c) {
    df %>%
        dplyr::filter(((d1_Gene == h) & (domain2 == c)) | ((domain1 == c) & (d2_Gene == h))) -> qdf
    return(qdf)
}

carth = cart$hits
cartc = cart$ctr_labs
cart_dfs = purrr::map(
    seq(length(carth)),
    function(x) {return(get_score_df(carth[x], cartc[x]))}
)
cart_df = bind_rows(cart_dfs)
cart_df %>%
    dplyr::mutate(
        control = if_else(
            str_detect(domain1, "control"),
            domain1,
            domain2
        ),
        hit_gene = if_else(
            str_detect(domain1, "control"),
            d2_Gene,
            d1_Gene
        ),
        ctrl_character = if_else(
            str_detect(domain1, "control"),
            d2_Description,
            d1_Description
        ),
        order = if_else(
            str_detect(domain1, "control"),
            "C-H",
            "H-C"
        )
    ) -> cart_df

cart_df$hit_gene = as.factor(cart_df$hit_gene)
cart_df %>%
    dplyr::select(c(hit_gene, control, ctrl_character, order,
                    # hit_score_d2, hit_score_d5, hit_score_d10,
                    avg_enrichment_d2, avg_enrichment_d5)) %>%
    tidyr::pivot_longer(
        cols = c(avg_enrichment_d2, avg_enrichment_d5),
        names_to = "day",
        names_prefix = "avg_enrichment_d",
        values_to = "avg_pair_enrichment"
    ) %>%
    dplyr::distinct() -> order_df
order_df %>%
    dplyr::group_by(
        hit_gene, order, day
    ) %>%
    tidyr::nest() %>%
    tidyr::spread(key = order, value = data) %>%
    dplyr::mutate(
        t_test = map2(`C-H`, `H-C`, ~{ t.test(.x$avg_pair_enrichment, 
                                              .y$avg_pair_enrichment,
                                              na.action = na.omit()
        ) %>% broom::tidy() })
    ) %>%
    tidyr::unnest(cols = t_test) -> pdf
pdf %>%
    dplyr::mutate(
        significant = p.value <= 0.05/nrow(pdf),
        dlab = if_else(day == 2, "Act.", "Rep.")
    ) -> pdf

p = ggplot(data = pdf)
p = p + geom_point(
    aes(
        x = estimate,
        y = -log10(p.value),
        color = significant
    )
)
p = p + geom_text_repel(
    aes(
        x = estimate,
        y = -log10(p.value),
        label = ifelse(significant, paste(hit_gene, ", ", dlab, sep=""), ''),
    ),
    size = 3,
    xlim = c(-2, 2),
    ylim = c(4, 7.75),
    # force_pull = 0.005,
    force = 10,
    # segment.size = 0.375,
    # segment.curvature = -0.1,
    # segment.ncp = 7,
    # segment.square = FALSE,
    # segment.inflect = TRUE,
    # min.segment.length = 0,
    # point.padding = 0.25,
    # box.padding = 1.0,
    # max.overlaps = Inf,
    max.time = 5,
    max.iter = 100000
)
p = p + geom_hline(
    yintercept=-log10(0.05/nrow(pdf)),
    linetype='longdash',
    color = "#676767"
)
p = p + coord_fixed(ratio = 4/8, xlim = c(-2, 2), ylim = c(0, 8))
p = p + scale_color_manual(values = c("#999999", "#f54278"), guide='none')
p = p + ylab(bquote(-log[10](p)))
p = p + xlab("C:H mean score - H:C mean score")
p = p + theme_bw() + theme(panel.grid = element_blank())
ggsave("./p_values_orientation_deltas.pdf", p, height = 3, width = 3)


pdf %>%
    dplyr::filter(significant) %>%
    dplyr::select(hit_gene, day) %>%
    dplyr::mutate(hd = paste(hit_gene, day, sep = "_")) %>%
    dplyr::select(hd) -> droplist_df
droplist = droplist_df$hd
df %>%
    dplyr::filter(type == "1 control Pair") -> c1df
c1df %>%
    dplyr::filter(
        paste(d1_Gene, 2, sep = "_") %in% droplist | paste(d2_Gene, 2, sep = "_") %in% droplist
    ) %>%
    dplyr::mutate(
        score = avg_enrichment_d2,
        order = if_else(str_detect(domain1, "control"), "C-H", "H-C"),
        hit_domain = if_else(order == "C-H", d2_Gene, d1_Gene),
        hit_char = if_else(order != "C-H", d1_Description, d2_Description),
        day = "Act."
    ) -> d2_vdf
c1df %>%
    dplyr::filter(
        paste(d1_Gene, 5, sep = "_") %in% droplist | paste(d2_Gene, 5, sep = "_") %in% droplist
    ) %>%
    dplyr::mutate(
        score = avg_enrichment_d5,
        order = if_else(str_detect(domain1, "control"), "C-H", "H-C"),
        hit_domain = if_else(order == "C-H", d2_Gene, d1_Gene),
        hit_char = if_else(order != "C-H", d1_Description, d2_Description),
        day = "Repr."
    ) -> d5_vdf

dplyr::bind_rows(d2_vdf, d5_vdf) %>%
    dplyr::select(hit_domain, order, score, hit_char, day) %>%
    dplyr::mutate(dom_lab = paste(hit_domain, day, sep=", ")) -> sig_df

p = ggplot(sig_df)
p = p + geom_boxplot(
    aes(
        x = dom_lab,
        y = score,
        fill = order
    )
)
p = p + labs(y = bquote(log[2](ON:OFF)), x="")
p = p + theme_bw() + theme(panel.grid = element_blank())
p = p + coord_fixed(ratio = 1/3)
ggsave("./orientation_boxes.pdf", p, height = 3, width = 6)

# Now, the filtering
df %>%
    dplyr::mutate(
        avg_enrichment_d2 = if_else(d2_Gene == "KMT2B", NA_real_, avg_enrichment_d2), # drop C-H
        avg_enrichment_d5 = if_else(d1_Gene == "HERC2", NA_real_, avg_enrichment_d5), # drop H-C
        avg_enrichment_d5 = if_else(d2_Gene %in% c("CBX7", "CHD3", "CREM"), NA_real_, avg_enrichment_d5), # drop C-H
    ) -> df


# Computing baseline scores
oligos = read_csv('../01_raw_counts/csvs/base_oligo_library.csv')

get_df = function (o) { 
    df %>%
        dplyr::filter(domain1 == o | domain2 == o) %>%
        dplyr::filter(type == "1 control Pair") %>%
        dplyr::summarise(
            domain = o,
            avg_d2 = mean(avg_enrichment_d2, na.rm = TRUE),
            med_d2 = median(avg_enrichment_d2, na.rm = TRUE),
            sd_d2 = sd(avg_enrichment_d2, na.rm = TRUE),
            avg_d5 = mean(avg_enrichment_d5, na.rm = TRUE),
            med_d5 = median(avg_enrichment_d5, na.rm = TRUE),
            sd_d5 = sd(avg_enrichment_d5, na.rm = TRUE),
        ) -> odf
    return(odf)
}

baseline_df = bind_rows(purrr::map(oligos$label, get_df))
baseline_df$description = oligos$Description
baseline_df %>%
    dplyr::mutate(
        description = case_when(
            str_detect(description, "Control") ~ "Control",
            str_detect(description, "Act") ~ "Activator",
            str_detect(description, "Repr") ~ "Repressor"
        )
    ) -> baseline_df

# call hits
baseline_df %>%
    dplyr::filter(description == "Control") %>%
    dplyr::summarise(
        mean_d2 = mean(med_d2, na.rm = TRUE),
        sd_d2 = sd(med_d2, na.rm = TRUE),
        mean_d5 = mean(med_d5, na.rm = TRUE),
        sd_d5 = sd(med_d5, na.rm = TRUE),
    ) -> threshold_df
d2_baseline_threshold = threshold_df$mean_d2[1] + 2 * threshold_df$sd_d2[1]
d5_baseline_threshold = threshold_df$mean_d5[1] - 2 * threshold_df$sd_d5[1]

# label hits
baseline_df %>%
    dplyr::mutate(
        baseline_type = case_when(
            description == "Control" ~ "Control",
            (med_d2 >= d2_baseline_threshold) & (med_d5 > d5_baseline_threshold) ~ "Activator",
            (med_d2 < d2_baseline_threshold) & (med_d5 <= d5_baseline_threshold) ~ "Repressor",
            TRUE ~ "Non-hit"
        )
    ) -> baseline_df
baseline_df$baseline_type = as.factor(baseline_df$baseline_type)
# levels(baseline_df$baseline_type) = c("Control", "Activator", "Repressor", "Non-hit")

# now we draw them
# acr_colors = c('#DBC60D', "#999999", "#23BEDB")
acr_colors = c("Control" = "#999999", 
               "Activator" = "#DBC60D", 
               "Repressor" = "#23BEDB",
               "Non-hit" = "#DEDEDE")
p = ggplot(data = baseline_df)
p = p + coord_fixed(ratio = 1, xlim = c(-3, 3.5), ylim = c(-3.5, 3))
p = p + geom_hline(yintercept = d5_baseline_threshold,
                   linetype = "dashed")
p = p + geom_vline(xintercept = d2_baseline_threshold,
                   linetype = "dashed")
p = p + theme_bw()

p = p + geom_point(aes(x = med_d2,
                       y = med_d5,
                       color = baseline_type))
p = p + labs(x = "Act. Control-paired log(ON:OFF)",
             y = "Rep. Control-paired log(ON:OFF)")
p = p + scale_color_manual(name = "", values = acr_colors)
p = p + guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave(
    "./baseline_d2_d5_scatter.pdf",
    p, 
    height = 3,
    width = 4)

write_csv(baseline_df, "./baseline_scores.csv")

# correlation with prior data
pfam = read_csv('../01_raw_counts/csvs/NucPfam_ReprActStbl_data.csv')
pfam %>%
    dplyr::select(
        label, "Avg ReprD5", "Avg Act"
    ) %>%
    dplyr::rename(
        # new = old
        avg_d5_pri = "Avg ReprD5",
        avg_d2_pri = "Avg Act"
    ) -> pfam
baseline_df %>%
    dplyr::left_join(
        pfam,
        by = c("domain" = "label")
    ) -> baseline_priored_df


p = ggplot(baseline_priored_df,
           aes(x = med_d2, y = -avg_d2_pri, color = baseline_type))
p = p + geom_point(size = 1.5)
p = p + geom_text(data = data.frame(x = -1.85, y = 9, label = "Pearson r = 0.81"),
          mapping = aes(x = x, y = y, label = label),
          inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-7, 9), ratio = 7/16)
p = p + xlab("Med. Control-Paired Act. log(ON:OFF)") + ylab("Prior Pfam Screen Act. log(ON:OFF)")
p = p + theme_linedraw() + theme(panel.grid = element_blank(), legend.position = c(0.775, 0.3))
p = p + guides(color = guide_legend(override.aes = list(size = 3)))
p = p + scale_color_manual(name = "", values = acr_colors)
ggsave("./nucpfam_correlation_act.pdf", p, width = 3.2, height = 3.2)


p = ggplot(baseline_priored_df,
           aes(x = med_d5, y = -avg_d5_pri, color = baseline_type))
p = p + geom_point(size = 1.5)
p = p + geom_text(data = data.frame(x = -1.85, y = 1.5, label = "Pearson r = 0.86"),
          mapping = aes(x = x, y = y, label = label),
          inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-3.5, 3.5), ylim = c(-7, 1.5), ratio = 7/8.5)
p = p + xlab("Med. Control-Paired Rep. log(ON:OFF)") + ylab("Prior Pfam Screen Rep. log(ON:OFF)")
p = p + theme_linedraw() + theme(panel.grid = element_blank(), legend.position = c(0.775, 0.3))
p = p + guides(color = guide_legend(override.aes = list(size = 3)))
p = p + scale_color_manual(name = "", values = acr_colors)
ggsave("./nucpfam_correlation_rep.pdf", p, width = 3.2, height = 3.2)
