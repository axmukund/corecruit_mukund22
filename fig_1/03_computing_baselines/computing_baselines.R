
# imports -------------------------------------------------------------------------------------



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
library(RColorBrewer)
library(ggExtra)

d2_color = "#28bedb"
d5_color = "#dbc60d"
d10_color = "#db9d0d"

act_score_thresh = -0.0236
rep_score_thresh = 0.771

acr_colors = c("Control" = "#999999",
               "Activator" = "#DBC60D",
               "Repressor" = "#23BEDB",
               "Dual" = "#0FA616",
               "Non-hit" = "#DEDEDE")
adr_colors = c("Activator" = "#DBC60D",
               "Repressor" = "#23BEDB",
               "Dual" = "#0FA616",
               "Non-hit" = "#DEDEDE")
ar_colors = c("Activator" = "#DBC60D",
              "Repressor" = "#23BEDB")

# cartesian dataframe construction ------------------------------------------------------------



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


# orientation ---------------------------------------------------------------------------------



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


# control variability plot --------------------------------------------------------------------

genes = c("Control", "FOXO3", "ZNF10", "CRTC2")
d2_threshold  = 0.639 # i know this because i am from the future
d5_threshold  = 0.103 # just kidding, this plot was made after the thresholds were computed
vcdf = df %>%
    dplyr::filter(type == "2 control Pair") %>%
    dplyr::mutate(hit = "Control",
                  d2 = avg_enrichment_d2,
                  d5 = avg_enrichment_d5,
                  cat = "Control") %>%
    dplyr::select(pair, hit, d2, d5, cat)
vfdf = df %>%
    dplyr::filter(type == "1 control Pair") %>%
    dplyr::mutate(
        hit = if_else(d1_Gene %in% c("DMD", NA), d2_Gene, d1_Gene),
        ctrl = if_else(d1_Gene %in% c("DMD", NA), d1_Gene, d2_Gene)
    ) %>%
    dplyr::filter(hit %in% genes) %>%
    dplyr::mutate(d2 = avg_enrichment_d2,
                  d5 = avg_enrichment_d5,
                  cat = hit) %>%
    dplyr::select(pair, hit, d2, d5, cat)
vdf = rbind(vcdf, vfdf)

p = ggplot(data = vdf) +
    geom_point(aes(x = d2, y = d5, color = cat)) +
    scale_color_manual(name = "", values = c("#CCCCCC", "#DBC60D", "#0FA616", "#23BEDB"),
                       labels=c("Control + Control", "CRTC2 + Control", "FOXO3 + Control", "ZNF10 + Control")) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    coord_fixed(ratio = 1,
                xlim = c(-3.5, 5),
                ylim = c(-5, 3.5)) +
    geom_hline(yintercept = 0.4, #d5_baseline_threshold,
               linetype = "dashed") +
    geom_vline(xintercept = 0.0, #d2_baseline_threshold,
               linetype = "dashed") +
    theme_linedraw() +
    theme(
        panel.grid = element_blank(),
        legend.position = c(0.75, 0.225),
        legend.text = element_text(size = 8),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(3, "mm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
    ) +
    labs(x = "Act. log(ON:OFF)", y = "Repr. log(ON:OFF)")
p
ggsave("./baseline_variability_plot.pdf", p, width = 3, height = 3)

# Computing baseline scores -----------------------------------------------
oligos = read_csv('../01_raw_counts/csvs/base_oligo_library.csv')

act_thr = 0.19
rep_thr = 0.19

get_df = function (o) {
    df %>%
        dplyr::filter(domain1 == o | domain2 == o) %>%
        dplyr::filter(type == "1 control Pair") %>%
        dplyr::summarise(
            domain = o,
            avg_d2 = mean(avg_enrichment_d2, na.rm = TRUE),
            med_d2 = median(avg_enrichment_d2, na.rm = TRUE),
            sd_d2 = sd(avg_enrichment_d2, na.rm = TRUE),
            min_d2 = min(avg_enrichment_d2, na.rm = TRUE),
            max_d2 = max(avg_enrichment_d2, na.rm = TRUE),
            delta_min_d2 = med_d2 - min_d2,
            delta_max_d2 = max_d2 - min_d2,
            avg_d5 = mean(avg_enrichment_d5, na.rm = TRUE),
            med_d5 = median(avg_enrichment_d5, na.rm = TRUE),
            sd_d5 = sd(avg_enrichment_d5, na.rm = TRUE),
            min_d5 = min(avg_enrichment_d5, na.rm = TRUE),
            max_d5 = max(avg_enrichment_d5, na.rm = TRUE),
            delta_min_d5 = med_d5 - min_d5,
            delta_max_d5 = max_d5 - min_d5,
            # frac_d2 = mean(avg_enrichment_d2 >= -0.9379871, na.rm = TRUE), # d2_baseline_threshold
            # frac_d5 = mean(avg_enrichment_d5 <= -0.3797294, na.rm = TRUE) # d5_baseline_threshold
            frac_d2 = mean(avg_enrichment_d2 >= act_score_thresh, na.rm = TRUE), # d2_baseline_threshold
            frac_d5 = mean(avg_enrichment_d5 <= rep_score_thresh, na.rm = TRUE), # d5_baseline_threshold
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

baseline_df %>%
    dplyr::filter(description == "Control") %>%
    dplyr::summarise(
        mean_d2 = mean(med_d2, na.rm = TRUE),
        sd_d2 = sd(med_d2, na.rm = TRUE),
        mean_d5 = mean(med_d5, na.rm = TRUE),
        sd_d5 = sd(med_d5, na.rm = TRUE),
    ) -> threshold_df
d2_baseline_threshold = threshold_df$mean_d2[1] + 2 * threshold_df$sd_d2[1] # we don't use these
d5_baseline_threshold = threshold_df$mean_d5[1] - 2 * threshold_df$sd_d5[1]

get_act_thresholded_doms_df = function(o) {
    df %>%
        dplyr::filter(domain1 == o | domain2 == o) %>%
        dplyr::filter(type == "1 control Pair") %>%
        dplyr::filter(avg_enrichment_d2 >= 0) %>%
        dplyr::summarise(
            domain = o,
            act_avg_d2 = mean(avg_enrichment_d2, na.rm = TRUE),
            act_med_d2 = median(avg_enrichment_d2, na.rm = TRUE),
            act_sd_d2 = sd(avg_enrichment_d2, na.rm = TRUE),
            act_min_d2 = min(avg_enrichment_d2, na.rm = TRUE),
            act_max_d2 = max(avg_enrichment_d2, na.rm = TRUE),
            act_delta_min_d2 = act_med_d2 - act_min_d2,
            act_delta_max_d2 = act_max_d2 - act_min_d2,
        ) -> odf
    return(odf)
}

get_rep_thresholded_doms_df = function(o) {
    df %>%
        dplyr::filter(domain1 == o | domain2 == o) %>%
        dplyr::filter(type == "1 control Pair") %>%
        dplyr::filter(avg_enrichment_d5 <= 0) %>%
        dplyr::summarise(
            domain = o,
            rep_avg_d5 = mean(avg_enrichment_d5, na.rm = TRUE),
            rep_med_d5 = median(avg_enrichment_d5, na.rm = TRUE),
            rep_sd_d5 = sd(avg_enrichment_d5, na.rm = TRUE),
            rep_min_d5 = min(avg_enrichment_d5, na.rm = TRUE),
            rep_max_d5 = max(avg_enrichment_d5, na.rm = TRUE),
            rep_delta_min_d5 = rep_med_d5 - rep_min_d5,
            rep_delta_max_d5 = rep_max_d5 - rep_min_d5,
        ) -> odf
    return(odf)
}

baseline_act_df = bind_rows(purrr::map(oligos$label, get_act_thresholded_doms_df))
baseline_rep_df = bind_rows(purrr::map(oligos$label, get_rep_thresholded_doms_df))

baseline_df %>%
    dplyr::left_join(
        baseline_act_df,
        by = c("domain"),
        suffix = c("", "_act_threshold")
    ) %>%
    dplyr::left_join(
        baseline_rep_df,
        by = c("domain"),
        suffix=c("", "_rep_threshold")
    ) -> baseline_df

baseline_df %>%
    dplyr::filter(str_detect(domain, "control")) %>%
    dplyr::select(avg_d2, avg_d5) %>%
    dplyr::summarise(
        m2 = mean(avg_d2, na.rm = TRUE),
        s2 = sd(avg_d2, na.rm = TRUE),
        m5 = mean(avg_d5, na.rm = TRUE),
        s5 = mean(avg_d5, na.rm = TRUE)
    ) -> control_threshold_df

baseline_df %>%
    dplyr::mutate(
        is_activator = (frac_d2 >= act_thr) & ((act_avg_d2 >= 0.95) | (avg_d2 >= -0.5)),
        is_repressor = (frac_d5 >= rep_thr) & ((rep_avg_d5 <= 0.20) | (avg_d5 <= -0.3)),
        baseline_type = case_when(
            description == "Control" ~ "Control",
            is_activator & !is_repressor ~ "Activator",
            !is_activator & is_repressor ~ "Repressor",
            is_activator & is_repressor ~ "Dual",
            TRUE ~ "Non-hit"
        ),
        prior_type = case_when(
            str_detect(description, "Act") ~ "Activator",
            str_detect(description, "Rep") ~ "Repressor",
            str_detect(description, "Control") ~ "Control",
            TRUE ~ "None"
        )
    ) -> baseline_df
baseline_df$baseline_type = as.factor(baseline_df$baseline_type)
get_gene = function(d) {
    return(str_replace(str_split(d, ";")[[1]][2], "_HUMAN", ""))
}
baseline_df$gene = purrr::map(baseline_df$domain, get_gene)


# label and plot hits w/ fractional thresholds  ------------------- ------------------- ------------

bdf = dplyr::filter(baseline_df, baseline_type != "Control")

p = ggplot(data = bdf, aes(x = frac_d2, y = act_avg_d2),)
p = p + geom_point(aes(color = baseline_type))
# p = p + geom_smooth(method="lm", se=FALSE)
p = p + geom_vline(xintercept=act_thr)
p = p + geom_hline(yintercept=0.95)
# p = p + xlab("Frac. Activating Control Pairs") + ylab("Density")
p = p + scale_color_manual(
    name = "Prior Screen Effect",
    breaks = c("Control", "Activator", "Repressor", "Dual", "Non-hit"),
    values = acr_colors
)
p = p + coord_fixed(ratio = 1/3)
p = p + theme_linedraw() + theme(
    panel.grid = element_blank(),
    # legend.position = c(0.77, 0.75),
    legend.text = element_text(size = 10),
)
p1 = p

p = ggplot(data = bdf, aes(x = avg_d5, y = rep_avg_d5),)
p = p + geom_point(aes(color = baseline_type))
# p = p + geom_smooth(method="lm", se=FALSE)
p = p + geom_vline(xintercept=rep_thr)
p = p + geom_hline(yintercept=-0.2)
# p = p + xlab("Frac. Activating Control Pairs") + ylab("Density")
p = p + scale_color_manual(
    name = "Prior Screen Effect",
    breaks = c("Control", "Activator", "Repressor", "Dual", "Non-hit"),
    values = acr_colors
)
p = p + coord_fixed(ratio = 2)
p = p + theme_linedraw() + theme(
    panel.grid = element_blank(),
    # legend.position = c(0.77, 0.75),
    legend.text = element_text(size = 10),
)
p2 = p

p = p1 + p2
p
ggsave("./avg_frac_thresholds.pdf", p, height=5, width=10)

# now densit plots ---------------

p = ggplot(data=dplyr::filter(baseline_df, description != "Control"))
p = p + geom_density(
    aes(
        x = frac_d2,
        fill = fct_rev(prior_type)
    ),
    bw = 0.025,
    alpha = 0.75,
)
p = p + geom_vline(aes(xintercept = act_thr))
p = p + xlab("Frac. Activating Control Pairs") + ylab("Density")
p = p + scale_fill_manual(name = "Prior Screen Effect", breaks = c("Activator", "Repressor"), values = ar_colors)
p = p + coord_fixed(ratio = 0.5/6.1, xlim=c(0, 1), ylim=c(0, 6.1))
p = p + scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0))
p = p + scale_y_continuous(breaks=c(0, 3, 6))
p = p + theme_linedraw() + theme(
    panel.grid = element_blank(),
    legend.position = c(0.77, 0.75),
    legend.text = element_text(size = 10),
)
h = 2.25
ggsave("./frac_d2_density.pdf", p, height=h, width = 2 * h)
p1 = p

p = ggplot(data=dplyr::filter(baseline_df, description != "Control"))
p = p + geom_density(
    aes(
        x = frac_d5,
        fill = fct_rev(prior_type)
    ),
    bw = 0.025,
    alpha = 0.75,
)
p = p + geom_vline(aes(xintercept = rep_thr))
p = p + xlab("Frac. Repressing Control Pairs") + ylab("Density")
p = p + scale_fill_manual(name = "Prior Screen Effect", breaks = c("Activator", "Repressor"), values = ar_colors)
p = p + coord_fixed(ratio = 0.5/10, xlim=c(0, 1), ylim=c(0, 10))
p = p + scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0))
p = p + scale_y_continuous(breaks=c(0, 5, 10))
p = p + theme_linedraw() + theme(
    panel.grid = element_blank(),
    legend.position = c(0.77, 0.75),
    legend.text = element_text(size = 10),
)
h = 2.25
ggsave("./frac_d5_density.pdf", p, height=h, width = 2 * h)
p2 = p

h = 4.75
p = p1 / p2
ggsave("./frac_both_density.pdf", p, height=h, width=h)

baseline_df %>%
    dplyr::mutate(
        nlab = if_else(
            (frac_d2 >= 0.25) & (frac_d5 >= 0.25),
            domain,
            ""
        )
    ) %>%
    dplyr::filter(description != "Control") -> tdf

p = ggplot(data = tdf)
p = p + geom_vline(aes(xintercept = act_thr), linetype = "dashed")
p = p + geom_hline(aes(yintercept = rep_thr), linetype = "dashed")
p = p + geom_point(aes(x = frac_d2,
                       y = frac_d5,
                       color = prior_type,),)
p = p + xlab("Frac. Activating Control Pairs") + ylab("Frac. Repressing Control Pairs")
p = p + scale_color_manual(
    name = "Prior Screen Effect",
    breaks = c("Activator", "Repressor"),
    values = ar_colors
)
p = p + coord_fixed(ratio = 1,
                    xlim = c(0, 1),
                    ylim = c(0, 1))
p = p + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1.0)))
p = p + theme_linedraw() + theme(
    panel.grid = element_blank(),
    legend.position = c(0.55, 0.83),
    legend.text = element_text(size = 8),
    legend.box.background = element_rect(color="black"),
    legend.key.height = unit(3.0, "mm"),
    legend.title = element_text(size = 9)
)
# p = p + ggrepel::geom_label_repel(aes(x = frac_d2,y= frac_d5,label = nlab))
p = ggMarginal(
    p,
    groupColour = TRUE,
    groupFill = TRUE,
    type = "density",
    bw = 0.025,
    size = 5
)
p


ggsave("./frac_scatterplot.pdf",
       p,
       height = 3,
       width = 3)
p

# drawing baseline correlations ------------------------------------------------------------


# now we draw them
# acr_colors = c('#DBC60D', "#999999", "#23BEDB")
adr_label_colors = c(
    "Activator" = "#DBC60D",
    "Repressor" = "#23BEDB",
    "Dual" = "#0FA616",
    "Non-hit" = "#999999"
)

p = ggplot(data = bdf) # dplyr::filter(baseline_df, baseline_type != "Control"))
p = p + coord_fixed(ratio = 1,
                    xlim = c(-3, 3.5),
                    ylim = c(-3.5, 3))
p = p + geom_point(aes(x = med_d2, #frac_d2,
                       y = med_d5, #frac_d5,
                       color = baseline_type))
# p = p + geom_hline(yintercept = -0.248, #rep_thr, #d5_baseline_threshold,
#                    linetype = "dashed")
# p = p + geom_vline(xintercept = -0.536, #act_thr, #d2_baseline_threshold,
#                    linetype = "dashed")
p = p + xlab(label = "Med. Ctrl-Paired Act. log(ON:OFF)")
p = p + ylab(label = "Med. Ctrl-Paired Rep. log(ON:OFF)")
p = p + scale_color_manual(
    name = "Effector Type",
    values = adr_label_colors,
    breaks = c("Activator", "Repressor", "Dual", "Non-hit")
)
p = p + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1.0)))
p = p + theme_linedraw() + theme(
    panel.grid = element_blank(),
    legend.position = c(0.78, 0.19),
    legend.text = element_text(size = 7),
    legend.key.height = unit(3.0, "mm"),
    legend.box.background = element_rect(color = "black"),
    legend.title = element_text(size = 9)
)
p = p + geom_text_repel(
    data = dplyr::filter(bdf, gene %in% c("CBX1", "FOXO3", "CRTC2", "U2AF4")),
    aes(
        x = med_d2,
        y = med_d5,
        label = gene,
        color = baseline_type
    ),
    force = 3,
    box.padding = 3,
    max.overlaps = 8,
    # segment.color = "black",
    ylim = c(-2.9, 4.75)#, xlim=c(-3, 3)
)
p

ggsave("./baseline_d2_d5_scatter.pdf",
       p,
       height = 3,
       width = 3)
write_csv(baseline_df, "./baseline_scores.csv")

# pfam  act correlation ----------------------------------------------------------------------------


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

# f%>%
#     dplyr::filter(baseline_type != "Control") -> baseline_priored_df

p = ggplot(baseline_priored_df,
           aes(x = med_d2, y = -avg_d2_pri, color = baseline_type))
# p = p + geom_vline(xintercept = -0.5, color = '#777777', linetype = "dashed")
p = p + geom_hline(yintercept = -1.94, color = '#777777', linetype = "dashed") # estimate from HT-recruit data
p = p + geom_point(size = 1.5)
p = p + geom_text(data = data.frame(x = -2, y = 8, label = "Pearson\nR = 0.82"),
          mapping = aes(x = x, y = y, label = label),
          inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-3, 3.5), ylim = c(-7, 9), ratio = 7/16)
p = p + xlab("Med. Ctrl-Paired Act. log(ON:OFF)") + ylab("Prior Pfam Screen Act. log(ON:OFF)")
p = p + theme_linedraw() + theme(panel.grid = element_blank(), legend.position = c(0.775, 0.3))
p = p + guides(color = guide_legend(override.aes = list(size = 3)))
p = p + scale_color_manual(name = "", values = acr_colors)
ggsave("./nucpfam_correlation_act.pdf", p, width = 3.2, height = 3.2)
p

# pfam  act correlation ----------------------------------------------------------------------------

p = ggplot(baseline_priored_df,
           aes(x = med_d5, y = -avg_d5_pri, color = baseline_type))
# p = p + geom_vline(xintercept = -0.3, color = '#777777', linetype = "dashed")
p = p + geom_hline(yintercept = -1, color = '#777777', linetype = "dashed") # estimate from HT-recruit data
p = p + geom_point(size = 1.5)
p = p + geom_text(data = data.frame(x = 0.25, y = -6.9, label = "Pearson R = 0.85"),
          mapping = aes(x = x, y = y, label = label),
          inherit.aes = FALSE)
p = p + coord_fixed(xlim = c(-3.5, 3), ylim = c(-7, 1.5), ratio = 7/8.5)
p = p + xlab("Med. Ctrl-Paired Act. log(ON:OFF)") + ylab("Prior Pfam Screen Rep. log(ON:OFF)")
p = p + theme_linedraw() + theme(panel.grid = element_blank(), legend.position = c(0.2, 0.772),
    legend.text = element_text(size = 8),
    legend.key.height = unit(5.0, "mm"),
    # legend.box.background = element_rect(color="black"),
                                 legend.title = element_blank())
p = p + guides(color = guide_legend(override.aes = list(size = 3, alpha = 1.0)))
p = p + scale_color_manual(name = , values = acr_colors)
ggsave("./nucpfam_correlation_rep.pdf", p, width = 3.2, height = 3.2)
p


# violins -------------------------------------------------------------------------------------


ddf = df %>%
    dplyr::filter(type == "1 control Pair") %>%
    dplyr::mutate(hit = if_else(d1_Gene %in% c("DMD", NA), domain2, domain1),
                  hit_gene = if_else(d1_Gene %in% c("DMD", NA), d2_Gene, d1_Gene)) %>%
    dplyr::left_join(baseline_df,
                     by = c("hit" = "domain"),
                     suffix = c("", "_baseline")) %>%
    tidyr::pivot_longer(
        cols = c("avg_enrichment_d2", "avg_enrichment_d5"),
        names_to = "screen_day",
        names_prefix = "avg_enrichment_",
        values_to = "screen_score"
    )

make_plot = function(qdf) {
    ggplot(qdf) +
        facet_grid(rows = vars(baseline_type)) +
        geom_boxplot(aes(x = hit_gene, y = screen_score, fill = screen_day)) +
        geom_hline(yintercept = d5_baseline_threshold,
                   linetype = "dashed",
                   color = "#23BEDB") +
        geom_hline(yintercept = d2_baseline_threshold,
                   linetype = "dashed",
                   color = "#DBC60D") +
        scale_fill_manual(
            name = "",
            values = c("#DBC60D", "#23BEDB"),
            labels = c("Act.", "Repr.")
        ) +
        theme_linedraw() +
        theme(axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1
        )) +
        theme(
            panel.grid = element_blank(),
            # legend.position = c(0.775, 0.225),
            # legend.text = element_text(size = 8),
            # legend.key.height = unit(3, "mm"),
            # legend.key.width = unit(3, "mm"),
            # legend.background = element_blank(),
            # legend.key = element_blank(),
        ) +
        labs(x = "Gene Name", y = "log(ON:OFF)")
}

pa = make_plot(dplyr::filter(ddf, baseline_type == "Activator"))
pr = make_plot(dplyr::filter(ddf, baseline_type == "Repressor"))
pn = make_plot(dplyr::filter(ddf, baseline_type == "Non-hit"))
pd = make_plot(dplyr::filter(ddf, baseline_type == "Dual"))
px = pa / pr / pd / pn
ggsave(filename="./domain_boxplots.pdf", plot=px, height = 10, width = 10)
