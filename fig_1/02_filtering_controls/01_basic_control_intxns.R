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

# Computing an initial, unfiltered set of singlet scores
df = read_csv('../01_raw_counts/csvs/concatenation_filtered_scored_pairs.csv')
oligos = read_csv('../01_raw_counts/csvs/base_oligo_library.csv')

get_df = function (o) { 
    df %>%
        dplyr::filter(domain1 == o | domain2 == o) %>%
        dplyr::filter(type == "1 control Pair") %>%
        dplyr::summarise(
            domain = o,
            avg_d2 = mean(avg_enrichment_d2, na.rm = TRUE),
            sd_d2 = sd(avg_enrichment_d2, na.rm = TRUE),
            avg_d5 = mean(avg_enrichment_d5, na.rm = TRUE),
            sd_d5 = sd(avg_enrichment_d5, na.rm = TRUE),
            avg_d10 = mean(avg_enrichment_d10, na.rm = TRUE),
            sd_d10 = sd(avg_enrichment_d10, na.rm = TRUE),
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

b1df = copy(baseline_df)
b2df = copy(baseline_df)
names(b1df) = paste0("d1_", names(b1df))
names(b2df) = paste0("d2_", names(b2df))

df %>%
    dplyr::left_join(
        b1df,
        by = c("domain1" = "d1_domain")
    ) %>%
    dplyr::left_join(
        b2df,
        by = c("domain2" = "d2_domain")
    ) -> fdf
fdf %>% # add compositions
    dplyr::mutate(
        composition = case_when(
            str_detect(d1_Description, "Control") & str_detect(d2_Description, "Control") ~ "C-C",
            str_detect(d1_Description, "Control") & str_detect(d2_Description, "Rep") ~ "C-R",
            str_detect(d1_Description, "Control") & str_detect(d2_Description, "Act") ~ "C-A",
            str_detect(d1_Description, "Rep") & str_detect(d2_Description, "Control") ~ "R-C",
            str_detect(d1_Description, "Rep") & str_detect(d2_Description, "Rep") ~ "R-R",
            str_detect(d1_Description, "Rep") & str_detect(d2_Description, "Act") ~ "R-A",
            str_detect(d1_Description, "Act") & str_detect(d2_Description, "Control") ~ "A-C",
            str_detect(d1_Description, "Act") & str_detect(d2_Description, "Rep") ~ "A-R",
            str_detect(d1_Description, "Act") & str_detect(d2_Description, "Act") ~ "A-A",
        )
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
            str_detect(composition, "A-A") ~ "Activator"
        )
    ) -> fdf
fdf$character = factor(fdf$character, levels = c("Control", "Activator",
                                                 "Repressor", "Both"))

# Generating an inital dataframe of control-domain interactions
onectrls = dplyr::filter(fdf, type == "1 control Pair")
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
    fdf %>%
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
        ),
        hit_score_d2 = if_else(
            str_detect(domain1, "control"),
            d2_avg_d2,
            d1_avg_d2
        ),
        hit_score_d5 = if_else(
            str_detect(domain1, "control"),
            d2_avg_d5,
            d1_avg_d5
        ),
        hit_score_d10 = if_else(
            str_detect(domain1, "control"),
            d2_avg_d10,
            d1_avg_d10
        ),
    ) -> cart_df

cart_df$hit_gene = as.factor(cart_df$hit_gene)
write_csv(cart_df, "./01_cartesian_onecontrol_df.csv")

# Generate an initial DF with an order to the domains
df = read_csv("./01_cartesian_onecontrol_df.csv")
df %>%
    dplyr::select(c(hit_gene, control, ctrl_character, order,
                    hit_score_d2, hit_score_d5, hit_score_d10,
                    avg_enrichment_d2, avg_enrichment_d5, avg_enrichment_d10)) %>%
    tidyr::pivot_longer(
        cols = c(hit_score_d2, hit_score_d5, hit_score_d10),
        names_to = "day",
        names_prefix = "hit_score_d",
        values_to = "hit_score"
    ) %>%
    tidyr::pivot_longer(
        cols = c(avg_enrichment_d2, avg_enrichment_d5, avg_enrichment_d10),
        names_to = "day_",
        names_prefix = "avg_enrichment_d",
        values_to = "avg_pair_enrichment"
    ) %>%
    dplyr::filter(
        day == day_
    ) %>%
    dplyr::select(-day_) %>%
    dplyr::distinct() -> order_df
write_csv(order_df, "./01_one_control_ordered_df.csv")