colorblind_palette = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(reticulate)
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

csv_indices = str_pad(seq(1, 12), width=2, pad=0, "left")
csv_files = paste("csvs/", csv_indices, "_raw_counts.csv", sep="")
df = read_csv(csv_files, show_col_types = FALSE)
df$...1 <- NULL # drops an unused index column from the CSVs

screen_meta = list(
    list(2,  "activation", "bound",   1),
    list(2,  "activation", "unbound", 1),
    list(2,  "activation", "bound",   2),
    list(2,  "activation", "unbound", 2),
    list(5,  "repression", "bound",   1),
    list(5,  "repression", "unbound", 1),
    list(5,  "repression", "bound",   2),
    list(5,  "repression", "unbound", 2),
    list(10, "repression", "bound",   1),
    list(10, "repression", "unbound", 1),
    list(10, "repression", "bound",   2),
    list(10, "repression", "unbound", 2)
)

get_day = function(idx) {
    return(screen_meta[[strtoi(idx, base=10)]][[1]])
}
get_type = function(idx) {
    return(screen_meta[[strtoi(idx, base=10)]][[2]])
}
get_magsep = function(idx) {
    return(screen_meta[[strtoi(idx, base=10)]][[3]])
}
get_replicate = function(idx) {
    return(screen_meta[[strtoi(idx, base=10)]][[4]])
}

samplevec = unlist(map(df$sample, get_day))
typevec = unlist(map(df$sample, get_type))
magvec = unlist(map(df$sample, get_magsep))
repvec = unlist(map(df$sample, get_replicate))
df %>%
    dplyr::mutate(
        day = samplevec
    ) %>%
    dplyr::mutate(
        screen_type = typevec
    ) %>%
    dplyr::mutate(
        fraction = magvec
    ) %>%
    dplyr::mutate(
        replicate = repvec
    ) -> df

# some tiles have incorrect sequences 
wrong_labels = c(
    "Silencer_tiles;ENSG00000163602",
    "Silencer_tiles;ENSG00000174197",
    "Silencer_tiles;ENSG00000061273",
    "Silencer_tiles;ENSG00000188779",
    "Silencer_tiles;ENSG00000155090",
    "Silencer_tiles;ENSG00000141644",
    "DMD_control_tiles;ENSG00000198947;;31;" # drop the activator DMD tile
)

df %>%
    pivot_wider(
        id_cols = c("pair", "day", "screen_type", "replicate"),
        names_from = "fraction",
        values_from = "count"
    ) %>%
    dplyr::mutate(bound = replace_na(bound, 0)) %>%
    dplyr::mutate(unbound = replace_na(unbound, 0)) -> dfx
dfx = dfx[unlist(purrr::map(dfx$pair, function(x) {!any(str_detect(x, wrong_labels))})),]

write_csv(dfx, "./csvs/raw_counts.csv")


# counts filtering
dfx %>% 
    dplyr::filter(
        (bound >= 5 | unbound >= 5) & (bound + unbound >= 50)
    ) %>%
    dplyr:: mutate(
        bound = if_else(
            bound < 5,
            5,
            bound
        )
    ) %>%
    mutate(
        unbound = if_else(
            unbound < 5,
            5,
            unbound
        )
    ) -> dfx
dfx %>%
    group_by(day, screen_type, replicate) %>%
    mutate(norm_bound = bound/sum(bound)) %>%
    mutate(norm_unbound = unbound/sum(unbound)) -> dfx

# get normalized counts
dfx %>%
    mutate(enrichment_ratio = log2(norm_bound/norm_unbound)) %>%
    mutate(enrichment_fraction = norm_bound/(norm_bound + norm_unbound)) -> dfx
dfx %>%
    pivot_wider(
        id_cols = c("pair", "day"),
        names_from = "replicate",
        names_prefix = "r",
        values_from = c("enrichment_ratio", "enrichment_fraction",
                        "bound", "unbound",
                        "norm_bound", "norm_unbound")) %>%
    dplyr::mutate(
        day = as.factor(day)
    ) -> sdf
sdf %>%
    dplyr::mutate(
        type = case_when(
            str_count(pair, "control") == 0 ~ "0 control Pair",
            str_count(pair, "control") == 1 ~ "1 control Pair",
            str_count(pair, "control") == 2 ~ "2 control Pair"
        )
    ) %>%
    dplyr::arrange(day, type) %>%
    dplyr::ungroup() -> cdf
cdf %>%
    dplyr::mutate(avg_enrichment = (enrichment_ratio_r1 + enrichment_ratio_r2)/2) -> adf

# assembling the big dataframe
adf %>%
    pivot_wider(
        id_cols = c("pair", "type"),
        names_from = "day",
        names_prefix = "d",
        values_from = c("enrichment_ratio_r1", "enrichment_ratio_r2",
                        "enrichment_fraction_r1", "enrichment_fraction_r2",
                        "bound_r1", "bound_r2", "unbound_r1", "unbound_r2",
                        "norm_bound_r1", "norm_bound_r2", "norm_unbound_r1", "norm_unbound_r2",
                        "avg_enrichment")) -> xdf
xdf$domain1 = unlist(map(str_split(xdf$pair, ' --- '), first))
xdf$domain2 = unlist(map(str_split(xdf$pair, ' --- '), 
                         function (x) { return(nth(x, 2)) }))

# incorporate prior data
lib1 = read_csv('./csvs/base_oligo_library.csv')
lib1$...1 <- NULL
lib2 = copy(lib1)

names(lib1) <- paste0("d1_", names(lib1))
names(lib2) <- paste0("d2_", names(lib2))

xdf %>%
    dplyr::left_join(
        lib1,
        by = c("domain1" = "d1_label")
    ) %>%
    dplyr::left_join(
        lib2,
        by = c("domain2" = "d2_label")
    ) %>%
    dplyr::select(
        !ends_with("Extended Domain DNA sequence")
    ) -> ddf
write_csv(ddf, "./csvs/concatenation_filtered_scored_pairs.csv")
