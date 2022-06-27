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

df = read_csv("./repressors_fitted.csv")

p = ggplot(data = df, aes(x = med_d5, y = m))
p = p + geom_smooth(formula = y ~ x, method = "lm", size = 0.75,
                    se = FALSE, n = 100, color = "#e15759")
p = p + geom_point(size = 1, color = "#4e79a7")
p = p + geom_point(size = 0.005, color = "white")
p = p + geom_text(data = data.frame(x = -1.3, y = 0.975, label = bquote("Pearson R=0.89")), 
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE) 
p = p + coord_fixed(xlim = c(-3.5, 2.5), ylim=c(0, 1), ratio = 6)
p = p + scale_y_continuous(breaks = c(0, 0.5, 1))
p = p + scale_x_continuous(breaks = c(-3, 0, 3))
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Median Rep." ~ log[2](ON:OFF)))
p = p + ylab(bquote("Slope of Linear Fit"))
ggsave("./m_vs_med_d5.pdf", p, height = 2.5, width = 2.5)

p = ggplot(data = df, aes(x = med_d5, y = b))
p = p + geom_smooth(method = "lm", size = 0.75, se = FALSE, n = 100, color = "#e15759")
p = p + geom_point(size = 1, color = "#4e79a7")
p = p + geom_point(size = 0.005, color = "white")
p = p + geom_text(data = data.frame(x = -1.3, y = 0.975, label = bquote("Pearson R=0.94")), 
                  mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE) 
p = p + coord_fixed(xlim = c(-3.5, 2.5), ylim=c(-3, 1), ratio = 6/4)
p = p + theme_linedraw() + theme(panel.grid = element_blank())
p = p + xlab(bquote("Median Rep." ~ log[2](ON:OFF)))
p = p + ylab(bquote("Intercept of Linear Fit"))
ggsave("./b_vs_med_d5.pdf", p, height = 2.5, width = 2.5)
