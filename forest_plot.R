rm(list=ls())

## code adapted from https://www.khstats.com/blog/forest-plots/ (Katherine Hoffman)

## load up the packages we will need: 
library(tidyverse)
library(gt)
library(patchwork)
## ---------------------------
## load data
# load in results
res = read.csv("/Users/nikolayivanov/Desktop/Anesthesiology_copyFromCluster/SafavyniaLab/AF/res.csv")

## plotting
## ---------------------------
# create forest plot on log scale (middle section of figure)
p_mid <-
     res |>
     ggplot(aes(y = fct_rev(model))) +
     theme_classic() +
     geom_point(aes(x=estimate), shape=15, size=3) +
     geom_linerange(aes(xmin=conf.low, xmax=conf.high)) +
     labs(x="Odds Ratio") +
     coord_cartesian(ylim=c(1,5), xlim=c(1, 2.5))+
     geom_vline(xintercept = 1, linetype="dashed") +
     #annotate("text", x = 0.7, y = 5, label = "Protective") +
     #annotate("text", x = 1.3, y = 5, label = "Harmful") +
     theme(axis.line.y = element_blank(),
           axis.ticks.y= element_blank(),
           axis.text.y= element_blank(),
           axis.title.y= element_blank())

# wrangle results into pre-plotting table form
res_plot <- res |>
     # round estimates and 95% CIs to 2 decimal places for journal specifications
     mutate(across(
         c(estimate, conf.low, conf.high),
         ~ str_pad(
             round(.x, 2),
             width = 4,
             pad = "0",
             side = "right"
         )
     ),
     # add an "-" between OR estimate confidence intervals
     estimate_lab = paste0(estimate, " (", conf.low, "-", conf.high, ")")) |>
     # round p-values to two decimal places, except in cases where p < .001
     mutate(p.value = case_when(
         p.value < .001 ~ "<0.001",
         round(p.value, 2) == .05 ~ as.character(round(p.value,3)),
         p.value < .01 ~ str_pad( # if less than .01, go one more decimal place
             as.character(round(p.value, 3)),
             width = 4,
             pad = "0",
             side = "right"
         ),
         TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
             as.character(round(p.value, 2)),
             width = 4,
             pad = "0",
             side = "right"
         )
     )) |>
     # add a row of data that are actually column names which will be shown on the plot in the next step
     bind_rows(
         data.frame(
             model = "Model",
             estimate_lab = "Odds Ratio (95% CI)",
             conf.low = "",
             conf.high = "",
             p.value = "p-value"
         )
     ) |>
     mutate(model = fct_rev(fct_relevel(model, "Model")))

# left side of plot - hazard ratios
p_left <-
     res_plot  |>
     ggplot(aes(y = model)) + 
     geom_text(aes(x=0, label=model), hjust=0, fontface = "bold") +
     geom_text(aes(x=1, label=estimate_lab), hjust=0, fontface = ifelse(res_plot$estimate_lab == "Odds Ratio (95% CI)", "bold", "plain")) +
     theme_void() +
     coord_cartesian(xlim=c(0,4))

# right side of plot - pvalues
p_right <-
     res_plot  |>
     ggplot() +
     geom_text(aes(x=0, y=model, label=p.value), hjust=0, fontface = ifelse(res_plot$p.value == "p-value", "bold", "plain")) +
     theme_void() 
# layout design (top, left, bottom, right)
layout <- c(
     area(t = 0, l = 0, b = 30, r = 3),
     area(t = 1, l = 4, b = 30, r = 9),
     area(t = 0, l = 9, b = 30, r = 11)
 )
# final plot arrangement
p_left + p_mid + p_right + plot_layout(design = layout)

## save final figure
ggsave("/Users/nikolayivanov/Desktop/Anesthesiology_copyFromCluster/SafavyniaLab/AF/forest-plot.pdf", width=9, height=4)

#NAI

