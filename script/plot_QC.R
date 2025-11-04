#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(stringr)
library(qualpalr)
require(cowplot)

integer_breaks <- function(limits) {
  # limits is a vector with min and max of the axis
  floor(limits[1]):ceiling(limits[2])
  }

plot_replicate_cor <- function(filename, title, graphname){
  df <- read_csv(filename) %>%
          filter(mut_class!='WT')
  print (filename)
  print (cor(df$minus_delta_log_Kd_rep1, df$minus_delta_log_Kd_rep2, use = "complete.obs"))
  palette <- c(qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  textsize <- 7
  p <- ggplot(df,aes(x=minus_delta_log_Kd_rep1, y=minus_delta_log_Kd_rep2, color=mut_class)) +
         geom_point(alpha=0.5, pch=16, size=0.5) +
         scale_color_manual('',values=palette,drop=FALSE) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         scale_y_continuous(breaks = integer_breaks) +
         labs(x=expression(bold('Replicate 1')),y=expression(bold('Replicate 2')), title = title)
  ggsave(graphname, p, height=1.5, width=1.5)
  }

plot_replicate_cor('result/HC_GL_KD_table_summary_kabat.csv', 'GL vs H1', 'graph/rep_cor_GL_H1.png')
plot_replicate_cor('result/HC_WT_H1_KD_table_summary_kabat.csv', 'WT vs H1', 'graph/rep_cor_WT_H1.png')
plot_replicate_cor('result/HC_WT_H3_KD_table_summary_kabat.csv', 'WT vs H3', 'graph/rep_cor_WT_H3.png')
plot_replicate_cor('result/HC_WT_fluB_KD_table_summary_kabat.csv', 'WT vs fluB', 'graph/rep_cor_WT_fluB.png')
