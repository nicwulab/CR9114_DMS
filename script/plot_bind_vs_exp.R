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

plot_bind_vs_exp <- function(df, graphname){
  print (paste('correlation for:', graphname, cor(df$minus_delta_log_Kd, df$exp_score, use = "complete.obs")))
  textsize <- 7
  palette <- c(qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  p <- ggplot(filter(df, mut_class != 'NA'),aes(x=exp_score, y=minus_delta_log_Kd, color=mut_class)) +
         geom_point(alpha=0.5, pch=16, size=0.5) +
         scale_color_manual('',values=palette,drop=FALSE,labels=c('Missense','Nonsense','Silent')) + 
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold",margin = margin(l = -3, unit = "pt")),
               legend.position='none') +
         scale_y_continuous(breaks = integer_breaks) +
         labs(x = bquote(bold("Expression score")), y = bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM)"))) +
         xlim(-1,2)
   ggsave(graphname, p, height=1.5, width=1.5)
  }

read_df <- function(infile){
  aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))
  df <- read_csv(infile) %>%
    filter(ipt_freq >= 0.0001) %>%
    mutate(Mutation=mut) %>%
    filter(Mutation != 'WT') %>%
    mutate(resi=str_sub(Mutation,1,-2)) %>%
    mutate(aa=str_sub(Mutation,-1,-1)) %>%
    filter(aa %in% aa_level) %>%
    mutate(aa=factor(aa,levels=aa_level)) %>%
    complete(resi, aa) %>%
    mutate(Pos=str_sub(resi,2,-1)) %>%
    mutate(Pos=as.numeric(as.character(Pos))) %>%
    arrange(Pos) %>%
    mutate(resi=factor(resi,levels=unique(resi))) %>%
    mutate(Mutation=paste(resi,aa,sep='')) %>%
    mutate(var = paste(Pos,aa,sep=''))
  return (df)
  }


df_exp  <- read_csv('result/HC_WT_expression_score_kabat.csv') %>%
             select(mut, exp_score)

df_GL   <- read_df('result/HC_GL_KD_table_summary_kabat.csv')
df_H1   <- read_df('result/HC_WT_H1_KD_table_summary_kabat.csv') %>%
             inner_join(df_exp, by='mut')
df_H3   <- read_df('result/HC_WT_H3_KD_table_summary_kabat.csv') %>%
             inner_join(df_exp, by='mut')
df_fluB <- read_df('result/HC_WT_fluB_KD_table_summary_kabat.csv') %>%
             inner_join(df_exp, by='mut')

plot_bind_vs_exp(df_GL, 'graph/bind_vs_exp_GL_H1.png')
plot_bind_vs_exp(df_H1, 'graph/bind_vs_exp_WT_H1.png')
plot_bind_vs_exp(df_H3, 'graph/bind_vs_exp_WT_H3.png')
plot_bind_vs_exp(df_fluB, 'graph/bind_vs_exp_WT_fluB.png')
