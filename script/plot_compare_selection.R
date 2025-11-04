#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggtext)
library(gridExtra)
library(stringr)
library(qualpalr)
require(cowplot)

integer_breaks <- function(limits) {
  # limits is a vector with min and max of the axis
  floor(limits[1]):ceiling(limits[2])
  }

plot_selection_cor <- function(df, graphname, xlab, ylab){
  colorscale  <- c(brewer.pal(9,"Set1"))
  df <- filter(df, mut_class=='missense')
  print (graphname)
  print (paste('correlation for (FR):', graphname, cor(filter(df, color=='FR')$x, filter(df, color=='FR')$y,
         use = "complete.obs", method="pearson")))
  print (paste('correlation for (CDR):', graphname, cor(filter(df, color=='CDR')$x, filter(df, color=='CDR')$y,
         use = "complete.obs", method="pearson")))
  textsize <- 7
  p <- ggplot(df,aes(x=x, y=y, color=color)) +
         geom_point(alpha=0.5, pch=16, size=0.5) +
         scale_color_manual('',values=colorscale,drop=FALSE) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text=element_text(size=textsize,face="bold"),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         scale_y_continuous(breaks = integer_breaks) +
         labs(x = bquote(bold(.(xlab))),
              y = bquote(bold(.(ylab))),
              title = bquote(bold(-"Δ"*log["10"]*' '*K[D]*" (nM)")))
   ggsave(graphname, p, height=1.5, width=1.5, dpi=300)
  }

coloring <- function(region){
  if (grepl('FR',region)){return ("FR")}
  else {return ('CDR')}
  }

read_df <- function(infile, nt_dist_min, nt_dist_max){
  aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))
  df <- read_csv(infile) %>%
    filter(ipt_freq >= 0.0001) %>%
    rename(Mutation=mut) %>%
    filter(Mutation != 'WT') %>%
    filter(codon_dist >= nt_dist_min) %>%
    filter(codon_dist <= nt_dist_max) %>%
    mutate(resi=str_sub(Mutation,1,-2)) %>%
    mutate(aa=str_sub(Mutation,-1,-1)) %>%
    filter(aa %in% aa_level) %>%
    mutate(aa=factor(aa,levels=aa_level)) %>%
    complete(resi, aa) %>%
    mutate(Pos=str_sub(resi,2,-1)) %>%
    mutate(Pos=as.numeric(as.character(Pos))) %>%
    arrange(Pos) %>%
    mutate(resi=factor(resi,levels=unique(resi))) %>%
    mutate(minus_delta_log_Kd=case_when(str_sub(resi,1,1)==aa ~ 0, TRUE ~ minus_delta_log_Kd)) %>%
    mutate(Mutation=paste(resi,aa,sep='')) %>%
    mutate(var = paste(Pos,aa,sep='')) %>%
    mutate(color=mapply(coloring, region)) %>%
    filter(ipt_freq >= 0.0001) %>%
    filter(mut_class=='missense') %>%
    select(Mutation, mut_class, mut_kabat, region, color, resi, Pos, aa, var, minus_delta_log_Kd, ipt_freq)
  return (df)
  }

# Combine region assignments

wrapper <- function(nt_dist_min, nt_dist_max, suffix){
  df_GL   <- read_df('result/HC_GL_KD_table_summary_kabat.csv', nt_dist_min, nt_dist_max) %>%
		rename(minus_delta_log_Kd_GL = minus_delta_log_Kd, ipt_freq_GL = ipt_freq)
  df_H1   <- read_df('result/HC_WT_H1_KD_table_summary_kabat.csv', nt_dist_min, nt_dist_max) %>%
		rename(minus_delta_log_Kd_H1 = minus_delta_log_Kd, ipt_freq_H1 = ipt_freq) 
  df_H3   <- read_df('result/HC_WT_H3_KD_table_summary_kabat.csv', nt_dist_min, nt_dist_max) %>%
		rename(minus_delta_log_Kd_H3 = minus_delta_log_Kd, ipt_freq_H3 = ipt_freq) 
  df_fluB <- read_df('result/HC_WT_fluB_KD_table_summary_kabat.csv', nt_dist_min, nt_dist_max) %>%
		rename(minus_delta_log_Kd_fluB = minus_delta_log_Kd, ipt_freq_fluB = ipt_freq) 

  df <- inner_join(select(df_H1, mut_class, region, color, var, minus_delta_log_Kd_H1), 
		   select(df_GL, var, minus_delta_log_Kd_GL), by='var') %>%
	  rename(x=minus_delta_log_Kd_GL, y=minus_delta_log_Kd_H1)
  plot_selection_cor(df, paste0('graph/Kd_cor_GL_H1',suffix,'.png'), 'Germline vs H1', 'Somatic vs H1')

  df <- inner_join(select(df_H3, mut_class, region, color, var, minus_delta_log_Kd_H3), 
		   select(df_GL, var, minus_delta_log_Kd_GL), by='var') %>%
	  rename(x=minus_delta_log_Kd_GL, y=minus_delta_log_Kd_H3)
  plot_selection_cor(df, paste0('graph/Kd_cor_GL_H3',suffix,'.png'), 'Germline vs H1', 'Somatic vs H3')

  df <- inner_join(select(df_fluB, mut_class, region, color, var, minus_delta_log_Kd_fluB), 
		   select(df_GL, var, minus_delta_log_Kd_GL), by='var') %>%
	  rename(x=minus_delta_log_Kd_GL, y=minus_delta_log_Kd_fluB)
  plot_selection_cor(df, paste0('graph/Kd_cor_GL_fluB',suffix,'.png'), 'Germline vs H1', 'Somatic vs BHA')

  df <- inner_join(select(df_H1, mut_class, region, color, var, minus_delta_log_Kd_H1),
		   select(df_H3, var, minus_delta_log_Kd_H3), by='var') %>%
	  rename(x=minus_delta_log_Kd_H1, y=minus_delta_log_Kd_H3)
  plot_selection_cor(df, paste0('graph/Kd_cor_H1_H3',suffix,'.png'), 'Somatic vs H1', 'Somatic vs H3')

  df <- inner_join(select(df_H1, mut_class, region, color, var, minus_delta_log_Kd_H1),
		   select(df_fluB, var, minus_delta_log_Kd_fluB), by='var') %>%
	  rename(x=minus_delta_log_Kd_H1, y=minus_delta_log_Kd_fluB)
  plot_selection_cor(df, paste0('graph/Kd_cor_H1_fluB',suffix,'.png'), 'Somatic vs H1', 'Somatic vs BHA')

  df <- inner_join(select(df_H3, mut_class, region, color, var, minus_delta_log_Kd_H3),
		   select(df_fluB, var, minus_delta_log_Kd_fluB), by='var') %>%
	  rename(x=minus_delta_log_Kd_H3, y=minus_delta_log_Kd_fluB)
  plot_selection_cor(df, paste0('graph/Kd_cor_H3_fluB',suffix,'.png'), 'Somatic vs H3', 'Somatic vs BHA')
  }

wrapper(0,99, '_all')
wrapper(1,1, '_nt1')
