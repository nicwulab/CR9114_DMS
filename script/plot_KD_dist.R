#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(stringr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(ggbeeswarm)
require(cowplot)
library(grid)

integer_breaks <- function(limits) {
  floor(limits[1]):ceiling(limits[2])
  }

plot_KD_dist_by_exp <- function(df, graphname){
  colorscale <- c(qualpal(n = 4, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  textsize <- 7
  p <- ggplot(df,aes(x=exp, y=minus_delta_log_Kd, fill=exp)) +
         geom_violin(linewidth=0.5,width=1.2) +
         geom_boxplot(width=0.05,outlier.shape = NA) +
         scale_fill_manual('',values=colorscale,drop=FALSE) +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text=element_text(size=textsize,face="bold",colour = 'black'),
               axis.text.x=element_text(angle=45,hjust=1,vjust=1.05,colour = 'black'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         #scale_y_continuous(breaks = integer_breaks) +
         scale_x_discrete(labels = c('Germline vs H1', 'Somatic vs H1', 'Somatic vs H3', 'Somatic vs BHA')) +
         labs(x = bquote(bold("")),
              y = bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM)")),
              title = bquote(bold("Missense mutations")))
   ggsave(graphname, p, height=2.5, width=4, dpi=300)
  }

plot_KD_dist_by_codon <- function(df, graphname){
  colorscale <- c(qualpal(n = 4, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  textsize <- 7
  p <- ggplot(df,aes(x=exp, y=minus_delta_log_Kd, fill=codon_dist)) +
         geom_violin(linewidth=0.5,position=position_dodge(0.7)) +
         geom_boxplot(width=0.05,outlier.shape = NA,position=position_dodge(0.7)) +
         scale_fill_manual('',values=colorscale,drop=FALSE) +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text=element_text(size=textsize,face="bold",colour = 'black'),
               axis.text.x=element_text(angle=45,hjust=1,vjust=1.05,colour = 'black'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='top') +
         #scale_y_continuous(breaks = integer_breaks) +
         scale_x_discrete(labels = c('Germline vs H1', 'Somatic vs H1', 'Somatic vs H3', 'Somatic vs BHA')) +
         labs(x = bquote(bold("")),
              y = bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM)")),
              title = bquote(bold("Missense mutations by nucleotide distance to WT")))
   ggsave(graphname, p, height=2.5, width=4, dpi=300)
  }

read_file <- function(infile, exp){
  df <- read_csv(infile) %>%
    filter(mut_class=='missense') %>%
    mutate(exp=exp) %>%
    filter(ipt_freq >= 0.0001) %>%
    mutate(codon_dist=as.character(codon_dist)) %>%
    select(exp, codon_dist, mut_class, minus_delta_log_Kd) %>%
    filter(!is.na(minus_delta_log_Kd)) %>%
    mutate(minus_delta_log_Kd = pmax(minus_delta_log_Kd, -2))
  return (df)
  }
  
# Run the wrapper function (uncomment and adjust paths as needed)
df1 <- read_file('result/HC_GL_KD_table_summary_kabat.csv', 'GL_H1')
df2 <- read_file('result/HC_WT_H1_KD_table_summary_kabat.csv', 'WT_H1')
df3 <- read_file('result/HC_WT_H3_KD_table_summary_kabat.csv', 'WT_H3')
df4 <- read_file('result/HC_WT_fluB_KD_table_summary_kabat.csv', 'WT_fluB')
df <- rbind(df1, df2, df3, df4) %>%
         mutate(exp=factor(exp,levels=c('GL_H1','WT_H1','WT_H3','WT_fluB')))
plot_KD_dist_by_exp(df, 'graph/KD_distribution_by_exp.png') 
print (paste0("median for GL vs H1: ", median(df1$minus_delta_log_Kd)))
print (paste0("median for WT vs H1: ", median(df2$minus_delta_log_Kd)))
print (paste0("median for WT vs H3: ", median(df3$minus_delta_log_Kd)))
print (paste0("median for WT vs B: ", median(df4$minus_delta_log_Kd)))
plot_KD_dist_by_codon(df, 'graph/KD_distribution_by_codon.png') 
