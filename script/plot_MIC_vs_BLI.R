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
library(data.table)
library(ggrepel)
require(cowplot)

plot_BLI_vs_neut <- function(df, graphname){
  print (df)
  colorscale <- c(qualpal(n = 4, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  print (paste('correlation between BLI and MIC:', cor(-df$log_Kd_BLI, df$MIC)))
  textsize <- 7
  p <- ggplot(df,aes(x=-log_Kd_BLI, y=MIC, color=background)) +
         geom_point(alpha=0.7, pch=16, size=1.5) +
         geom_hline(yintercept = -0.3, linetype = "dashed", color = "black", size = 0.5) +
         scale_color_manual('',values=colorscale,drop=FALSE,
                            labels = c('Germline', 'Somatic')) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5,family="Arial"),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5,family="Arial"),
               axis.text=element_text(size=textsize,face="bold",family="Arial"),
               legend.key.size=unit(0.12,'in'),
               legend.spacing.x=unit(0, 'in'),
               legend.spacing.y=unit(0.3, 'in'),
               legend.title=element_text(size=textsize,face="bold",family="Arial"),
               legend.text=element_text(size=textsize,face="bold",family="Arial", margin = margin(r = 0, l = 0, t=-4, b=-4)),
               #legend.position='none') +
               legend.position=c(0.01,0.9)) +
         scale_y_continuous(limit=c(-1,6.5),breaks=c(0,1,2,3,4,5,6),labels=c('100','50','25','12.5','6.3','3.1','1.6')) +
         scale_x_continuous(limit=c(4.8,10.2),breaks=c(5,6,7,8,9,10),labels=c(expression(bold('10'^'4')),
                                                                              expression(bold('10'^'3')),
                                                                              expression(bold('10'^'2')),
                                                                              expression(bold('10'^'1')),
                                                                              expression(bold('10'^'0')),
                                                                              expression(bold('10'^'-1')))) +
         labs(x = bquote(bold(K[D]*" (nM) vs H1 stem")),
              y = expression(bold('MIC (μg/mL) vs H1N1 virus')))
  ggsave(graphname, p, height=2.3, width=2.3, dpi=300)
  }

df <- read_tsv('result/KD_compare_BLI_DMS.tsv') %>%
        mutate(exp = str_replace(exp, pattern = "GL_miniH1", replacement = "GL")) %>%
        mutate(exp = str_replace(exp, pattern = "WT_miniH1", replacement = "WT")) %>%
        rename(background=exp) %>%
        rename(mutation=mut) %>%
        select(background, mutation, log_Kd_BLI) %>%
        inner_join(read_tsv('data/CR9114_mut_neut.tsv'),by=c('background','mutation')) %>%
        mutate(MIC=log2(100/MIC))
plot_BLI_vs_neut(df, 'graph/CR9114_mut_neut_vs_BLI.png')
