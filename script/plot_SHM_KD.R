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
library(ggrepel)
require(cowplot)

plot_SHM_KD <- function(df, ab_list, start, end){
  textsize <- 7
  ab_list_tmp <- ab_list[start:end]
  df <- df %>%
          filter(ab_name %in% ab_list_tmp) %>%
          mutate(ab_name=factor(ab_name, levels=ab_list_tmp)) %>%
          arrange(desc(in_CR9114))
  p <- ggplot() +
         geom_point(data=df, aes(x=ab_name, y=`minus_delta_log_KD (GL)`, color=in_CR9114), alpha=0.7, pch=16, size=1) +
         scale_color_manual('',values=c('red','grey20'),drop=FALSE) +
         geom_hline(yintercept = -1, color = "black", linetype = "dotted", size = 0.5) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text.x=element_text(size=textsize-1,face="bold",angle=45, hjust=1, vjust=1),
               axis.text.y=element_text(size=textsize,face="bold"),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         scale_y_continuous(limit = c(-3.7,0.6)) +
         guides(color = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                      label.theme=element_text(size=7,face="bold",colour='black'),
                                      frame.colour="black",
                                      frame.linewidth = 0.5,
                                      ticks = TRUE,
                                      ticks.colour = "black",
                                      barwidth = 10, barheight = 0.5, title=bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM)")))) +
         labs(x = bquote(bold("")),
              y = bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM)")))
  return (p)
  }

plot_SHM_freq_vs_KD <- function(df, output_graph){
  textsize <- 8
  p <- ggplot() +
         geom_point(data=df, aes(x=freq, y=`minus_delta_log_KD (GL)`, color=in_CR9114), alpha=0.7, pch=16, size=1) +
         geom_label_repel(data=filter(df, (freq>0.3) | (freq > 0.07 & `minus_delta_log_KD (GL)` < -1)),
                          aes(x=freq, y=`minus_delta_log_KD (GL)`,label=SHM),
                          color="black", min.segment.length=0, segment.size=0.2, size=2.5, force=15, force_pull=1,
                          seed=5, max.overlaps = Inf, label.size=0, label.padding=0.05) +
         scale_color_manual('',values=c('red','grey20'),drop=FALSE) +
         geom_hline(yintercept = -1, color = "black", linetype = "dotted", size = 0.5) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text.x=element_text(size=textsize-1,face="bold",angle=0, hjust=0.5, vjust=1),
               axis.text.y=element_text(size=textsize,face="bold"),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         scale_y_continuous(limit = c(-3.7,0.6)) +
         scale_x_continuous(breaks = c(0.0,0.1,0.2,0.3,0.4,0.5), labels=c('0%','10%','20%','30%','40%','50%')) +
         guides(color = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                      label.theme=element_text(size=7,face="bold",colour='black'),
                                      frame.colour="black",
                                      frame.linewidth = 0.5,
                                      ticks = TRUE,
                                      ticks.colour = "black",
                                      barwidth = 10, barheight = 0.5, title=bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM)")))) +
         labs(x = bquote(bold("Frequency among IGHV1-69 HA stem antibodies")),
              y = bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM)")))
  ggsave('graph/Ab_SHM_KD_vs_freq.png', p, height=3, width=3, dpi=300)
  }

df <- read_tsv('result/IGHV1-69_HAstem_Ab_SHM_KD.tsv') %>%
        mutate(in_CR9114=factor(in_CR9114, levels=c('yes','no'))) 
print (range(df$`minus_delta_log_KD (GL)`))
ab_list <- df %>%
             group_by(ab_name) %>%
             summarise(mean = min(`minus_delta_log_KD (GL)`)) %>%
             arrange(mean) %>%
             .$ab_name %>%
             rev()
p1 <- plot_SHM_KD(df, ab_list, 1, 40)
p2 <- plot_SHM_KD(df, ab_list, 41, 80)
p3 <- plot_SHM_KD(df, ab_list, 81, 120)
p4 <- plot_SHM_KD(df, ab_list, 121, 160)
ggsave('graph/Ab_SHM_KD_1.png', p1, height=1.5, width=6.5, dpi=300)
ggsave('graph/Ab_SHM_KD_2.png', p2, height=1.44, width=6.5, dpi=300)
ggsave('graph/Ab_SHM_KD_3.png', p3, height=1.35, width=6.5, dpi=300)
ggsave('graph/Ab_SHM_KD_4.png', p4, height=1.52, width=6.5, dpi=300)

ab_count <- length(unique(df$ab_name))
df_SHM <- select(df, SHM, in_CR9114, `minus_delta_log_KD (GL)`) %>%
            group_by(SHM, in_CR9114, `minus_delta_log_KD (GL)`) %>%
            summarize(count = n()) %>%
            mutate(freq=count/ab_count)
plot_SHM_freq_vs_KD(df_SHM, "graph/IGHV169_SHM_vs_DMS.png")
