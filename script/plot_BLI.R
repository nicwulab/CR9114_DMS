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

integer_breaks <- function(limits) {
  floor(limits[1]):ceiling(limits[2])
  }

plot_BLI_vs_DMS <- function(df, graphname){
  print (df)
  df <- df %>%
         filter(!is.na(minus_delta_log_Kd_DMS)) 
  colorscale <- c(qualpal(n = 4, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  #colorscale <- c('black','black','black','black')
  print (paste('correlation between BLI and DMS:', cor(df$minus_delta_log_Kd_DMS, df$minus_delta_log_Kd_BLI)))
  textsize <- 7
  p <- ggplot(df,aes(x=minus_delta_log_Kd_BLI, y=minus_delta_log_Kd_DMS, color=exp)) +
         #geom_point(alpha=0.7, pch=16, size=1.5) +
         geom_point(alpha=0.7, pch=16, size=1.5) +
         scale_color_manual('',values=colorscale,drop=FALSE,
                            labels = c('Germline vs H1', 'Somatic vs H1', 'Somatic vs H3', 'Somatic vs BHA')) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.12,'in'),
               legend.spacing.x=unit(0, 'in'),
               legend.spacing.y=unit(0.3, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold", margin = margin(r = 0, l = 0, t=-4, b=-4)),
               legend.position='none') +
               #legend.position=c(0.5,0.2)) +
         scale_y_continuous(breaks = integer_breaks) +
         labs(x = bquote(bold(-"Δ"*log["10"]*' '*K[D]*" (nM) from BLI")),
              y = bquote(bold(-"Δ"*log["10"]*' '*K["D,app"]*" (nM) from Tite-Seq")))
  ggsave(graphname, p, height=2.3, width=2.3, dpi=300)
  #ggsave(graphname, p, height=1.5, width=1.5, dpi=300)
  }

plot_BLI_vs_DMS_raw <- function(df, graphname){
  print (df)
  df <- df %>%
         filter(!is.na(log_Kd_DMS))
  colorscale <- c(qualpal(n = 4, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex)
  print (paste('correlation between BLI and DMS (raw):', cor(df$log_Kd_DMS, df$log_Kd_BLI)))
  textsize <- 9
  p <- ggplot(df,aes(x=log_Kd_BLI, y=log_Kd_DMS, color=exp)) +
         geom_point(alpha=0.7, pch=16, size=1.5) +
         geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
         scale_color_manual('',values=colorscale,drop=FALSE,
                            labels = c('Germline vs H1', 'Somatic vs H1', 'Somatic vs H3', 'Somatic vs BHA')) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.12,'in'),
               legend.spacing.x=unit(0, 'in'),
               legend.spacing.y=unit(0.3, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold", margin = margin(r = 0, l = 0, t=-4, b=-4)),
               #legend.position='none') +
               legend.position='right') +
         scale_y_continuous(breaks = integer_breaks) +
         labs(x = bquote(bold(log["10"]*' '*K[D]*" (nM) from BLI")),
              y = bquote(bold(log["10"]*' '*K["D,app"]*" (nM) from Tite-Seq")))
  ggsave(graphname, p, height=2.5, width=3.7, dpi=300)
  }

plot_BLI_heatmap <- function(df, graphname){
  textsize <- 7
  p <-  ggplot() +
          geom_tile(data=df,aes(x=variable,y=mut,fill=value)) +
          scale_fill_gradientn(colours=c("purple","white","white","orange"),
                limits=c(-4,2),
                values=rescale(c(-4, -0.1, 0.1, 2)),
                #breaks=c(-2,-1,0,1),
                #labels=c('-2','-1','0','1'),
                guide="colorbar",
                na.value="grey") +
          scale_x_discrete(labels = c('Germline vs H1', 'Somatic vs H1', 'Somatic vs H3', 'Somatic vs BHA')) +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=45,hjust=1,vjust=1.05,colour = 'black'),
                axis.text.y=element_text(hjust=1,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=7,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1)) +
          guides(fill = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                       label.theme=element_text(size=7,face="bold",colour='black'),
                                       frame.colour="black",
                                       frame.linewidth = 0.5,
                                       ticks = TRUE,
                                       ticks.colour = "black",
                                       barwidth = 0.7, barheight = 5, title=bquote(bold(-"Δ"*log["10"]*' '*K[D]*" (nM)")))) +
          xlab("") +
          ylab("Mutation")
  ggsave(graphname, p, height=2.5, width=2.4, dpi=300)
  }

df <- read_tsv('result/KD_compare_BLI_DMS.tsv') %>%
        filter(mut != 'WT') %>%
        mutate(exp=factor(exp,levels=c('GL_miniH1','WT_miniH1','WT_H3','WT_fluB'))) %>%
        data.table()
plot_BLI_vs_DMS(df, 'graph/KD_compare_BLI_DMS_delta.png')
plot_BLI_vs_DMS_raw(df, 'graph/KD_compare_BLI_DMS_raw.png')

df[mut == 'R83E', mut := 'R/T83E']
df[mut == 'T83E', mut := 'R/T83E']
df[mut == 'K73W', mut := 'K/I73W']
df[mut == 'I73W', mut := 'K/I73W']
df[mut == 'I52V', mut := 'I/S52V']
df[mut == 'S52V', mut := 'I/S52V']
df[mut == 'I52G', mut := 'I/S52G']
df[mut == 'S52G', mut := 'I/S52G']

print (df)
df <- df %>%
        select(mut, exp, minus_delta_log_Kd_BLI) %>%
        pivot_wider(names_from=exp, values_from=minus_delta_log_Kd_BLI) %>% 
        data.table() %>%
        melt(id='mut') %>%
        mutate(variable=factor(variable,levels=c('GL_miniH1','WT_miniH1','WT_H3','WT_fluB'))) %>%
        mutate(mut=factor(mut,levels=rev(c("S24F","S24L","D46E","I/S52G","I/S52V","F54H","K/I73W","R/T83E","V102L","I73W/F74S"))))
plot_BLI_heatmap(df, 'graph/BLI_heatmap.png')
