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
require(cowplot)

plot_score_heatmap <- function(fitness_table, WTresibox, start_resi, end_resi, legend_title){
  textsize <- 5
  fitness_table <- fitness_table %>%
    filter(Pos >= start_resi & Pos <= end_resi)
  WTresibox     <- WTresibox %>%
    filter(Pos >= start_resi & Pos <= end_resi) %>%
    mutate(x=x-min(x)+1)
  p <-  ggplot() +
    geom_tile(data=fitness_table,aes(x=resi,y=aa,fill=parameter),color='black') +
    scale_fill_gradientn(colours=c("white","pink","red"),
                         limits=c(0,1000),
                         values=rescale(c(0,50,1000)),
                         #breaks=c(-1,0,1,2,3,4),
                         #labels=c('-1','0','1','2','3','4'),
                         guide="colorbar",
                         na.value="grey50") +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=textsize,face="bold",colour = 'black'),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
          axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
          axis.title=element_text(size=7,face="bold"),
          axis.line = element_line(colour = 'black', size = 0),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    guides(fill = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                 label.theme=element_text(size=7,face="bold",colour='black'),
                                 frame.colour="black",
                                 frame.linewidth = 0.5,
                                 ticks = TRUE,
                                 ticks.colour = "black",
                                 barwidth = 0.5, barheight = 6, title=legend_title)) +
    geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.2) +
    xlab("") +
    ylab("amino acid")
}

wrapper <- function(df_plot, graphname){
  df_plot <- df_plot %>%
    mutate(parameter=case_when(str_sub(resi,1,1)==aa ~ 1, TRUE ~ parameter)) %>% #Set WT parameter (usually 0 or 1)
    mutate(Mutation=paste(resi,aa,sep='')) %>%
    #select(Mutation, avg_ipt_freq, resi, Pos, aa, parameter) #Variable name for input freq may need to chance
    select(Mutation, resi, Pos, aa, parameter) #Variable name for input freq may need to chance
  print (df_plot)
  print (range(df_plot$parameter,na.rm=T))
  p1 <- plot_score_heatmap(df_plot, WTresibox, 1, 60, legend_title)
  p2 <- plot_score_heatmap(df_plot, WTresibox, 61, 120, legend_title)
  p <- grid.arrange(p1, p2, nrow=2)
  ggsave(graphname,p,width=4.5, height=4, dpi=300) #Output file may need to adjust
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','*'))

start_pos <- 1
end_pos   <- 120

df <- read_tsv('result/CR9114HC_WT_count_aa.tsv') %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  mutate(resi=str_sub(mut,1,-2))
residues <- unique(df$resi)

df <- df %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  mutate(resi=str_sub(mut,1,-2)) %>%
  mutate(aa=str_sub(mut,-1,-1)) %>%
  mutate(aa=ifelse(aa=='_', '*', aa)) %>%
  filter(aa %in% aa_level) %>%
  mutate(aa=factor(aa,levels=aa_level)) %>%
  mutate(resi=factor(resi,levels=unique(residues))) %>%
  complete(resi, aa) %>%
  mutate(Pos=str_sub(resi,2,-1)) %>%
  mutate(Pos=as.numeric(as.character(Pos))) %>%
  arrange(Pos)

WTresibox  <- df %>%
  select(resi,Pos) %>%
  unique() %>%
  mutate(WT_resi=str_sub(resi,1,1)) %>%
  mutate(x=seq(1,end_pos-start_pos+1)) %>%
  mutate(y=match(WT_resi,aa_level)) %>%
  select(resi,WT_resi,Pos,x, y)

legend_title <- "count"

df_plot <- df %>% mutate(parameter=input) 
wrapper(df_plot, 'graph/CR9114_WT_ipt_count_heatmap.png')
