#R code for Delta Thermal-shift Mos99 variants
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(data.table)
library(reshape)
library(stringr)
library(dplyr)
require(cowplot)

jitter <- position_jitter(width = 0.2, height = 0.1)

Ab_levels <- c('GL','WT')
df <- read_csv('data/CR9114_thermal_shift.csv') %>%
        mutate(Ab=factor(Ab,levels=Ab_levels))
colorscale  <- c(brewer.pal(8,"Accent"))
textsize <- 7
p <- ggplot() +
        geom_line(data=df, aes(x=Temperature,y=mean,color=Ab), linewidth=0.8) +
  scale_color_manual(values=colorscale, labels=c(bquote(bold("Germline CR9114" ~ '('*T[m]*' = 77.5°C)')),
                                                 bquote(bold("Somatic CR9114" ~ '('*T[m]*' = 78.5°C)')))) +
  theme_cowplot(12) +
  theme(plot.title=element_blank(),
        plot.background = element_rect(fill = "white"),
	axis.title.x=element_text(size=textsize,face="bold",hjust=0.5),
        axis.title.y=element_text(size=textsize,face="bold",vjust=0.5),
        axis.text.x=element_text(size=textsize,face="bold"),
	axis.text.y=element_text(size=textsize,face="bold"),
	legend.title=element_blank(),
        legend.text=element_text(size=textsize,face="bold", margin = margin(r = 0, l = 0, t=-4, b=-4)),
        legend.position=c(0.05,0.2))+
  labs(y=expression(bold('-d(RFU/dT)')),x=expression(bold('Temperature')))
ggsave('graph/thermal_shift.png',p,height=2,width=3, bg='white')
