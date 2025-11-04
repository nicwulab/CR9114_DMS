#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(stringr)
library(dplyr)
library(gridExtra)
require(cowplot)
library(grid)

extract_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[legend_index]]
}

plot_freq_heatmap <- function(df, WTresibox, SHMresibox, start_resi, end_resi, xlabels, highlight_list) {

  textsize <- 5
  df <- df %>%
    filter(Pos >= start_resi & Pos <= end_resi)
  
  WTresibox <- WTresibox %>%
    filter(Pos >= start_resi & Pos <= end_resi) %>%
    mutate(x = x - min(x) + 1)
  
  SHMresibox <- SHMresibox %>%
    filter(x >= start_resi & x <= end_resi) %>%
    mutate(x = x - start_resi + 1)
  
  # Combine highlights using the shared list with labels
  highlight_tiles <- df %>%
    inner_join(highlight_list, by = c("Pos", "aa")) %>%
    mutate(aa = factor(aa, levels = levels(df$aa)),
           y = as.numeric(aa))
  
  p <- ggplot() +
    geom_tile(data = df, aes(x = x, y = aa, fill = ipt_freq)) +
    
    # Green highlight outlines
    scale_fill_gradientn(colours = c("white", "red", "red"),
                         limits = c(0, 0.0022),
                         values = rescale(c(0,0.0011,0.0022)),
                         breaks = c(0,0.001,0.002),
                         labels = c('0%','0.1%','0.2%'),
                         guide = "colorbar",
                         na.value = "white") +
    
    guides(
      color = guide_legend(
        override.aes = list(
          size = 0.5,
          fill = NA
        ),
        label.theme = element_text(size = 6, face = "bold", colour = "black"),
        title.theme = element_text(size = 7, face = "bold", colour = "black"),
        keywidth = unit(0.15, "cm"),
        keyheight = unit(0.15, "cm")
      ),
      fill = guide_colorbar(
        title = bquote(bold('')),
        title.theme = element_text(size = 6, face = "bold", colour = 'black', hjust = 0.5),
        label.theme = element_text(size = 6, face = "bold", colour = 'black', hjust = 0.5),
        frame.colour = "black",
        frame.linewidth = 0.5,
        ticks = TRUE,
        ticks.colour = "black",
        barwidth = unit(3.0, "inches"),
        barheight = unit(0.05, "inches"),
        direction = "horizontal"
      )
    ) +
    
    theme_void() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = textsize, face = "bold", colour = 'black'),
      axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = textsize, face = "bold", colour = 'black'),
      axis.title.y = element_text(size = 7, face = "bold", margin = margin(r = 6), angle = 90),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.1, "cm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      legend.position = "right",
      legend.box = "horizontal",
      legend.text.align = 0,
      legend.spacing.y = unit(0.0, "cm"),
      legend.margin = margin(t = -20, r = 0, b = 0, l = 0)
    ) +
    
    geom_point(data = WTresibox, aes(x = x, y = y), color = 'black', size = 0.2) +
    scale_x_continuous(breaks = 1:length(xlabels), labels = xlabels, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    xlab("") +
    ylab("amino acid")
  
  return(p)
}

wrqpper <- function(infile, graphname, c1, c2, c3, c4, c5){
  aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','*'))
  
  df <- read_csv(infile) %>%
    rename(Mutation = mut) %>%
    filter(Mutation != 'WT') %>%
    mutate(resi = str_sub(Mutation, 1, -2),
           aa = str_sub(Mutation, -1, -1),
           aa = ifelse(aa == "_", "*", aa)) %>%
    filter(aa %in% aa_level) %>%
    mutate(aa = factor(aa, levels = aa_level)) %>%
    complete(resi, aa) %>%
    ungroup() %>%
    mutate(Pos = as.numeric(str_sub(resi, 2, -1))) %>%
    arrange(Pos) %>%
    mutate(resi = factor(resi, levels = unique(resi)),
           minus_delta_log_Kd = case_when(str_sub(resi, 1, 1) == aa ~ 0, TRUE ~ minus_delta_log_Kd),
           minus_delta_log_Kd = pmin(minus_delta_log_Kd, 0.5),
           minus_delta_log_Kd = pmax(minus_delta_log_Kd, -2),
           Mutation = paste(resi, aa, sep = ''))
  
  # DEBUG: print available columns
  print(colnames(df))
  
  # Use only the columns that exist in df
  expected_cols <- c("Mutation", "mut_class", "mut_kabat", "resi", "Pos", "aa", "minus_delta_log_Kd", "ipt_freq", "exp_score")
  existing_cols <- expected_cols[expected_cols %in% colnames(df)]
  df <- dplyr::select(df, dplyr::all_of(existing_cols))
  
  WTresibox <- df %>%
    dplyr::select(resi, Pos) %>%
    unique() %>%
    mutate(
      WT_resi = str_sub(resi, 1, 1),
      x = seq(1, 120),
      y = match(WT_resi, aa_level)
    ) %>%
    dplyr::select(resi, WT_resi, Pos, x, y)
  # === ADD THIS LINE RIGHT HERE ===
  df <- df %>% left_join(WTresibox %>% select(resi, x), by = "resi")
  
  kabat_num <- unique(str_sub(df$mut_kabat, 2, -2))
  kabat_num <- kabat_num[!is.na(kabat_num)]
  kabat_resi <- c()
  for (i in seq(1, 120, 1)) {
    kabat_resi <- c(kabat_resi, paste(WTresibox$WT_resi[i], kabat_num[i], sep = ''))
  }
  
  # SHMresibox <- read_tsv('data/CR9114_SHMheatmap.tsv')
  SHMresibox <- tibble(x = numeric(), y = numeric())
  
  if ("mut_class" %in% colnames(df) & "ipt_freq" %in% colnames(df)) {
    df <- df %>%
      mutate(minus_delta_log_Kd = case_when(
        (mut_class == 'missense' & ipt_freq < 0) ~ 99,
        TRUE ~ minus_delta_log_Kd
      ))
  }
  
  print(range(df$ipt_freq, na.rm = TRUE))
  
  validation_highlight <- if (grepl("GL", infile)) {
    tibble(
      Pos = c(52, 52, 55, 74, 87, 110),
      aa  = c("G", "V", "H", "W", "E", "L")
    )
  } else {
    tibble(
      Pos = c(24, 24, 46, 52, 52, 55, 74, 87, 110),
      aa  = c("F", "L", "E", "G", "V", "H", "W", "E", "L")
    )
  }
  
  nonkabat_highlight <- tibble(
    Pos = c(24, 29, 30, 31, 46, 52, 57, 59, 71, 74, 75, 77, 84, 87, 95, 106),
    aa  = c("S", "S", "N", "N", "D", "S", "S", "A", "S", "I", "F", "N", "N", "T", "F", "S")
  )
  
  # Combine with a new column to label the groups
  combined_highlight_list <- bind_rows(
    nonkabat_highlight %>% mutate(label = "SHM"),
    validation_highlight %>% mutate(label = "mutants tested for validation")
  )
  
  p1 <- plot_freq_heatmap(
    df, WTresibox, SHMresibox, 1, 120, kabat_resi[1:120],
    highlight_list = combined_highlight_list
  )
  
  # Extract legend and plot without it
  legend_grob <- extract_legend(p1)
  p_main <- p1 + theme(legend.position = "none")
  
  # Save main heatmap plot (no legend)
  ggsave(graphname, p_main, width = 7.2, height = 1.7, dpi = 300, bg='white')
  
  # Save legend as separate PNG
  png_filename <- gsub(".png", "_legend.png", graphname)
  png(png_filename, width = 6.5, height = 0.5, units = "in", res = 300, bg='white')
  grid.newpage()
  grid.draw(legend_grob)
  dev.off()
  
}

# Run the wrapper function (uncomment and adjust paths as needed)
wrqpper('result/HC_GL_KD_table_summary_kabat.csv', 'graph/CR9114_GL_ipt_freq_heatmap.png')
wrqpper('result/HC_WT_H1_KD_table_summary_kabat.csv', 'graph/CR9114_WT_H1_ipt_freq_heatmap.png')
wrqpper('result/HC_WT_H3_KD_table_summary_kabat.csv', 'graph/CR9114_WT_H3_ipt_freq_heatmap.png')
wrqpper('result/HC_WT_fluB_KD_table_summary_kabat.csv', 'graph/CR9114_WT_fluB_ipt_freq_heatmap.png')
