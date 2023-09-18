library(qqman)
library(ggrepel)
data <- read.table('/ix/djishnu/Javad/ForPriyamvada/gwas_filtered.txt', header = T)
data$log_10_p_val <- -log10(data$P)
plot_data <- data[, c('SNP', 'CHR', 'BP', 'P')]
snps_of_interset <- data[data$log_10_p_val >= 5, 'SNP']
plot_data <- plot_data %>% mutate(is_highlight=ifelse(SNP %in% snps_of_interset, "yes", "no"))
plot_df <- plot_data %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(plot_data, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
plot_df <- plot_df %>% mutate(is_highlight=ifelse(-log10(P) > 8, "yes", "no"))
plot_df <- plot_df %>% mutate(is_annotate=ifelse(-log10(P) > 8, "yes", "no"))
axisdf = plot_df %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
ggplot(plot_df, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#C1C0C8", "black"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 16, 2), limits = c(0, 15)) +  # remove space between plot area and x axis
  #add the geomlines 
  geom_hline(yintercept = 8, linetype = "dashed", color = "#a70000", size = 0.5) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "#ff5252", size = 0.5) +
  #add label 
  geom_point(data=subset(plot_df, is_annotate=="yes"), color="#ff5252", size=1) +
  geom_label_repel( data=subset(plot_df, is_annotate=="yes"), aes(label=SNP), size=6) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.title = element_text(face = 'bold', size = 10, color = 'black')
  ) + 
  xlab("CHROMOSOMES") + ylab("NEGATIVE LOG OF P-VALUE") 
manhattan(plot_data, chr="CHR", bp="BP", snp="SNP", p="P")
#plot only snps with p-val above 
p.adjust(0.05, method = 'BH', n = 14)

