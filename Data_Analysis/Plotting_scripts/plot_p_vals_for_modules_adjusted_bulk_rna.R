#Plot p-values#
library(tidyr)
library(ggplot2)
library(stringr)
cell_dirs <- list.dirs('/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/LASSO_logistic_new/Updated_LDAK_module/500_permute_results', recursive = F)
cell_names <- str_remove(cell_dirs, '/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/LASSO_logistic_new/Updated_LDAK_module/500_permute_results/')
colors <- c("B cell" = '#EF8318', "Fibro" ='#E60D7A',"Mono" = '#009900',"T cell" = '#9999FF')
#modules_to_keep <- c("2", "3", "5", "7", "8", "12", "13")
#reordering the modules from largest to smallest
p_val_df_plot <- data.frame("module" = character(), "p_val" = numeric(),"gene_overlap" = numeric(),"cell_type" = character())
for (cell in 1:length(cell_names)) {
  #get p-value table
  p_val_df <- read.table(paste0(cell_dirs[cell], '/', 'module_significance.txt'), header = T)
  #p_val_df$p_val <- p.adjust(p_val_df$p_val, method = "BH", n = length(p_val_df$p_val))
  p_val_df$module <- as.character(p_val_df$module)
  p_val_df$log_10_p_vals <- 10^-(p_val_df$p_val)
  p_val_df$p_val <- p.adjust(p_val_df$log_10_p_vals, method = 'BH') #adjust p-values 
  p_val_df$p_val <- -log10(p_val_df$p_val)
  p_val_df$cell_type <- cell_names[cell]
  p_val_df <- p_val_df[, c('module', 'p_val', 'gene_overlap', 'cell_type')]
  p_val_df_plot <- rbind(p_val_df_plot, p_val_df)
}
  #plot histogram
cell_type_labels <- c("MONOCYTES", "FIBROBLASTS", "B CELLS", "T CELLS")
names(cell_type_labels) <- c("Mono", "Fibro", "B cell", "T cell")

#png(paste0('/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/LASSO_logistic_new/Updated_LDAK_module/500_permute_results', '/', 'p_value_hist_adjusted.png'), width = 513, height = 366)
plot <- ggplot(p_val_df_plot) +
          geom_bar(aes(x = module, y = p_val, fill = cell_type), stat = 'identity', color = 'black', size = 0.75) +
          facet_wrap(.~cell_type, labeller = labeller(cell_type = cell_type_labels)) +
          scale_fill_manual(values = c("B cell" = '#EF8318', "Fibro" ='#E60D7A',"Mono" = '#009900',"T cell" = '#9999FF')) + 
          scale_x_discrete(limits = c("14", "13", "12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"),
                           labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")) + #changed the order of modules
          geom_hline(yintercept = 1.3, linetype = "dashed", color = "black", size = 0.75) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.75), 
          plot.title = element_text(face = 'bold'), axis.text = element_text(face = 'bold', size = 10, color = 'black'), 
          axis.title = element_text(face = 'bold', size = 10, color = 'black'), legend.position = 'none', strip.background = element_blank(),
          strip.text.x = element_text(size = 10, color = 'black', face = 'bold')) +
          xlab("MODULES") + ylab("NEGATIVE LOG OF P-VALUE") +
          ggtitle(paste0('SIGNIFICANT MODULES FOR', ' ', 'DIFFERENT CELL TYPES'))
dev.off()
ggsave(filename = 'p_value_hist_adjusted_bulk_rna.png',
       plot = plot,
       device = 'png',
       path = '/ix/djishnu/Priyamvada/Immport/final_plots',
       width = 5.75,
       height = 4,
       units = 'in'
       
  )
