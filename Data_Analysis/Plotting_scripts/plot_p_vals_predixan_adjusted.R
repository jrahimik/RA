#Plot p-values#
library(tidyr)
library(ggplot2)
library(stringr)
library(tidyverse)
library(reshape2)
module_p_val <- read.table('/ix/djishnu/Priyamvada/Immport/predixcan/logistic_regression/tissue_results_df_500_permute.tsv', header = T)
module_p_val$Adipose_Subcutaneous <- -log10(p.adjust(module_p_val$Adipose_Subcutaneous, method = 'BH'))
module_p_val$Muscle_Skeletal <- -log10(p.adjust(module_p_val$Muscle_Skeletal, method = 'BH'))
module_p_val$Whole_Blood <- -log10(p.adjust(module_p_val$Whole_Blood, method = 'BH'))
x <- melt(module_p_val, id = 'module')
x$module <- as.character(x$module)
labels <- c(`Adipose_Subcutaneous` = 'SUBCUTANEOUS ADIPOSE TISSUE', `Muscle_Skeletal` = 'SKELETAL MUSCLE TISSUE', `Whole_Blood` = 'WHOLE BLOOD')
#png(paste0('/ix/djishnu/Priyamvada/Immport/predixcan/logistic_regression', '/', 'p_value_hist_all_updated_ldak_500_permute_v2.png'), width = 800, height = 366)
plot <- ggplot(x) +
  geom_bar(aes(x = module, y = value, fill = variable), stat = 'identity', color = 'black', linewidth = 0.75) +
  facet_wrap(.~variable, labeller = as_labeller(labels)) +
  scale_fill_manual(values = c("Adipose_Subcutaneous" = "#063970", "Muscle_Skeletal" ='#c576b5',"Whole_Blood" = '#b5c576')) + 
  scale_x_discrete(limits = c("14", "13", "12", "11", "10", "9", "8", "7", "6", "5", "4", "3", "2", "1"),
                   labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black", size = 0.75) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.75), 
        plot.title = element_text(face = 'bold'), axis.text = element_text(face = 'bold', size = 10, color = 'black'), 
        axis.title = element_text(face = 'bold', size = 10, color = 'black'), legend.position = 'none', strip.background = element_blank(),
        strip.text.x = element_text(size = 10, color = 'black', face = 'bold')) +
  xlab("MODULES") + ylab("NEGATIVE LOG OF P-VALUE") +
  ggtitle(paste0('SIGNIFICANT MODULES FOR', ' ', 'DIFFERENT TISSUE TYPES'))
ggsave(file = 'p_vals_predixan_adjusted.png',
       plot = plot,
       device = 'png',
       path = '/ix/djishnu/Priyamvada/Immport/final_plots',
       width = 9,
       height = 3,
       units = 'in'
)

  

