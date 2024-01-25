##Check if the modules obtained have signals using Bulk-RNA and plot au roc curve for only the significant modules for each cell type##
##Started 23.08.2023##
#modified from /ix/djishnu/Priyamvada/Immport/RA_project_git/get_auc_sig_modules_bulk_rna.R#
##Run logistic regression to check if the tpm values of the genes in the modules can classify sample as phenotype 1 vs 0##
##For gene modules that have length(genes) > length(sample), perform LASSO regularized logistic regression##
library(biomaRt)
library(tidyverse)
library(dplyr)
library(pROC)
#source file with function to create cell_specific dataset#
source('/ix/djishnu/Priyamvada/Immport/RA_project_git/functions_auc_logistic.R')
##Input Data##
output_directory <- '/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/LASSO_logistic_new/Updated_LDAK_module/500_permute_results/new_results'
module_dir <- '/ix/djishnu/Priyamvada/Immport/Network_analysis/Network_modules_updated/gene_scores_ldak_500_permute/homosapiens_binary_co_complex_feb2023-ldak_score_ra_prg65'
module_file <- 'significant_modules_members.RDS'
expression_mt <- read.table('/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/low_input_gene_sample_tpm_matrix.725714.tsv')
metadata_df <- read.table('/ix/djishnu/Priyamvada/Immport/RA_inflammatory_paper_analysis/metadata_for_bulk_RNA_seq.tsv')
module_list <- readRDS(paste0(module_dir, '/', module_file))
len <- sapply(module_list, length)
module_list <- module_list[order(len, decreasing = T)]
##For each cell type, subset to the specific genes and compare expression values and use lasso regularized logistic regression
cell_types <- c('B cell', 'Fibro', 'Mono', 'T cell')
#get hgnc gene names 
hsapiens_genes <- getBM(attributes = c("ensembl_gene_id", 
                                           "hgnc_symbol", "gene_biotype"),
                            mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
set.seed(1)
for (cell in 1:length(cell_types)) {
  print(paste0('Processing', ' ', cell_types[cell]))
  create_cell_data(cell = cell_types[cell],
                   metadata_df,
                   expression_mt,
                   hsapiens_genes)
  p_val_pos_df <- data.frame("module" = numeric(), "p_val" = numeric(), "module_len" = numeric(), 'log_p_val' = numeric())
  roc_data <- list()
  for (mod in 1:length(module_list)) {
    module <- mod
    print(paste0('Processing module', ' ', mod))
    create_module_specific_dataset(module = module_list[[mod]],
                                   expression_data = cell_expression_data,
                                   cell_type_metadata = cell_metadata, cell = cell_types[cell])
    if (gene_overlap >= 1) {
      if (nrow(cell_metadata) <= gene_overlap) {
        expression_data_wo_phenotype <- module_data[, !(colnames(module_data) %in% 'Phenotype')] 
        phenotype_labels <- module_data$Phenotype
        k <- length(phenotype_labels)
        module_p_val <- lasso_logistic_pval(X = expression_data_wo_phenotype, 
                                            y = phenotype_labels,
                                            k = k)
        lasso_logistic_predict(X = expression_data_wo_phenotype, 
                       y = phenotype_labels,
                       k = k)
        roc_data[[module]] <- roc(lasso_module_info$Y_true_module, lasso_module_info$Y_pred_module)
        names(roc_data)[module] <- as.character(module)
        #save p-value
        p_val_pos_module <- data.frame("module" = module, "p_val" = module_p_val, 
                                       "gene_overlap" = (ncol(module_data) - 1), 
                                       'log_p_val' = -log10(module_p_val))
      } else {
        expression_data_wo_phenotype <- module_data[, !(colnames(module_data) %in% 'Phenotype')] 
        phenotype_labels <- module_data$Phenotype
        k <- length(phenotype_labels)
        module_p_val <- logit_run_pval(X = expression_data_wo_phenotype, 
                                            y = phenotype_labels)
        logit_run_predict(X = expression_data_wo_phenotype, 
                       y = phenotype_labels)
        roc_data[[module]] <- roc(logit_module_info$Y_true_module, logit_module_info$Y_pred_module)
        names(roc_data)[module] <- as.character(module)
        #save p-value
        p_val_pos_module <- data.frame("module" = module, "p_val" = module_p_val, 
                                       "gene_overlap" = (ncol(module_data) - 1), 'log_p_val' = -log10(module_p_val))
      } 
      print(paste0("running module", " ", module))
      p_val_pos_df <- rbind(p_val_pos_df, p_val_pos_module)
    } else {
      print(paste0("running module", " ", module))
    }
  } 
  #save table with p-value for data set
  #create cell_specific directory if not present
  dir.create(file.path(output_directory, cell_types[cell]), showWarnings = FALSE)
  p_val_pos_df$p_val <- p.adjust(p_val_pos_df$p_val, method = 'BH')
  write.table(p_val_pos_df, file = paste0(output_directory, '/', cell_types[cell], '/', 'module_significance.txt'), col.names = T, row.names = F, sep = '\t')
  #create plot showing distribution of p-value
  sig_modules <- p_val_pos_df[p_val_pos_df$p_val < 0.05, ]
  sig_modules <- sig_modules$module
  roc_data_sig <- roc_data[sig_modules]
  # extract auc
  # Extract AUC values and group names
  auc_values <- roc_data_sig %>%
    map(~ .x$auc)
  group_names <- names(roc_data_sig)
  
  # Create a data frame
  data.auc <- data.frame(
    AUC = unlist(auc_values),
    name = rep(group_names, lengths(auc_values)),
    stringsAsFactors = FALSE
  )
  
  # generate labels labels
  data.auc %>% 
    mutate(label_long=paste0('Module ', name,", AUC = ",paste(round(AUC,2)))) -> data.labels
  plot <- ggroc(roc_data_sig) + 
    scale_color_discrete(labels=data.labels$label_long) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.75), 
          plot.title = element_text(face = 'bold'), axis.text = element_text(face = 'bold', size = 10, color = 'black'), 
          axis.title = element_text(face = 'bold', size = 10, color = 'black'),
          strip.text.x = element_text(size = 10, color = 'black', face = 'bold')) +
    xlab("SPECIFICITY") + ylab("SENSITIVITY") +
    ggtitle(paste0('AUC FOR SIGNIFICANT MODULES FOR', ' ', cell_types[cell]))
  
  # Save the plot to a file
  ggsave(filename = paste0(output_directory, '/', cell_types[cell], '/', 'significant_modules_auc.pdf'), plot = plot, device = 'pdf')
  
  } #closing bracket for module loop
  