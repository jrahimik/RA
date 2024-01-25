#compare auc for size and distrbution matched modules 
#sig modules - 1, 2, 5, 6, 10 
library(biomaRt)
library(tidyverse)
library(dplyr)
library(doParallel)
library(foreach)
source('/ix/djishnu/Priyamvada/Immport/RA_project_git/functions_auc_logistic_coefficient_regularization.R')
output_directory <- '/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/LASSO_linear_regression/v1'
expression_mt <- read.table('/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/low_input_gene_sample_tpm_matrix.725714.tsv')
metadata_df <- read.table('/ix/djishnu/Priyamvada/Immport/RA_inflammatory_paper_analysis/metadata_for_bulk_RNA_seq.tsv')
cell_specific_gene_dir <- '/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/Cell_specific_gene_list'
#get the distribution of gene scores for the significant modules 
module_dir <- '/ix/djishnu/Priyamvada/Immport/Network_analysis/Network_modules_updated/gene_scores_ldak_500_permute/homosapiens_binary_co_complex_feb2023-ldak_score_ra_prg65'
module_file <- 'significant_modules_members.RDS'
module_set <- readRDS(paste0(module_dir, '/', module_file))
len <- sapply(module_set, length)
module_set <- module_set[order(len, decreasing = T)]
#sig modules are 1, 2, 6, 10 
module_set <- module_set[c(1, 2, 5, 6, 10)]
module_set_genes <- unlist(module_set)
#for all genes get the high, mid and low scores based on quartiles 
ldak_genes_all <- read.table('/ix/djishnu/Javad/ForPriyamvada/LDAK_score.csv')
ldak_genes <- ldak_genes_all[!ldak_genes_all$V1 %in% module_set_genes,]
ldak_genes_top <- ldak_genes_all[ldak_genes_all$V2 >= 2, ]
gene_score_summary <- quantile(ldak_genes_all$V2, probs = c(0.5, 0.75, 1))
ldak_low <- ldak_genes[ldak_genes$V2 <= gene_score_summary[1],]
ldak_med <- ldak_genes[ldak_genes$V2 <= gene_score_summary[2] & ldak_genes$V2 > gene_score_summary[1], ]
ldak_high <- ldak_genes[ldak_genes$V2 <= gene_score_summary[3] & ldak_genes$V2 > gene_score_summary[2], ]
##For each cell type, subset to the specific genes and compare expression values and use lasso regularized logistic regression
cell_types <- c('B cell', 'Fibro', 'Mono', 'T cell')
#set.seed(123)
#get module lists
dist_mod_set <- list()
rand_mod_set <- list()
top_mod_set <- list()
for (m in 1:length(module_set)) {
  module <- module_set[[m]]
  #get distribution of scores 
  ldak_genes_mod <- ldak_genes_all[ldak_genes_all$V1 %in% module,]
  dist <- data.frame('Low' = length(which(ldak_genes_mod$V2 <= gene_score_summary[1])),
                     'Med' = length(which(ldak_genes_mod$V2 <= gene_score_summary[2] & ldak_genes_mod$V2 > gene_score_summary[1])),
                     'High' = length(which(ldak_genes_mod$V2 <= gene_score_summary[3] & ldak_genes_mod$V2 > gene_score_summary[2])))
  dist_mod_set[[m]] <- c(ldak_high[sample(1:nrow(ldak_high), dist$High), "V1"], 
                        ldak_low[sample(1:nrow(ldak_low), dist$Low), "V1"],
                        ldak_med[sample(1:nrow(ldak_med), dist$Med), "V1"])
  rand_mod_set[[m]] <- ldak_genes[sample(1:nrow(ldak_genes), length(module)), 'V1']
  top_mod_set[[m]] <- ldak_genes_top[1:length(module), 'V1']
  
}
module_all <- list('real_module' = module_set, 'dist_module' = dist_mod_set, 'random_module' = rand_mod_set, 'top_mod' = top_mod_set)
module_types <- c('real_module', 'dist_module', 'random_module', 'top_mod')
all_auc_df <- data.frame('AUC_val' = numeric(),'Module_size' = numeric(), "Module_type" = character(), 'Cell_type' = character())
for (cell in cell_types) {
  print(paste0('Processing', ' ', cell))
  cell_gene_list <- read.table(paste0(cell_specific_gene_dir, '/', cell, '_gene_list_top_75.txt'), fill = T, header = T)
  create_cell_data(cell = cell,
                   metadata_df,
                   expression_mt,
                   cell_gene_list)
  roc_data <- list()
  for (m in module_types) {
    module_list <- module_all[[m]]
    for (mod in 1:length(module_list)) {
      module <- mod
      print(paste0('Processing module', ' ', mod))
      create_module_specific_dataset(module = module_list[[mod]],
                                   expression_data = cell_expression_data,
                                   cell_type_metadata = cell_metadata, cell = cell)
      if (nrow(cell_metadata) <= gene_overlap) {
        expression_data_wo_phenotype <- module_data[, !(colnames(module_data) %in% 'Phenotype')] 
        phenotype_labels <- module_data$Phenotype
        k <- length(phenotype_labels)
        lasso_logistic_predict(X = expression_data_wo_phenotype, 
                             y = phenotype_labels,
                             k = k)
        AUC_module_df$Module_size <- length(module_list[[mod]])
        AUC_module_df$Module_type <- m
        AUC_module_df$Cell_type <- cell
        AUC_module_df <- AUC_module_df %>% top_n(1, AUC_val)
        AUC_module_df <- AUC_module_df[, c('AUC_val','Module_size', "Module_type", 'Cell_type')]
        all_auc_df <- rbind(all_auc_df, AUC_module_df)
      } else {
        expression_data_wo_phenotype <- module_data[, !(colnames(module_data) %in% 'Phenotype')] 
        phenotype_labels <- module_data$Phenotype
        k <- length(phenotype_labels)
        true_model <- logit_run_predict(X = expression_data_wo_phenotype, 
                                         y = phenotype_labels)
        roc_data[[module]] <- roc(true_model$Y_true_module, true_model$Y_pred_module)
        names(roc_data)[module] <- as.character(module)
        true_auc <- ifelse(roc_data[[module]]$auc < 0.5, 1 - roc_data[[module]]$auc, roc_data[[module]]$auc)
        AUC_module_df <- data.frame('AUC_val' = true_auc, 'Module_size' = length(module_list[[mod]]), "Module_type" = m,  'Cell_type' = cell)
        all_auc_df <- rbind(all_auc_df, AUC_module_df)
        #save p-value
      } 
    }
  }
  }
write.table(all_auc_df, file = paste0(output_directory, '/', 'three_compare_df.txt'), 
            col.names = T, row.names = F, sep = '\t', quote = F)