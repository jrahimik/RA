#Get p-value and AUC from linear regression 
##For gene modules that have length(genes) > length(sample), perform LASSO regularized logistic regression##
library(biomaRt)
library(tidyverse)
library(dplyr)
library(doParallel)
library(foreach)
#source file with function to create cell_specific dataset#
source('/ix/djishnu/Priyamvada/Immport/RA_project_git/functions_auc_logistic_coefficient_regularization.R')
##Input Data##
output_directory <- '/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/LASSO_logistic_permute/v2'
module_dir <- '/ix/djishnu/Priyamvada/Immport/Network_analysis/Network_modules_updated/gene_scores_ldak_500_permute/homosapiens_binary_co_complex_feb2023-ldak_score_ra_prg65'
module_file <- 'significant_modules_members.RDS'
expression_mt <- read.table('/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/low_input_gene_sample_tpm_matrix.725714.tsv')
metadata_df <- read.table('/ix/djishnu/Priyamvada/Immport/RA_inflammatory_paper_analysis/metadata_for_bulk_RNA_seq.tsv')
cell_specific_gene_dir <- '/ix/djishnu/Priyamvada/Immport/Bulk_rna_nov_29/Cell_specific_gene_list'
module_list <- readRDS(paste0(module_dir, '/', module_file))
len <- sapply(module_list, length)
module_list <- module_list[order(len, decreasing = T)]
##For each cell type, subset to the specific genes and compare expression values and use lasso regularized logistic regression
cell_types <- c('B cell', 'Fibro', 'Mono', 'T cell')
#set.seed(123)
num_cores <- detectCores() # Set the desired number of cores

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)
all_module_genes <- unlist(module_list)
get_auc_p_vals <- function(cell, metadata_df, expression_mt, cell_specific_gene_dir, all_module_genes, module_list) {
  all_auc_df <- data.frame('AUC_val' = numeric(), 'Variables' = character(), 'Reg_coef' = character(), 'Module' = numeric(), 'Cell_type' = character(), 'Sample' = numeric())
  p_val_df <- data.frame('AUC_val' = numeric(), 'p_val' = numeric(), 'Module' = numeric(), 'Cell_type' = character())
  print(paste0('Processing', ' ', cell))
  cell_gene_list <- read.table(paste0(cell_specific_gene_dir, '/', cell, '_gene_list_top_75.txt'), fill = T, header = T)
  create_cell_data(cell = cell,
                   metadata_df,
                   expression_mt,
                   cell_gene_list)
  create_permute_data_set(gene_list = all_module_genes,
                          expression_data = cell_expression_data,
                          cell_type_metadata = cell_metadata)
  roc_data <- list()
  for (mod in 1:length(module_list)) {
    module <- mod
    #print(paste0('Processing module', ' ', mod))
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
      AUC_module_df$Module <- mod
      AUC_module_df$Cell_type <- cell
      AUC_module_df$Sample <- k
      #roc_data[[module]] <- roc(true_model$Y_true_module, true_model$Y_pred_module)
      #names(roc_data)[module] <- as.character(module)
      #true_auc <- ifelse(roc_data[[module]]$auc < 0.5, 1 - roc_data[[module]]$auc, roc_data[[module]]$auc)
      #save p-value
      reg_coef_permute <- AUC_module_df[AUC_module_df$AUC_val == max(AUC_module_df$AUC_val), 'Reg_coef']
      print(reg_coef_permute)
      auc_mod <-  max(AUC_module_df$AUC_val)
      permute_auc <- lasso_logistic_permute(X = permute_data, 
                                          y = phenotype_labels, k = k,
                                          module_gene_n = gene_overlap, 
                                          reg_coef = reg_coef_permute) 
      permute_p_val <- (length(permute_auc[permute_auc >= auc_mod])/length(permute_auc))
      p_val_df <- rbind(p_val_df, data.frame('AUC_val' = auc_mod, 'p_val' = permute_p_val, 'Module' = mod, 'Cell_type' = cell))
      all_auc_df <- rbind(all_auc_df, AUC_module_df)
    } else {
      expression_data_wo_phenotype <- module_data[, !(colnames(module_data) %in% 'Phenotype')] 
      phenotype_labels <- module_data$Phenotype
      k <- length(phenotype_labels)
      true_model <- logit_run_predict(X = expression_data_wo_phenotype, 
                                       y = phenotype_labels)
      roc_data[[module]] <- roc(true_model$Y_true_module, true_model$Y_pred_module)
      names(roc_data)[module] <- as.character(module)
      auc_mod <- ifelse(roc_data[[module]]$auc < 0.5, 1 - roc_data[[module]]$auc, roc_data[[module]]$auc)
      permute_auc <- logit_run_permute(X = permute_data, 
                                        y = phenotype_labels,
                                        module_gene_n = gene_overlap) 
      permute_p_val <- (length(permute_auc[permute_auc >= auc_mod])/length(permute_auc))
      p_val_df <- rbind(p_val_df, data.frame('AUC_val' = auc_mod, 'p_val' = permute_p_val, 'Module' = mod, 'Cell_type' = cell))
      AUC_module_df <- data.frame('AUC_val' = auc_mod, 'Variables' = gene_overlap, 'Reg_coef' = 0, 'Module' = mod, 'Cell_type' = cell, 'Sample' = k)
      all_auc_df <- rbind(all_auc_df, AUC_module_df)
      #save p-value
    } 
  }
  cell_auc_data <- list('p_val_df' = p_val_df,  'all_auc_df' = all_auc_df)
  return(cell_auc_data)
}
auc_cell_data <- foreach(cell = cell_types, .packages = c("pROC", "foreach", "caret", "tidyverse"), .combine = 'list', .multicombine=TRUE) %dopar% {
  get_auc_p_vals(cell, metadata_df, expression_mt, cell_specific_gene_dir, all_module_genes, module_list)
}
combined_p_val_df <- do.call(rbind, lapply(auc_cell_data, `[[`, 'p_val_df'))
combined_auc_df <- do.call(rbind, lapply(auc_cell_data, `[[`, 'all_auc_df'))
write.table(combined_p_val_df, file = paste0(output_directory, '/', 'permute_p_val.txt'), sep = '\t', col.names = T, row.names = F, quote = F)
write.table(combined_auc_df, file = paste0(output_directory, '/', 'combined_auc.txt'), sep = '\t', col.names = T, row.names = F, quote = F)