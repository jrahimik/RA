#functions for auc and logistic regression 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
create_cell_data <- function(cell, metadata, expression_matrix, cell_specific_genes) {
  #cell_specific_genes <- read.table(paste0(input_dir, '/', cell_file), fill = T, header = T)
  #cell_specific_genes_entrez <- cell_specific_genes[,c(1]
  #subset expression matrix first to the cell type of interest 
  cell_type_metadata <- metadata %>% subset(metadata$type == cell)
  cell_type_expression_matrix <- expression_matrix[,which(colnames(expression_matrix) %in% cell_type_metadata$Tube.label)]
  #cell_type_expression_matrix <- cell_type_expression_matrix[cell_specific_genes[,1],]
  cell_type_metadata <- cell_type_metadata[which(cell_type_metadata$Tube.label %in% colnames(cell_type_expression_matrix)),]
  #get number of 0 = controls and 1 = cases
  n_controls <- length(which(cell_type_metadata$Phenotype == 0))
  n_cases <- length(which(cell_type_metadata$Phenotype == 1))
  sample_size <- n_cases + n_controls
  #match indices to add phenotype labels to expression matrix
  reorder_id <- match(cell_type_metadata$Tube.label, colnames(cell_type_expression_matrix)) ##vector of index of colnames(expression_matrix) that match index in cell_type_metadata$Tube.label
  cell_type_expression_matrix <- cell_type_expression_matrix[, reorder_id]
  ##Transpose the matrix 
  expression_data <- data.frame(t(cell_type_expression_matrix))
  #Get hgnc gene names for the columns
  ##Download dataset for mapping transcript ID to gene ID
  gene_names <- cell_specific_genes[cell_specific_genes$ensembl_gene_id %in% colnames(expression_data),]
  gene_names <- gene_names[gene_names$hgnc_symbol != '',]
  gene_names <- gene_names[!duplicated(gene_names$ensembl_gene_id),]
  expression_data <- expression_data[, colnames(expression_data) %in% gene_names$ensembl_gene_id]
  order_id <- match(gene_names$ensembl_gene_id, colnames(expression_data))
  expression_data <- expression_data[, order_id]
  colnames(expression_data) <- gene_names$hgnc_symbol
  assign("cell_expression_data", expression_data, envir = .GlobalEnv)
  assign("cell_metadata", cell_type_metadata, envir = .GlobalEnv)
}
#create module list specific dataset and size matched negative control#
create_module_specific_dataset <- function(module, expression_data, cell_type_metadata, cell) {
  genes <- module
  gene_overlap <- length(which((colnames(expression_data) %in% genes)))
  #subset expression data to only genes of interest and a random sized matched gene set
  ##all genes might not be present in dataset##
  if (gene_overlap > 1) {
    data <- expression_data[, (colnames(expression_data) %in% genes)]
    print(paste0("no of genes for ", cell, " for module length", " ",  length(genes), " is ", gene_overlap))
    data$Phenotype <- cell_type_metadata$Phenotype
    assign("module_genes", genes, envir = .GlobalEnv)
    assign("module_data", data, envir = .GlobalEnv)
    assign("gene_overlap", gene_overlap, envir = .GlobalEnv)
  } else {
    gene_overlap <- 0
    print(paste0("No gene overlap for", " ", cell))
    assign("gene_overlap", gene_overlap, envir = .GlobalEnv)
  }
}
create_permute_data_set <- function(gene_list, expression_data, cell_type_metadata) {
  #subset expression data to remove all genes across modules to prevent them from being in the permute data 
    data <- expression_data[, !(colnames(expression_data) %in% gene_list)]
    assign("permute_data", data, envir = .GlobalEnv)
}
#run lasso and predict using linear regression using leave one out 
lasso_linear_predict <- function(X, y, k) {
  X_mt <- as.matrix(X)
  require(caret)
  require(glmnet)
  require(pROC)
  reg_coef <- c(0.4, 0.5, 0.6, 0.7, 0.9)
  roc_data <- list()
  auc_df <- data.frame('AUC_val' = numeric(), 'Variables' = character(), 'Reg_coef' = character())
  for (co_ef in 1:length(reg_coef)) {
  Y_true_all <- c()
  Y_pred_all <- c()
  p_val_all <- c()
  selected_vars <- c()
  for (num in 1:nrow(X)) {
    X_train <- X_mt[-num,]
    Y_train <- y[-num]
    Y_test <- y[num]
    variables <- c()
    #p_val <- c()
    #run ten rounds of 5 fold cross-validation to get
    for (num_FS in 1:10) {
      glmnet1 <- cv.glmnet(X_train, Y_train, alpha=1, nfolds = k)
      # tuning lambda
      lambda <- glmnet1$lambda.min 
      lambda <- lambda*reg_coef[co_ef]
      #print(lambda)
      ###################### Using Selected Lambda to Run Lasso ############
      glmnet2 <- glmnet(X_train, Y_train, alpha=1, lambda = lambda)
      c <- coef(glmnet2)
      ##################### Saving Selected Features From Lasso ############
      inds<-which(c!=0)
      #remove intercept from the list
      tmp_variables <- tail(row.names(c)[inds], -1)
      #save variables
      variables<-c(variables,tmp_variables)
    }
    #get number of times each variable was selected and filter
    freq <- sort(table(variables),decreasing=TRUE)/num_FS
    print(freq)
    stable_Var <- names(which(freq>=0.6))  #### the 0.6 can be changed ######
    selected_vars <- c(selected_vars, length(stable_Var))
    print(selected_vars)
    X_train_log <- data.frame(X[-num, names(X) %in% stable_Var])
    names(X_train_log) <-  names(X)[names(X) %in% stable_Var]
    df_Train <- cbind.data.frame(Y_train, X_train_log)
    colnames(df_Train)[1] <- "Phenotype"
    X_test_log <- data.frame(X[num, names(X) %in% stable_Var])
    names(X_test_log) <-  names(X)[names(X) %in% stable_Var]
    linear_model <- lm(Phenotype~., data = df_Train)
    Y_pred <- predict(linear_model, newdata = X_test_log, type = 'response')
    Y_true_all <- c(Y_true_all, Y_test)
    Y_pred_all <- c(Y_pred_all, Y_pred)
  } 
  roc_data[[co_ef]] <- roc(Y_true_all, Y_pred_all)
  names(roc_data)[co_ef] <- reg_coef[co_ef]
  auc_val <- ifelse(roc_data[[co_ef]]$auc < 0.5, 1 - roc_data[[co_ef]]$auc, roc_data[[co_ef]]$auc)
  auc_df <- rbind(auc_df, data.frame('AUC_val' = as.numeric(auc_val), 'Variables' = paste(selected_vars, collapse = ", "), 'Reg_coef' = reg_coef[co_ef]))
  }
  assign('ROC_module_list', roc_data, envir = .GlobalEnv)
  assign('AUC_module_df', auc_df, envir = .GlobalEnv)
}
#get random samples of genes and run lasso followed by linear regression but with the regularization coefficient of max AUC
lasso_linear_permute <- function(X, y, k, module_gene_n, reg_coef) {
  X_mt <- as.matrix(X)
  n_permute <- 1000
  require(caret)
  require(glmnet)
  require(pROC)
  roc_data <- list()
  auc_df <- c()
  for (n in 1:n_permute) {
    print(n)
    X_perm <- X_mt[, sample(colnames(X_mt), size = module_gene_n)]
    Y_true_all <- c()
    Y_pred_all <- c()
    p_val_all <- c()
    selected_vars <- c()
    for (num in 1:nrow(X)) {
      X_train <- X_perm[-num,]
      Y_train <- y[-num]
      Y_test <- y[num]
      variables <- c()
      #p_val <- c()
      #run ten rounds of 5 fold cross-validation to get
      for (num_FS in 1:10) {
        glmnet1 <- cv.glmnet(X_train, Y_train, alpha=1, nfolds = k)
        # tuning lambda
        lambda <- glmnet1$lambda.min 
        lambda <- lambda*reg_coef
        #print(lambda)
        ###################### Using Selected Lambda to Run Lasso ############
        glmnet2 <- glmnet(X_train, Y_train, alpha=1, lambda = lambda)
        c <- coef(glmnet2)
        ##################### Saving Selected Features From Lasso ############
        inds<-which(c!=0)
        #remove intercept from the list
        tmp_variables <- tail(row.names(c)[inds], -1)
        #save variables
        variables<-c(variables,tmp_variables)
      }
      #get number of times each variable was selected and filter
      freq <- sort(table(variables),decreasing=TRUE)/num_FS
      print(freq)
      stable_Var <- names(which(freq>=0.6))  #### the 0.6 can be changed ######
      selected_vars <- c(selected_vars, length(stable_Var))
      print(selected_vars)
      X_train_log <- data.frame(X[-num, names(X) %in% stable_Var])
      names(X_train_log) <-  names(X)[names(X) %in% stable_Var]
      df_Train <- cbind.data.frame(Y_train, X_train_log)
      colnames(df_Train)[1] <- "Phenotype"
      X_test_log <- data.frame(X[num, names(X) %in% stable_Var])
      names(X_test_log) <-  names(X)[names(X) %in% stable_Var]
      linear_model <- lm(Phenotype~., data = df_Train)
      Y_pred <- predict(linear_model, newdata = X_test_log, type = 'response')
      Y_true_all <- c(Y_true_all, Y_test)
      Y_pred_all <- c(Y_pred_all, Y_pred)
    } 
    roc_data[[n]] <- roc(Y_true_all, Y_pred_all)
    auc_val <- ifelse(roc_data[[n]]$auc < 0.5, 1 - roc_data[[n]]$auc, roc_data[[n]]$auc)
    auc_df <- c(auc_df, auc_val)
  }
  return(auc_df)
}
#run linear regression 
linear_run_predict <- function(X, y) {
  Y_true_all <- c()
  Y_pred_all <- c()
  p_val_all <- c()
  for (num in 1:nrow(X)) {
    Y_train <- y[-num]
    Y_test <- y[num]
    variables <- c()
    X_train_log <- data.frame(X[-num, ])
    df_Train <- cbind.data.frame(Y_train, X_train_log)
    colnames(df_Train)[1] <- "Phenotype"
    X_test_log <- data.frame(X[num, ])
    linear_model <- lm(Phenotype~., data = df_Train)
    Y_pred <- predict(linear_model, newdata = X_test_log, type = 'response')
    Y_true_all <- c(Y_true_all, Y_test)
    Y_pred_all <- c(Y_pred_all, Y_pred)
  } 
  Y_labels <- list('Y_pred_module' = Y_pred_all, 'Y_true_module' = Y_true_all)
  return(Y_labels)
}
#run linear regression for permuted data 
linear_run_permute <- function(X, y, module_gene_n) {
  n_permute <- 1000 
  roc_data <- list()
  auc_df <- c()
  for (n in 1:n_permute) {
    X_perm <- X[, sample(colnames(X), size = module_gene_n)]
    Y_true_all <- c()
    Y_pred_all <- c()
    p_val_all <- c()
    for (num in 1:nrow(X)) {
      Y_train <- y[-num]
      Y_test <- y[num]
      variables <- c()
      X_train_log <- data.frame(X_perm[-num, ])
      df_Train <- cbind.data.frame(Y_train, X_train_log)
      colnames(df_Train)[1] <- "Phenotype"
      X_test_log <- data.frame(X_perm[num, ])
      linear_model <- lm(Phenotype~., data = df_Train)
      Y_pred <- predict(linear_model, newdata = X_test_log, type = 'response')
      Y_true_all <- c(Y_true_all, Y_test)
      Y_pred_all <- c(Y_pred_all, Y_pred)
    }
  roc_data[[n]] <- roc(Y_true_all, Y_pred_all)
  auc_val <- ifelse(roc_data[[n]]$auc < 0.5, 1 - roc_data[[n]]$auc, roc_data[[n]]$auc)
  auc_df <- c(auc_df, auc_val)
  }
  return(auc_df)
}
