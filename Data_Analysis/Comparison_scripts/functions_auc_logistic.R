#functions for auc and logistic regression 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
create_cell_data <- function(cell, metadata, expression_matrix, hsapiens_genes_all) {
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
  gene_names <- hsapiens_genes_all[hsapiens_genes_all$ensembl_gene_id %in% colnames(expression_data),]
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
#run lasso and get p-value for logistic regression
lasso_logistic_pval <- function(X, y, k) {
  X_mt <- as.matrix(X)
  require(caret)
  require(glmnet)
  selected_vars <- c()
  variables <- c()
    #run ten rounds of 5 fold cross-validation to get
    for (num_FS in 1:10) {
      glmnet1 <- cv.glmnet(X_mt, y, alpha=1, nfolds = k)
      # tuning lambda
      lambda <- glmnet1$lambda.min 
      lambda <- lambda*0.9
      #print(lambda)
      ###################### Using Selected Lambda to Run Lasso ############
      glmnet2 <- glmnet(X_mt, y, alpha=1, lambda = lambda)
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
    X_train_log <- data.frame(X[, names(X) %in% stable_Var])
    names(X_train_log) <-  names(X)[names(X) %in% stable_Var]
    df_Train <- cbind.data.frame(y, X_train_log)
    colnames(df_Train)[1] <- "Phenotype"
    logit_model <- glm(as.factor(Phenotype)~., data = df_Train, family = 'binomial')
    residualDeviance <- logit_model$deviance
    nullDeviance     <- logit_model$null.deviance
    devianceDiff     <- nullDeviance - residualDeviance
    df    <- ncol(df_Train) - 1
    model_pval <- pchisq(devianceDiff, df=df, lower.tail=FALSE)
    return(model_pval)
}
#run lasso and predict using logistic regression using leave one out 
lasso_logistic_predict <- function(X, y, k) {
  X_mt <- as.matrix(X)
  require(caret)
  require(glmnet)
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
      lambda <- lambda*0.5
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
    logit_model <- glm(as.factor(Phenotype)~., data = df_Train, family = 'binomial')
    Y_pred <- predict(logit_model, newdata = X_test_log, type = 'response')
    Y_true_all <- c(Y_true_all, Y_test)
    Y_pred_all <- c(Y_pred_all, Y_pred)
  } 
  Y_labels <- list('Y_pred_module' = Y_pred_all, 'Y_true_module' = Y_true_all)
  return(Y_labels)
  #assign("lasso_module_info", Y_labels, envir = .GlobalEnv)
}
#run logistic regression to get model p-val
logit_run_pval <- function(X, y) {
    df_Train <- cbind.data.frame(y, X)
    colnames(df_Train)[1] <- "Phenotype"
    logit_model <- glm(as.factor(Phenotype)~., data = df_Train, family = 'binomial')
    residualDeviance <- logit_model$deviance
    nullDeviance     <- logit_model$null.deviance
    devianceDiff     <- nullDeviance - residualDeviance
    df    <- ncol(df_Train) - 1
    model_pval <- pchisq(devianceDiff, df=df, lower.tail=FALSE)
    return(model_pval)
}
#run leave one out logistic regression to get auc
logit_run_predict <- function(X, y) {
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
    logit_model <- glm(as.factor(Phenotype)~., data = df_Train, family = 'binomial')
    Y_pred <- predict(logit_model, newdata = X_test_log, type = 'response')
    Y_true_all <- c(Y_true_all, Y_test)
    Y_pred_all <- c(Y_pred_all, Y_pred)
} 
  Y_labels <- list('Y_pred_module' = Y_pred_all, 'Y_true_module' = Y_true_all)
  return(Y_labels)
  #assign("logit_module_info", Y_labels, envir = .GlobalEnv)
}