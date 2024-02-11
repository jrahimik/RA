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
  # assign("cell_expression_data", expression_data, envir = .GlobalEnv)
  # assign("cell_metadata", cell_type_metadata, envir = .GlobalEnv)
  
  list(cell_expression_data=expression_data,cell_metadata=cell_type_metadata)
  
}

create_module_specific_dataset <- function(module, expression_data, cell_type_metadata, cell) {
  genes <- module
  gene_overlap <- length(which((colnames(expression_data) %in% genes)))
  #subset expression data to only genes of interest and a random sized matched gene set
  ##all genes might not be present in dataset##
  if (gene_overlap > 1) {
    data <- expression_data[, (colnames(expression_data) %in% genes)]
    print(paste0("no of genes for ", cell, " for module length", " ",  length(genes), " is ", gene_overlap))
    data$Phenotype <- cell_type_metadata$Phenotype
    return(list(genes=genes,data=data,overalp=gene_overlap)) 
  } else {
    gene_overlap <- 0
    print(paste0("No gene overlap for", " ", cell))
    return(0)
  }
}


  # X: input data
  # Y: labels
  # K: Number of k folds
  
library(glmnet)
library(pROC)

# Function to perform k-fold cross-validation for logistic regression with Lasso penalty
k_fold_cv <- function(X, y, nfolds = 5, alpha = 1) {
  # Initialize vectors to store predictions and actual labels
  y_pred <- matrix(-1000, nrow=dim(y)[1],ncol=1)
  y_true <- matrix(-1000, nrow=dim(y)[1],ncol=1)
  n <- dim(y)[1]
  
  
  
  # Create fold indices
  #foldid <- sample(rep(seq_len(nfolds), length = dim(y)[1]))
  # Fit the model on the training set
  cv_fit <- cv.glmnet(x = as.matrix(X), y = as.matrix(y), 
                      family = "binomial", alpha = alpha,type.measure = "auc")
  
  cat(paste0("Print the dim of X is ",dim(as.matrix(X))))
  best_lambda <- cv_fit$lambda.min
  
  print(sprintf("The best lambda is %s",best_lambda))
  # Perform k-fold cross-validation
  for (i in 1:n) {
    # Select indices for the current fold
    fold_indices <- i
    
    # Split the data into training and validation sets
    X_train <- X[-fold_indices, , drop = FALSE]
    y_train <- y[-fold_indices,,drop=F]
    X_valid <- X[fold_indices, , drop = FALSE]
    y_valid <- y[fold_indices,,drop=F]
    
   
    
    # Get the optimal lambda value selected by cross-validation
    
    
    # Fit the final model with the optimal lambda
    final_model <- glmnet(x = as.matrix(X_train), y = as.matrix(y_train), 
                          family = "binomial", alpha = alpha, lambda = best_lambda)
    
    # Make predictions on the validation set
    y_pred[fold_indices,1] <- predict(final_model, newx = as.matrix(X_valid), 
                                    s = best_lambda, type = "response")
    y_true[fold_indices,1] <- as.matrix(y_valid)
  }
  
  # Calculate AUC
  auc <- roc(y_true, y_pred)
  
  # Return AUC
  return(auc$auc)
}



  
  
  
  
  
  
  
  
  
  






