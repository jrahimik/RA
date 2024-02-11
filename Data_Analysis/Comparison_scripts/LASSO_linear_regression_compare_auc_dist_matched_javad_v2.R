#library(biomaRt)
library(tidyverse)
library(dplyr)
library(doParallel)
library(foreach)
library(caret)


################################################################################
#                             Reading Data
################################################################################

setwd("~/RA/")
source("Data_Analysis/Comparison_scripts/Funcs.R")



output_directory <- './Data_Analysis/results'


expression_mt <- read.table('./Data/low_input_gene_sample_tpm_matrix.725714.tsv')

metadata_df <- read.table('./Data/metadata_for_bulk_RNA_seq.tsv')

cell_specific_gene_dir <- './Data/Cell_specific_gene_list/'

module_dir <- './Data'

module_file <- 'significant_modules_members.RDS'

module_set <- readRDS(paste0(module_dir, '/', module_file))

len <- sapply(module_set, length)

module_set <- module_set[order(len, decreasing = T)]

output_directory <- './RA/Data_Analysis/results'

################################################################################
#                              Forming LDAK genes
################################################################################
module_set_genes <- unlist(module_set)

genes_all        <- read.table('./Data/LDAK_score.csv')

noModule         <- genes_all[!genes_all$V1 %in% module_set_genes,]

ldak_genes_top   <- noModule[noModule$V2 >= 2, ]

moduleandtopLDAK <- c(module_set_genes,ldak_genes_top$V1)

randomgeneall    <- genes_all[!genes_all$V1 %in% moduleandtopLDAK,]

################################################################################
#                              Running Lasso
################################################################################


cell_types          <- c('B cell', 'Fibro', 'Mono', 'T cell')
#module_types        <- c("real_module","dist_module","random_module","top_mod" )
module_types <- c("real_module","dist_module","random_module")

#names(module_types) <- c("real_module","dist_module","random_module","top_mod" )
names(module_types) <- c("real_module","dist_module","random_module")
#module_types        <- c("top_mod" )


for( k in 1:length(module_set)){

report <- NULL

for(cell in cell_types){

for(rep in 1:100){


     cell_gene_list <- read.table(paste0(cell_specific_gene_dir, '/', cell, '_gene_list_top_75.txt'), fill = T, header = T)

     cell_data  <- create_cell_data(cell = cell,
                                    metadata=metadata_df,
                                    expression_matrix=expression_mt,
                                    cell_specific_genes=cell_gene_list)

     for(typ in module_types){


     if(typ=="real_module"){



       cat(paste0("processing module ",k,"\n\n"))
       cat(paste0("processing cell type",cell,"\n"))



    cellmoData      <- create_module_specific_dataset(module = module_set[[k]],
                                                   expression_data = cell_data$cell_expression_data,
                                                cell_type_metadata = cell_data$cell_metadata,
                                                              cell = cell) ## This is the data for a specific cell and specific module
    if(class(cellmoData)=="numeric"){

      print("Yes")

      next

      }

     COLS <- colnames(cellmoData$data)
     X    <- cellmoData$data[,COLS!="Phenotype"]

     modulesgene <- colnames(X)

     y    <- cellmoData$data[,COLS=="Phenotype",drop=F]

     data <- data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)


     myControl <- caret::trainControl(
       method = "LOOCV",
       summaryFunction = twoClassSummary,
       classProbs = TRUE # IMPORTANT!
     )

     # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

     # Train glm with custom trainControl: model
     model <- train(y ~., data, method = "glmnet", trControl = myControl)

     print(max(model$results$ROC))

     cat(paste0(typ,"\n"))

     report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,1-max(model$results$ROC),max(model$results$ROC))))


     }


       if(typ=="random_module"){


         if(class(cellmoData)=="numeric"){

           print("Yes")

           next

         }

         ss <- sample(intersect(colnames(cell_data$cell_expression_data),noModule$V1),size=cellmoData$overalp)


         RandomcellmoData     <-      create_module_specific_dataset(module = ss,
                                                           expression_data = cell_data$cell_expression_data,
                                                           cell_type_metadata = cell_data$cell_metadata,
                                                           cell = cell)



         COLS <- colnames(RandomcellmoData$data)

         X    <- RandomcellmoData$data[,COLS!="Phenotype"]

         y    <- RandomcellmoData$data[,COLS=="Phenotype",drop=F]


         cat(paste0(typ,"\n\n"))

         data <-data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)
            myControl <- trainControl(
              method = "LOOCV",
           number =10,
           summaryFunction = twoClassSummary,
           classProbs = TRUE # IMPORTANT!
           # verboseIter = TRUE
         )

         # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

         model <- train(y ~., data, method = "glmnet", trControl = myControl)


         report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,1-max(model$results$ROC),max(model$results$ROC))))



       }


       if(typ=="dist_module"){


         if(class(cellmoData)=="numeric"){

           print("Yes")

           next

         }

         ss <- summary(noModule$V2)

         high_threshold    <- ss[5]
         medium_threshold  <- ss[3]
         # Categorize the vector elements into High and Low
         high <- noModule[noModule$V2 >= high_threshold,]
         medium <- noModule[noModule$V2 < high_threshold & noModule$V2 >= medium_threshold,]
         low <- noModule[ noModule$V2 < medium_threshold,]

         actual_mod_score <- genes_all[genes_all$V1 %in% modulesgene,]

           highgenesLen <-   sum(actual_mod_score$V2 > high_threshold )
           mediumgenesLen <- sum(actual_mod_score$V2 < high_threshold & actual_mod_score$V2 >= medium_threshold )
           lowLen <-  sum( actual_mod_score$V2 < medium_threshold )


           ss <- NULL
           if(highgenesLen!=0){ss1 <- sample(intersect(colnames(cell_data$cell_expression_data),high$V1),size=highgenesLen)
           ss <- c(ss,ss1)
           }
           if(mediumgenesLen!=0){ss2 <- sample(intersect(colnames(cell_data$cell_expression_data),medium$V1),size=mediumgenesLen)
           ss <- c(ss,ss2)
           }
           if(lowLen!=0){ss3 <- sample(intersect(colnames(cell_data$cell_expression_data),medium$V1),size=lowLen)
           ss <- c(ss,ss3)
           }





         #ss <- sample(intersect(colnames(cell_data$cell_expression_data),noModule$V1),size=cellmoData$overalp)




         DistcellmoData     <-      create_module_specific_dataset(module = ss,
                                                                     expression_data = cell_data$cell_expression_data,
                                                                     cell_type_metadata = cell_data$cell_metadata,
                                                                     cell = cell)



         COLS <- colnames(DistcellmoData$data)

         X    <- DistcellmoData$data[,COLS!="Phenotype"]

         y    <- DistcellmoData$data[,COLS=="Phenotype",drop=F]


         cat(paste0(typ,"\n\n"))

         data <-data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)
         myControl <- trainControl(
           method = "LOOCV",
           number =10,
           summaryFunction = twoClassSummary,
           classProbs = TRUE # IMPORTANT!
           # verboseIter = TRUE
         )

         # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

         model <- train(y ~., data, method = "glmnet", trControl = myControl)


         report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,1-max(model$results$ROC),max(model$results$ROC))))



       }



if(typ=="top_mod"){


         if(class(cellmoData)=="numeric"){

           print("Yes")

           next

         }


        topNoMod <- intersect(noModule$V1,ldak_genes_top$V1)
        ss <- sample(topNoMod,size=cellmoData$overalp)


         topcellmoData     <-      create_module_specific_dataset(module = ss,
                                                                     expression_data = cell_data$cell_expression_data,
                                                                     cell_type_metadata = cell_data$cell_metadata,
                                                                     cell = cell)



         COLS <- colnames(topcellmoData$data)

         X    <- topcellmoData$data[,COLS!="Phenotype"]

         y    <- topcellmoData$data[,COLS=="Phenotype",drop=F]


         cat(paste0(typ,"\n\n"))

         data <-data.frame(y=factor(y$Phenotype,levels=c(0,1),labels = c("X0","X1")),X)
         myControl <- trainControl(
           method = "LOOCV",
           number =10,
           summaryFunction = twoClassSummary,
           classProbs = TRUE # IMPORTANT!
           # verboseIter = TRUE
         )

         # After creating a custom trainControl object, it's easy to fit caret models that use AUC rather than accuracy to tune and evaluate the model. You can just pass your custom trainControl object to the train() function via the trControl argument.

         model <- train(y ~., data, method = "glmnet", trControl = myControl)


         report <-  rbind(report,data.frame(typ=typ,k=k,cell=cell,auc=ifelse(max(model$results$ROC)<0.5,1-max(model$results$ROC),max(model$results$ROC))))



       }











       }

     }

}

  # Define file name based on k and cell
  file_name <- paste("plot_k", k, ".pdf", sep = "")

  # Add p-values to your plot
  p <- ggplot(report, aes(x = typ, y = auc, fill = typ)) +
    geom_boxplot() +
    facet_wrap(~ cell, scales = "free_y", ncol = 2) +
    labs(x = "Module Type", y = "AUC", title = paste0("Comparison of different settings by Cell Type for Moduel",k),
         fill = "Type") +
    theme_minimal()

  ggsave(file = file_name, plot = p, width = 6, height = 4, dpi = 300)
  saveRDS(report,file=paste0(file_name,".rds"))


  }



## Show this for the power of the test
n <- dim(X)[1]
predictions <- c()

for (i in 1:n) {
  model <- glmnet(X[-i,], y$Phenotype[-i], family="binomial",alpha = 1,lambda = 0.01430420)
  predictions <- c(predictions, predict(model, newx=X[i,]))
}

library(pROC)

roc(y$Phenotype, predictions)






























