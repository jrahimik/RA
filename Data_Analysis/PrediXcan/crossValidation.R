library(caret)
library(dplyr)
library(rlist)
library(randomForest)
library(pROC)

# Example data frames
# The results of this scripts is saved at : "/ix/djishnu/Javad/RA/Data_Analysis/PrediXcan/results/aucresults.RDS"


data <- read.table("/ix/djishnu/Akash/RA/RA_CyTOF/en_Adipose_Subcutaneous_predicted_expression_gene_mapped.tsv",
                   sep = "\t",
                   header=T)
famFile = read.delim("/ix/djishnu/Javad/racer/racer-genetics/racer-imputed/racer-mega-imputed-f-matched-v3-QC.fam",sep = " ",header = FALSE)


y = famFile[,6]-1 # -1 to convert from {1,2} to {0,1} for glm

k=11

module_res <- readRDS("/ix/djishnu/Javad/racer/racer-genetics/racer-imputed/LDAK_reml_analysis/modules_08_07_2023/500_permutations/significant_modules_members.RDS")

#moduleRes <- list()
modules_out <- vector("list", length = 14)



for(module_number in 1:14){

  cat(paste0("Processing module ",module_number," ... \n"))

  module_x <- tryCatch({
    # Your code that might produce an error
    module_x <- data %>% select(colnames(data)[colnames(data) %in% module_res[[module_number]]])

    # Rest of your code within the loop
    # ...
  }, error = function(e) {
    # Handle the error if needed (or do nothing)
    print(paste("An error occurred in iteration", module_number))

    # Set the flag to indicate an error
    error_flag <- TRUE
    return(error_flag)
  })


  if(ncol(module_x)==0){

next
  }
# Continue with the rest of your code

module_data <- data.frame(y=y,module_x)


x <- module_x

# Specify the training control parameters for cross-validation
ctrl <- trainControl(
  method = "cv",         # Cross-validation
  number = 5,            # Number of folds
  verboseIter = TRUE
)

# Train a random forest model using train function
model <- train(
  x = x,
  y = as.factor(y),
  method = "glm",         # logistic regression
  trControl = ctrl,
  metric='Accuracy'
)



modules_out[[module_number]] <- list(result=model$results,finalModel=model$finalModel)


}


#### Manual implementation ######################################################

noEmptymod <- which(lapply(modules_out,length)!=0)
aucList <- vector("list",length=length(modules_out))


for(module_number in noEmptymod){

cat(paste0("Processing module ",module_number," ... \n"))

module_x <- tryCatch({
  # Your code that might produce an error
  module_x <- data %>% select(colnames(data)[colnames(data) %in% module_res[[module_number]]])

  # Rest of your code within the loop
  # ...
}, error = function(e) {
  # Handle the error if needed (or do nothing)
  print(paste("An error occurred in iteration", module_number))

  # Set the flag to indicate an error
  error_flag <- TRUE
  return(error_flag)
})


if(ncol(module_x)==0){

  next
}






## Cross validation

cvList <- createFolds(y,k = k)

valList <- list()

for(i in 1:k){

train_idx <- unlist(cvList[-i])
valid_idx <- cvList[[i]]


#Model     <- glm(y~.,data=module_data[train_idx,],family = "binomial")

model     <- modules_out[[module_number]]$finalModel
predictY  <- predict(model,module_x[valid_idx,,drop=F])

valList <- c(valList,list(cbind(predicted=predictY,actual=y[valid_idx])))

}

valDf  <- do.call(rbind,valList)
aucVal <- auc(valDf[,"actual"],valDf[,"predicted"])
aucList[[module_number]] <- aucVal



}



################################################################################
# Plotting the results
################################################################################
unlist(aucList)
