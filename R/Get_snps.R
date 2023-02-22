library(dplyr)


map_data <- read.table("C:/Users/javad/Documents/GitHub/RAGit/Data/snp_2_gene.map",header=T)




#' Title
#'
#' @param membership_Sets
#' @param gene_snp_map
#'
#' @return
#' @export
#'
#' @examples
#'
#'
Get_snps <- function(membership_Sets,gene_snp_map){

Sets_snps <- list()

Sets_snps <- sapply(1:length(membership_Sets),function(x){
  map_data[which(map_data$gene.name %in% membership_Sets[[x]]),"snp.name"]},
  simplify = F)


return(Sets_snps)


}
