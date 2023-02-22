#' A function to make adjcany matrix from edge list and find isolated modules in it
#'
#' @param adjancy_matrix: a matrix P X P
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
get_membership<-function(adjmat){

  require(igraph)
  require(rlist)




g <- graph_from_data_frame(adjmat,directed = F)

compList <- components(g)

Component_number <- max(compList$membership)

membership_Sets <- NULL

for(i in 1:Component_number){

membership_Sets[[i]] <- names(which(compList$membership==i))


}

return(membership_Sets)


}







