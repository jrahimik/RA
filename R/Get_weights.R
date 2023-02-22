## Get_weights

Get_weights <- function(edgeList){

  require(igraph)
  g <- igraph::graph_from_data_frame(edgeList,directed = F)
  degrees  <- degree(g)



  ## Get weight
  membership_Sets <- get_membership(edgeList)
  sprintf("Number of modules are %s",length(membership_Sets))
  weight_Sets <- lapply(membership_Sets,function(x){degrees[intersect(x,names(degrees))]})
  return(weight_Sets)


}

