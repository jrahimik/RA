#' Title
#'
#' @param g
#'
#' @return degree of nodes
#' @export
#'
#' @examples
Get_degree <- function(g){


compList <- components(g)
degree_vals <- degree(g)

return(degree_vals)

}
