#' @param X a parameter
#' @return the same thing
#' @export
get_fbesag <- function(graph, id, sd_gamma = 0.2, param = list(p1 = 1, p2 = 1e-5), initial = c(-999)){

  #getit <- paste(find.package("fbesag"), "/wrapper/subregions.R", sep="")
  #source(getit)
  .connection(graph=graph, id = id, sd_gamma = sd_gamma, param = param, initial = initial)
}




