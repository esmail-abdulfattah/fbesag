#' @param X a parameter
#' @return the same thing
#' @export
get_fbesag <- function(graph, id, sd_gamma = 0.2, param = list(p1 = 1, p2 = 1e-5), initial = c(-999)){

  #source("/home/abdulfe/R/x86_64-pc-linux-gnu-library/4.2/fbesag/wrapper/subregions.R")
  getit <- paste(find.package("fbesag"), "/wrapper/subregions.R", sep="")
  source(getit)
  print(getit)
  connection(graph=graph, id = id, sd_gamma = sd_gamma, param = param, initial = initial)
}




