
.wrapper_pbesag <- function(graph, id_s, sd_gamma = 0.2, param = list(p1 = 1, p2 = 1e-5)){

  lam <- - log(param$p2)/param$p1
  n <- dim(graph)[1]
  P <- length(unique(id_s))
  g <- inla.read.graph(graph)

  constr.inter <- list(A = matrix(1,1,dim(graph)[1]), e = rep(0, 1))
  scaled_graph = as.matrix(INLA:::inla.scale.model(graph,constr.inter))
  scaled_cnst = scaled_graph[1,1]/graph[1,1]

  sigm2 <- sd_gamma*sd_gamma
  e = eigen(diag(P) - (1/P)*matrix(1,P,P))
  D = diag(c(1.0/e$values[1:(P-1)]))
  inv_tilda_Sigma = (1/sigm2)*e$vectors[,1:(P-1)]%*%D%*%t(e$vectors[,1:(P-1)])
  val1 <- inv_tilda_Sigma[1,1]
  val2 <- inv_tilda_Sigma[1,2]

  invSig = c(val1,val2,scaled_cnst)

  INLA_CGENERIC_GRAPH <- function(g){
    k = 1
    ii <- c(); jj <- c()

    for(i in 1:g$n){

      ind <- which(g$nbs[[i]]>=i)
      if(length(ind)>0){
        new_n <- g$nbs[[i]][c(ind)]
        ii <- c(ii, i, rep(i,length(new_n)))
        jj <- c(jj, i, new_n)
      }else{
        ii <- c(ii, i)
        jj <- c(jj, i)
      }

    }
    return(c(g$n, length(ii), ii-1,jj-1))

  }
  get_one_vector_final <- function(g){

    prec <- rep(exp(4),2)
    y_scaled = numeric(0)
    k = 1
    ii <- c(); jj <- c()

    one_vector <- c()
    for(i in 1:g$n){

      num_nei_i <- g$nnbs[i]
      one_vector <- c(one_vector, num_nei_i)
      mas_prec <- prec[id_s[i]]
      one_vector <- c(one_vector, id_s[i]-1)
      size_neighbors <- length(g$nbs[[i]])
      one_vector <- c(one_vector, size_neighbors)
      g$nbs[[i]] <- sort(c(g$nbs[[i]]))

      tmp_neighbors <- c()
      sum_tmp_neighbors <- 0
      y_scaled[k] <- 0
      save_k = k
      k = k + 1
      for(j in 1:size_neighbors){
        tick <- g$nbs[[i]][j]
        tmp_neighbors[j] <- -0.5*prec[c(id_s[tick])]
        one_vector <- c(one_vector, id_s[tick]-1)

        one_vector <- c(one_vector, tick-1)
        if(tick>i){
          y_scaled[k] <- tmp_neighbors[j] - 0.5*mas_prec
          k = k + 1
        }
        sum_tmp_neighbors <- sum_tmp_neighbors + tmp_neighbors[j]
      }

      y_scaled[save_k] <- -sum_tmp_neighbors + 0.5*num_nei_i*mas_prec + 1e-5

    }

    return(one_vector)

  }

  VEC_CGENERIC_GRAPH <- INLA_CGENERIC_GRAPH(g)
  VEC_CGENERIC_GRAPH <- c(length(VEC_CGENERIC_GRAPH), INLA_CGENERIC_GRAPH(g))
  final_vec <- get_one_vector_final(g)

  return(list(v1 = VEC_CGENERIC_GRAPH, v2 = final_vec, lam = lam, P = P, n = n, invSig = invSig))
}


.connection <- function(graph, id, sd_gamma, param, initial, useINLAprecomp){

  num_p <- length(unique(id))
  res <- .wrapper_pbesag(graph, id_s = id, sd_gamma = sd_gamma, param = param)
  if(initial[1]==-999){
    initial <- rep(4, num_p)
  }else{
    if(length(initial) != num_p) {
      stop("Error: Initial vector of theta has wrong length")
    }
  }

  libpath <- NULL
  if (is.null(libpath)) {
    if (useINLAprecomp) {
      libpath <- INLA::inla.external.lib("fbesag")
    } else {
      if (Sys.info()["sysname"] == "Windows") {
        files <- list.files(find.package("fbesag"), pattern = "fbesag.dll", recursive = TRUE, full.names = TRUE)
        if (length(files) > 0) {
          file_dir <- dirname(files[1])
          print(file_dir)
        } else {
          print("File not found.")
        }
        libpath <- paste(file_dir, "/fbesag.dll", sep="")
      } else {
        files <- list.files(find.package("fbesag"), pattern = "fbesag.so", recursive = TRUE, full.names = TRUE)
        if (length(files) > 0) {
          file_dir <- dirname(files[1])
          print(file_dir)
        } else {
          print("File not found.")
        }
        libpath <- paste(file_dir, "/fbesag.so", sep="")
      }
    }
  }

  cmodel <- inla.cgeneric.define(model = "inla_cgeneric_pbesag_model",
                                 shlib = libpath, n = as.integer(res$n), npart = res$P, VEC_CGENERIC_GRAPH = as.integer(res$v1), debug = FALSE,  lam=c(res$lam),
                                 invSig = res$invSig, misc = as.integer(res$v2), initial = initial)

  return(cmodel)
}

#' @param X a parameter
#' @return the same thing
#' @export
get_fbesag <- function(graph, id, sd_gamma = 0.2, param = list(p1 = 1, p2 = 1e-5), initial = c(-999), useINLAprecomp = FALSE){

  return(.connection(graph=graph, id = id, sd_gamma = sd_gamma, param = param, initial = initial, useINLAprecomp = useINLAprecomp))
}




