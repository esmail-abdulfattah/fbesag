
# Load initialization script

#dyn.load(system.file("./src", "fbesag.so", package = "fbesag"))

# Define internal functions that use the shared object

'inla.rgeneric.pbesag.model' <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {

  pbesag <- function(R, prec, id_p, scaled_cnst, npart){

    get_ij <- function(g){

      ii <- c(); jj <- c()
      for(i in 1:g$n){
        ii <- c(ii,rep(i,g$nnbs[i]), i)
        jj <- c(jj,g$nbs[[i]], i)
      }
      return(list(i=ii,j=jj))
    }

    g <- inla.read.graph(R)

    x_scaled = c()
    for(i in 1:g$n){
      tmp <- -0.5*prec[c(id_p[g$nbs[[i]]])]
      mas_prec <- prec[id_p[i]]
      x_scaled <- c(x_scaled, tmp - 0.5*mas_prec, -sum(tmp) + 0.5*g$nnbs[i]*mas_prec + 1e-5)
    }
    r <- get_ij(g)
    return(sparseMatrix(i = r$i, j = r$j, x = scaled_cnst*x_scaled))
  }

  npart <- length(unique(id_p))
  dim_theta <- npart

  interpret.theta <- function() {
    res <- list()
    for(i in 1:dim_theta){
      res[[i]] = exp(theta[i])
    }
    return(res)
  }

  graph <- function(){
    require(Matrix)
    return(Matrix(W,sparse=TRUE))
  }

  Q <- function() {
    res <- interpret.theta()
    prec = c()
    for(i in 1:dim_theta){
      prec = c(prec, res[[i]])
    }

    myR <- pbesag(W, prec, id_p, scaled_cnst=scaled_cnst, npart)
    return(inla.as.sparse(myR))
  }

  mu <- function(){return(numeric(0))}

  log.norm.const <- function() {
    return (numeric(0))
  }

  log.prior <- function() {

    p1 = 1; p2 = 1e-5
    lam <- - log(p2)/p1
    res <- interpret.theta()
    prec = c()
    for(i in 1:dim_theta){
      prec = c(prec, res[[i]])
    }

    theta_p = log(prec)
    sigm2 <- sd_sim^2

    mean_theta <- mean(theta_p)
    P = npart
    e = eigen(diag(P) - (1/P)*matrix(1,P,P))
    D = diag(c(1.0/e$values[1:(P-1)]))
    inv_tilda_Sigma = (1/sigm2)*e$vectors[,1:(P-1)]%*%D%*%t(e$vectors[,1:(P-1)])

    res1 <- log(lam) - (lam)*exp(-0.5*mean_theta) -0.5*mean_theta
    res2 <- -0.5*(theta_p-mean_theta)%*%inv_tilda_Sigma%*%(theta_p-mean_theta)
    res <- drop(res1) + drop(res2)
    return(res)


  }

  initial <- function() {
    #return(initial_theta)
    return(rep(4,npart))
  }

  quit <- function() {
    return(invisible())
  }

  res <- do.call(match.arg(cmd), args = list())
  return(res)
}

.wrapper_pbesag <- function(graph, id, sd_gamma = 0.2, param = list(p1 = 1, p2 = 1e-5)){

  lam <- - log(param$p2)/param$p1
  n <- dim(graph)[1]
  P <- length(unique(id))
  g <- inla.read.graph(graph)

  constr.inter <- list(A = matrix(1,1,dim(graph)[1]), e = rep(0, 1))
  scaled_graph = as.matrix(inla.scale.model(graph,constr.inter))
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
      mas_prec <- prec[id_p[i]]
      one_vector <- c(one_vector, id_p[i]-1)
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
        tmp_neighbors[j] <- -0.5*prec[c(id_p[tick])]
        one_vector <- c(one_vector, id_p[tick]-1)

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

#' @param X a parameter
#' @return the same thing
#' @export
get_fbesag <- function(graph, id_p, sd_gamma, param = list(p1 = 1, p2 = 1e-5), initial = c(-999)){

  suppressPackageStartupMessages(library(INLA))

  num_p <- length(unique(id_p))
  res <- .wrapper_pbesag(graph, id = id_p, sd_gamma = sd_gamma, param = list(p1 = 1, p2 = 1e-5))
  if(initial[1]==-999){
    initial <- rep(4, num_p)
  }else{
    if(length(initial) != num_p) {
      stop("Error: Initial vector of theta has wrong length")
    }
  }

  txt_fbesag <- paste(system.file(package="fbesag"),"/fbesag.so",sep="")
  #dyn.load(paste(system.file(package="fbesag"),"/fbesag.so",sep=""))
  setwd(txt_fbesag)

  cmodel <- INLA::inla.cgeneric.define(model = "inla.rgeneric.pbesag.model",
                                 shlib = "fbesag.so", n = as.integer(res$n), npart = res$P, VEC_CGENERIC_GRAPH = as.integer(res$v1), debug = FALSE,  lam=c(res$lam),
                                 invSig = res$invSig, misc = as.integer(res$v2), initial = initial)

  return(cmodel)
}


