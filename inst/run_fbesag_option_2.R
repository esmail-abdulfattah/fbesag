
rm(list=ls())

##################################################################################-------> functions
sbesag <- function(R, prec, id_p, scaled_cnst, npart){

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
get_sample_from <- function(nsample, Q, rankdef = 1) {
  e <- eigen(Q, symmetric = TRUE)
  n <- dim(Q)[1]

  x <- matrix(0, nsample, n)
  for(j in 1:nsample) {
    for(i in 1:(n-rankdef)) {
      x[j, ] <- x[j, ] + rnorm(1) * sqrt(1/e$values[i]) * e$vector[, i]
    }
  }

  sam <- c()
  for(i in 1:nsample){
    sam <- c(sam, x[i, ])
  }

  return(sam)
}
##################################################################################-------> functions

#######################-------------------#######################
#-------------------fbesag package -----------------------
#######################-------------------#######################

source("wrapper.R")

###try your compiled file:
source = "fbesag.so" #if linux
#source = "pbesag.dll" #if windows

###if it didn't work try this
#source = INLA::inla.external.lib("fbesag")

#######################-------------------#######################
#-------------------------Models--------------------------------
#######################-------------------#######################

library("INLA")

#partitions ids, number
id_p <- sample(c(1,2,3,4),544, replace =TRUE)
P <- npar <- length(unique(id_p))
nobs <- 544

#standard deviations for the gammas: tau = mean_tau * exp(gamma)
sd_sim <- 0.2
mean_tau <- exp(2)
gamma_vector <- rnorm(npar, mean = 0, sd =sd_sim)
gamma_vector <- gamma_vector - mean(gamma_vector)
tau_vec <- mean_tau*exp(gamma_vector)
log(tau_vec)

#get graph
graph.name = system.file("demodata/germany.graph", package="INLA")
graph = inla.graph2matrix(inla.read.graph(system.file("demodata/germany.graph", package="INLA")))
diag(graph) <- 0
graph <- -graph
diag(graph) <- -rowSums(graph)
graph <- graph[1:nobs, 1:nobs]

#some needed parameters
constr.inter <- list(A = matrix(1,1, dim(graph)[1]), e = rep(0, 1))
scaled_graph = as.matrix(INLA:::inla.scale.model(graph,constr.inter))
scaled_cnst = scaled_graph[1,1]/graph[1,1]

#simulate data?
sim_new_graph <- as.matrix(INLA:::inla.as.sparse(sbesag(inla.read.graph(graph),prec = tau_vec, id_p, scaled_cnst = scaled_cnst)))
nsam <- 1
xsim <- get_sample_from(nsam,sim_new_graph,1)

new_data <- list()
new_data$id <- rep(c(Germany$region), nsam)
new_data$y <- rpois(nsam*nobs, lambda = exp(2 + xsim))
new_data$idx.rep <- rep(1:nsam, each = nobs)

formula1 = y ~ 1 + f(id, model="generic", Cmatrix = inla.as.sparse(scaled_graph),
                     constr = TRUE,
                     rankdef=1,
                     replicate = idx.rep,
                     hyper = list(theta = list("pc.prec", param=c(1,1e-5),
                                               initial = 4)))

ctrc <- list(dic=TRUE, waic=TRUE, cpo=TRUE)
ga.rin1 = inla(formula1,
               data = new_data,
               verbose = TRUE,
               control.compute = ctrc,
               control.inla=list(control.vb = list(enable = FALSE), int.strategy = "eb"),
               family = "poisson")

cmodel_pbesag <- get_fbesag(graph, id = id_p, sd_gamma = sd_sim, param = list(p1 = 1, p2 = 1e-5), source = source)
formula2 = y ~ 1 + f(id, model = cmodel_pbesag, constr= TRUE, rankdef=1, replicate = idx.rep)

ga.rin2 <- inla(formula2,
                data = new_data,
                verbose = TRUE,
                control.compute = ctrc,
                control.inla=list(control.vb = list(enable = FALSE), int.strategy = "eb"),
                family = "poisson")


ga.rin1$internal.summary.hyperpar$mean
ga.rin2$internal.summary.hyperpar$mean




