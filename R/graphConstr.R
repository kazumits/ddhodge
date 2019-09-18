## Graph constructor

# geomatric (median-of-ratios) scaling
gscale <- function(X) {
  if(any(X<0)) stop("All values should be non-negative.")
  scale(X,
    scale = apply(X/apply(X, 1, function(x) exp(mean(log(x)))), 2,
      function(x) stats::median(x,na.rm=TRUE)
    ), center=FALSE
  )
}

# k nearest neibour
knn <- function(x,k){
  r <- rank(x,ties.method="max")
  r <= k+1 & r > 1
}

# alt. mat. -> graph
graph_altmat <- function(A,tol=1e-7) {
  A[A<0] <- 0 # drop negatives
  igraph::graph_from_adjacency_matrix(
    Matrix::drop0(A,tol),
    mode="directed",
    weighted=TRUE,
    diag=FALSE
  )
}

# graph -> alt. mat.
as_altmat <- function(g,attr="weight") {
  if(igraph::ecount(g)==0)
    return(matrix(0,igraph::vcount(g),igraph::vcount(g)))
  A <- igraph::as_adj(g,attr=attr)
  A - Matrix::t(A)
}

# Construct random alternating matrix: just for test purpose
randomalt <- function(n,p=c(0.9,0.05,0.05)) {
  A <- matrix(sample(c(0,1,2),n*n,prob=p,replace=TRUE),n)
  A-t(A)
}

#' Construction of diffusion graph
#'
#' Construct graph from data matrix based on diffusion maps (Coifman, 2015) with adaptive scaling according to local density.
#' @param X data matrix (columns: n samples, rows: p variables)
#' @param roots indices of columns of X or logical vector; start points of diffusion (i.e., the initial density is 1/n at specified n points).
#' @param k integer; number of nearest neighbours to construct initial k-NN graph
#' @param npc integer; number of principal components in calculating euclid distance.
#' @param ndc integer; number of diffusion components in calculating diffusion distance.
#' @param s numeric; time parameter of diffusion maps.
#' @param j integer; use j-th nearest neighbouor in local scaling of sigma
#' @param lambda numeric; ridge penalty in pulling back divergence
#' @param sigma numeric; fixed sigma for isotropic diffusion
#' @export
#' @references "Diffusion maps"
#' @references "Self-tuning spectral clustering"
diffusionGraph <- function(X,roots,k=11,npc=min(100,dim(X)-1),ndc=40,s=1,j=7,lambda=1e-4,sigma=NULL){
  cat("Normalization: ",file=stderr())
  Y <- t(log(gscale(X+0.5)))
  cat("done.",file=stderr(),fill=TRUE)
  #Y <- Matrix::t(X)
  #pc <- stats::prcomp(Y,center=TRUE,scale=FALSE)$x[,seq(npc)]
  cat("Pre-PCA: ",file=stderr())
  pc <- irlba::prcomp_irlba(Y,npc,center=TRUE,scale=FALSE)$x
  cat("done.",file=stderr(),fill=TRUE)
  # Euclid distance matrix
  R <- stats::dist(pc)
  #M <- affinity(R,normalize = TRUE)
  #A <- M - t(M)
  cat("Diffusionmaps: ",file=stderr())
  d <- diffusionMaps(R,j,sigma)
  cat("done.",file=stderr(),fill=TRUE)
  # Diffusion distance matrix
  W <- with(d,as.matrix(dist(Psi%*%diag(eig^s)[,seq(2,ndc+1)])))
  # Transition matrix at t=s
  M <- with(d,Psi%*%diag(eig^s)%*%t(Phi))
  # Set potential as density at time t=s
  u <- colMeans(M[roots,,drop=FALSE])
  #u <- -log(d$Phi[,1])
  # Set potential u=-log(p) where p is density at time t=s
  #u <- -log(colMeans(M[roots,,drop=FALSE]))
  names(u) <- colnames(X)
  # approximated -grad (related to directional deriv.?)
  A <- outer(u,u,"-")#/W
  #P <- with(d,Psi%*%diag(colSums(outer(1:10,eig,`^`)))%*%t(Phi))
  #P <- with(d,Psi[,1]%*%t(Psi[1,]))
  #P <- with(d,Psi %*% diag(eig-1) %*% diag(sapply(eig,function(x) sum(x^seq(0,100)))) %*% t(Phi)) # totalflow
  #A <- P-t(P)
  # divergence of fully-connected diffusion graph
  g_o <- graph_altmat(A)
  #return(g_o)
  cat("Rewiring: ",file=stderr())
  div_o <- drop(div(g_o))
  # u_o <- drop(potential(g_o)) # time consuming
  # k-NN graph using diffusion distance
  nei <- apply(W,1,knn,k)
  A[!(nei|t(nei))] <- 0
  g <- graph_altmat(A)
  # Pulling back the original potential using pruned graph
  #Lgi <- MASS::ginv(as.matrix(laplacian0(g)))
  #div_s <- glmnet::glmnet(Lgi,u_o,alpha=0.5,lambda=lambda)$beta
  #igraph::E(g)$weight <- gradop(g)%*%Lgi%*%div_s
  # Pulling back the original divergence using pruned graph
  igraph::E(g)$weight <- Matrix::solve(
    Matrix::crossprod(divop(g))+lambda*diag(igraph::ecount(g)),
    -gradop(g)%*%div_o
  )
  igraph::E(g)$weight <- grad(g)
  # drop edges with 0 weights and flip edges with negative weights
  g <- graph_altmat(as_altmat(g))
  cat("done.",file=stderr(),fill=TRUE)
  igraph::V(g)$u <- potential(g)
  #igraph::V(g)$u_o <- u_o
  igraph::V(g)$div <- div(g)
  igraph::V(g)$div_o <- div_o
  return(g)
}

#' Multi-source multi-sink maximum flow
#'
#' Simple extension of `max_flow` for multi-source/sink.
#' @param g igraph object
#' @param source sourse(s)
#' @param sink sink(s)
#' @return igraph object with `E(g)$flow` and `V(g)$pass` attributes.
#' @export
multimaxflow <- function(g,source,sink)
{
  usrc <- igraph::vcount(g) + 1 # super-source
  usink <- usrc + 1 # super-sink
  h <- igraph::add_edges(
    igraph::add_vertices(g,2),
    t(rbind(
      cbind(usrc,source),
      cbind(sink,usink)
    )),
    # not to exceeds any maximum flow
    attr=list(weight=sum(igraph::E(g)$weight))
  )
  mf <- igraph::max_flow(h,usrc,usink,igraph::E(h)$weight)
  igraph::E(h)$flow <- mf$flow
  igraph::V(h)$pass <- igraph::strength(h,mode="in",weights=mf$flow)
  igraph::delete_vertices(h,c(usrc,usink))
}
