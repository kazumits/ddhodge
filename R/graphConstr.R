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
diffusionGraph <- function(X,roots,k=11,npc=min(1000,dim(X)),ndc=40,s=1,j=7,lambda=1e-4,sigma=NULL){
  Y <- t(log(gscale(X+1)))
  pc <- stats::prcomp(Y,center=TRUE,scale=FALSE)$x
  # Euclid distance matrix
  R <- stats::dist(pc[,1:npc])
  d <- diffusionMaps(R,j,sigma)
  # Diffusion distance matrix
  W <- with(d,as.matrix(dist(Psi%*%diag(eig^s)[,seq(2,ndc+1)])))
  # Transition matrix at t=s
  M <- with(d,Psi%*%diag(eig^s)%*%t(Phi))
  # Set potential as density at time t=s
  u <- colMeans(M[roots,,drop=FALSE])
  # Set potential u=-log(p) where p is density at time t=s
  #u <- -log(colMeans(M[roots,,drop=FALSE]))
  names(u) <- colnames(X)
  # approximated -grad (related to directional deriv.?)
  A <- outer(u,u,"-")
  # divergence of fully-connected diffusion graph
  div_o <- drop(div(graph_altmat(A)))
  # k-NN graph using diffusion distance
  A[!(apply(W,1,knn,k) | t(apply(W,1,knn,k)))] <- 0
  g <- graph_altmat(A)
  # Pulling back the original divergence using pruned graph
  igraph::E(g)$weight <- Matrix::solve(
    Matrix::crossprod(divop(g))+lambda*diag(igraph::ecount(g)),
    -gradop(g)%*%div_o
  )
  igraph::E(g)$weight <- grad(g)
  # drop edges with 0 weights and flip edges with negative weights
  g <- graph_altmat(as_altmat(g))
  igraph::V(g)$u <- potential(g)
  igraph::V(g)$div <- div(g)
  igraph::V(g)$div_o <- div_o
  return(g)
}
