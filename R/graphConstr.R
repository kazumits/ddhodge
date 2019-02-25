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
graph_altmat <- function(A) {
  A[A<0] <- 0 # drop negatives
  igraph::graph_from_adjacency_matrix(
    Matrix::drop0(A),mode="directed",weighted=TRUE,diag=FALSE
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
#' @param X data matrix (column: samples, rows: variables)
#' @param roots indices of columns of X or logical vector; start points of diffusion (i.e., the initial density is 1/n at specified n points).
#' @param j integer; use j-th nearest neighbouor in local scaling of sigma
#' @param k integer; number of nearest neighbours to construct initial k-NN graph
#' @param s numeric; time parameter of diffusion maps.
#' @param npc integer; number of principal components in calculating euclid distance.
#' @param ndc integer; number of diffusion components in calculating diffusion distance.
#' @param sigma numeric; fixed sigma for isotropic diffusion
#' @param lambda numeric; LASSO penalty of elasticnet
#' @param alpha numeric; mixing perameter of ridge:0 and LASSO:1 (defaut: 0.5)
#' @export
#' @references "Diffusion maps"
#' @references "Self-tuning spectral clustering"
diffusionGraph <- function(X,roots,j=7,k=11,s=1,ndc=10,npc=min(1000,dim(X)),lambda=1e-4,alpha=0.5,sigma=NULL){
  Y <- log(gscale(X+1))
  p <- stats::prcomp(t(Y))
  R <- stats::dist(p$x[,1:npc]) # Distance matrix
  d <- diffusionMaps(R,j,sigma)
  # Diffusion distance matrix
  W <- with(d,as.matrix(dist(Psi%*%diag(eig^s)[,seq(2,ndc+1)])))
  # Transition matrix at t=s
  M <- with(d,Psi%*%diag(eig^s)%*%t(Phi))
  # Potential energy u=-log(p) where p is density at time s
  u <- -log(colMeans(M[roots,,drop=FALSE]))
  names(u) <- colnames(X)
  # approximated gradient along diffusion coordinate
  A <- t(outer(u,u,"-"))/W
  # divergence of original complete graph
  div_o <- drop(scale(div(graph_altmat(A)),center=FALSE))
  A[!(apply(W,1,knn,k) | t(apply(W,1,knn,k)))] <- 0 # k-NN graph
  g <- graph_altmat(A)
  la <- glmnet::glmnet(divop(g),div_o,alpha=alpha,lambda=lambda)
  igraph::E(g)$weight <- drop(la$beta)
  # drop edges with 0 and flip edges with negatives
  g <- graph_altmat(as_altmat(g))
  igraph::V(g)$u <- potential(g)
  igraph::V(g)$div <- div(g)
  return(g)
}
