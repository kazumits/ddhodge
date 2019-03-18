## Hodge decomposition toolkit

#' Gradient operator on graph
#'
#' node -> edge
#' @param g igraph object
#' @export
gradop <- function(g) {
  e <- igraph::get.edges(g,igraph::E(g))
  ne <- igraph::ecount(g)
  Matrix::sparseMatrix(
    i = rep(seq(ne),2),
    j = c(e[,1],e[,2]),
    x = rep(c(-1,1),each=ne), # s_j - s_i
    dims = c(ne,igraph::vcount(g)),
    dimnames = list(NULL,igraph::V(g)$name)
  )
}

#' Divergence operator on graph
#'
#' edge -> node
#' @param g igraph object
#' @export
divop <- function(g) -Matrix::t(gradop(g))

#' Curl operator on graph
#'
#' edge -> triangle (2-simplex)
#' @param g igraph object: directed graph with edge weights
#' @export
curlop <- function(g) {
  triv <- igraph::triangles(g)
  ntri <- length(triv)/3
  if(ntri==0) return(matrix(NA,0,igraph::ecount(g)))
  trie <- apply(matrix(triv,3),2,function(x) igraph::E(g)[igraph::`%--%`(x,x)])
  cc <- apply(trie, 2, function(x) {
    e <- igraph::get.edges(g,x)
    # cyclic direction compared to the first edge
    # either `(tail#1) ->` or `-> (head#1)`
    ifelse(c(
      TRUE,
      e[1,2] == e[2,1] | e[1,1] == e[2,2],
      e[1,2] == e[3,1] | e[1,1] == e[3,2]
    ), +1, -1)
  })
  Matrix::sparseMatrix(
    i = rep(seq(ntri),each=3),
    j = c(trie), # edge ids
    x = c(cc),
    dims = c(ntri,igraph::ecount(g))
  )
}

#' 0-Laplacian of nodes (graph laplacian)
#'
#' node -> node
#' @param g igraph object
#' @export
laplacian0 <- function(g) Matrix::crossprod(gradop(g))

#' 1-Laplacian of edges
#'
#' edge -> edge
#' @param g igraph object
#' @export
laplacian1 <- function(g)
  Matrix::crossprod(curlop(g)) - Matrix::tcrossprod(gradop(g))

#' Calculate potential
#'
#' graph -> numeric vector
#' @param g igraph object
#' @param bias bias term in ridge regression
#' @export
potential <- function(g,bias=0) {
  L <- igraph::laplacian_matrix(igraph::as.undirected(g),weight=NA)
  # Solve the normal equation through QR decomposition
  p <- Matrix::solve(Matrix::qr(L),-div(g))
  #p <- Matrix::solve(Matrix::qr(L + bias*diag(ncol(L))),-div(g))
  #p <- drop(qr.fitted(qr(L + bias*diag(ncol(L))),-div(g)))
  #p <- drop(solve(qr(laplacian0(g)),-div(g)))
  p <- as.numeric(p)
  names(p) <- igraph::V(g)$name
  p - min(p)
}

#' Calculate gradient
#'
#' graph -> numeric vector
#' Extract gradient flow by orthogonal projection
#' @param g igraph object
#' @param bias bias term in ridge regression
#' @export
grad <- function(g,bias=0) as.numeric(gradop(g)%*%potential(g,bias))

#' Calculate divergence
#'
#' graph -> numeric vector
#' @param g igraph object
#' @export
div <- function(g) as.numeric(divop(g)%*%igraph::E(g)$weight)

#' Calculate curl of cliques (triangles)
#'
#' graph -> numeric vector
#' @param g igraph object
#' @export
curl <- function(g) drop(curlop(g)%*%igraph::E(g)$weight)
