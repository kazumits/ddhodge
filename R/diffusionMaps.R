# Diffusion maps

logsumexp <- function(x,a=max(x)) a + log(sum(exp(x - a)))

#' Gaussian affinity matrix constructor
#'
#' equivalent to diffusion maps (Coifman, 2015) with adaptive scaling according to local density.
#' @param R symmetric matrix; distance or any similarity matrix.
#' @param k integer; number of neighbours in adaptive-scaling. (ignored if sigma is specified)
#' @param sigma numeric; fixed sigma for isotropic diffusion
#' @param normalize logical; scale s.t. row summing to 1 (i.e. returns transition matrix)
#' @param log logical; if TRUE, log-transformed affinity matrix is returned.
#' @export
#' @references "Self-tuning spectral clustering"
affinity <- function(R,k=7,sigma=NULL,log=FALSE,normalize=FALSE) {
  R <- as.matrix(R)
  if(is.null(sigma)){
    s <- apply(R,2,function(x) sort(x)[k])
    S <- sqrt(outer(s,s))
  } else S <- sigma # isotropic
  logW <- -(R/S)^2
  if(normalize) logW <- sweep(logW,1,apply(logW,1,logsumexp))
  if(log) logW else exp(logW)
}

#' Diffusion maps
#'
#' equivalent to diffusion maps (Coifman, 2015) with adaptive scaling according to local density.
#' @param R distance or any similarity matrix.
#' @param k number of neighbours in adaptive-scaling. (ignored if sigma is specified)
#' @param sigma fixed sigma for isotropic diffusion
#' @export
#' @references "Self-tuning spectral clustering"
diffusionMaps <- function(R,k=7,sigma=NULL){
  logW <- affinity(R,k,sigma,log=TRUE,normalize=FALSE)
  rs <- exp(apply(logW,1,logsumexp))
  D <- diag(sqrt(rs))
  Dinv <- diag(1/sqrt(rs))
  Ms <- Dinv%*%exp(logW)%*%Dinv # symmetrise
  e <- eigen(Ms,symmetric=TRUE)
  s <- sum(sqrt(rs)*e$vec[,1]) # scaling
  # Phi is orthonormal under the weighted inner product
  # t(Phi)%*%diag(1/Phi[,1])%*%Phi == I
  list(Psi=s*Dinv%*%e$vec, Phi=(1/s)*D%*%e$vec, eig=e$val)
}
