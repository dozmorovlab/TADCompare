#' Compute the top eigenvalues/eigenvectors of a symmetric matrix
#'
#' Internal helper used to compute the leading eigenvalues and eigenvectors
#' of a symmetric matrix. A base R \code{eigen} decomposition is used, and
#' the eigenpairs corresponding to the \code{NEig} largest-magnitude
#' eigenvalues are returned.
#'
#' @param mat A symmetric matrix.
#' @param NEig Number of eigenvalues/eigenvectors to return. Defaults to 2.
#'
#' @return A list with elements \code{values} (a numeric vector of the
#'   \code{NEig} eigenvalues, largest magnitude first) and \code{vectors}
#'   (a matrix whose columns are the corresponding eigenvectors).
#'
#' @noRd

get_eigs = function(mat, NEig = 2) {

  full_eigen = eigen(mat, symmetric = TRUE)
  NEig = min(NEig, length(full_eigen$values))
  ord = order(-abs(full_eigen$values))[seq_len(NEig)]
  Eigen = list(values = full_eigen$values[ord],
               vectors = full_eigen$vectors[, ord, drop = FALSE])

  return(Eigen)

}
