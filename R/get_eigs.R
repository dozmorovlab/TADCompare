#' Compute the top eigenvalues/eigenvectors of a symmetric matrix
#'
#' Internal helper used to compute the leading eigenvalues and eigenvectors
#' of a symmetric matrix. If the \code{PRIMME} package is installed, it is
#' used via \code{PRIMME::eigs_sym} for fast, iterative eigendecomposition of
#' the largest-magnitude eigenpairs. If \code{PRIMME} is not available, a
#' base R \code{eigen} decomposition is used instead, and the eigenpairs
#' corresponding to the \code{NEig} largest-magnitude eigenvalues are
#' returned so that the output matches \code{PRIMME::eigs_sym}'s default
#' behavior (\code{which = "LM"}).
#'
#' @param mat A symmetric matrix.
#' @param NEig Number of eigenvalues/eigenvectors to return. Defaults to 2.
#'
#' @return A list with elements \code{values} (a numeric vector of the
#'   \code{NEig} eigenvalues, largest magnitude first) and \code{vectors}
#'   (a matrix whose columns are the corresponding eigenvectors), mirroring
#'   the structure returned by \code{PRIMME::eigs_sym}.
#'
#' @noRd

get_eigs = function(mat, NEig = 2) {

  if (requireNamespace("PRIMME", quietly = TRUE)) {

    Eigen = PRIMME::eigs_sym(mat, NEig = NEig)

  } else {

    # Fall back to base R eigendecomposition when PRIMME is unavailable

    full_eigen = eigen(mat, symmetric = TRUE)

    # PRIMME::eigs_sym defaults to the NEig largest-magnitude eigenvalues

    NEig = min(NEig, length(full_eigen$values))

    ord = order(-abs(full_eigen$values))[seq_len(NEig)]

    Eigen = list(values = full_eigen$values[ord],
                 vectors = full_eigen$vectors[, ord, drop = FALSE])

  }

  return(Eigen)

}
