# R/matrix_utils.R
# Utilities for assembling and validating pairwise dN/dS distance matrices.

#' Read a PAML pairwise distance matrix (2ML.dN or 2ML.dS format)
#'
#' @param path Path to the PAML matrix file.
#' @return A named numeric matrix.
read_paml_matrix <- function(path) {
  if (!file.exists(path)) stop("PAML matrix file not found: ", path)

  lines <- readLines(path)
  lines <- lines[nchar(trimws(lines)) > 0]

  n <- as.integer(trimws(lines[1]))
  if (is.na(n) || n < 2) stop("Could not parse matrix size from: ", path)

  mat <- matrix(0.0, nrow = n, ncol = n)
  labels <- character(n)

  for (i in seq_len(n)) {
    row_line <- trimws(lines[i + 1])
    parts    <- strsplit(row_line, "\\s+")[[1]]
    labels[i] <- parts[1]
    vals       <- as.numeric(parts[-1])
    if (length(vals) > 0) mat[i, seq_along(vals)] <- vals
  }

  # PAML lower-triangle format: mirror to symmetric
  mat <- mat + t(mat)

  rownames(mat) <- labels
  colnames(mat) <- labels
  mat
}


#' Stitch a list of pairwise matrices into one full distance matrix.
#'
#' Overlapping entries are averaged; missing pairs remain NA.
#'
#' @param matrices  List of named numeric matrices (output of read_paml_matrix).
#' @param expected_labels  Character vector of all expected taxon labels.
#' @return A symmetric named numeric matrix covering expected_labels.
matrix_stitch <- function(matrices, expected_labels) {
  n      <- length(expected_labels)
  full   <- matrix(NA_real_, nrow = n, ncol = n,
                   dimnames = list(expected_labels, expected_labels))
  counts <- matrix(0L,       nrow = n, ncol = n,
                   dimnames = list(expected_labels, expected_labels))

  for (mat in matrices) {
    if (is.null(mat)) next
    shared <- intersect(rownames(mat), expected_labels)
    if (length(shared) < 2) next
    sub <- mat[shared, shared, drop = FALSE]
    full[shared, shared]   <- rowSums(cbind(
      ifelse(is.na(full[shared, shared]), 0, full[shared, shared]),
      sub
    ), na.rm = FALSE)  # manual accumulation below
    counts[shared, shared] <- counts[shared, shared] + 1L
  }

  # Proper average accumulation
  full <- matrix(NA_real_, nrow = n, ncol = n,
                 dimnames = list(expected_labels, expected_labels))
  count <- matrix(0L, nrow = n, ncol = n,
                  dimnames = list(expected_labels, expected_labels))

  for (mat in matrices) {
    if (is.null(mat)) next
    shared <- intersect(rownames(mat), expected_labels)
    if (length(shared) < 2) next
    sub <- mat[shared, shared, drop = FALSE]
    cur <- full[shared, shared]
    cur[is.na(cur)] <- 0
    full[shared, shared]  <- cur + sub
    count[shared, shared] <- count[shared, shared] + 1L
  }

  # Average where we have values
  full <- full / ifelse(count == 0L, NA_real_, count)

  diag(full) <- 0
  full
}


#' Validate that a distance matrix is square, symmetric, and the right size.
#'
#' @param mat      Numeric matrix to validate.
#' @param expected_n  Expected number of rows/columns.
#' @return Invisibly returns mat; stops with error on failure.
validate_matrix <- function(mat, expected_n) {
  if (!is.matrix(mat))             stop("validate_matrix: input is not a matrix")
  if (nrow(mat) != ncol(mat))      stop("validate_matrix: matrix is not square")
  if (nrow(mat) != expected_n)     stop("validate_matrix: expected ", expected_n,
                                        " taxa, got ", nrow(mat))
  if (!isSymmetric(mat, tol = 1e-8, check.attributes = FALSE))
    warning("validate_matrix: matrix is not symmetric (tolerance 1e-8)")
  if (any(diag(mat) != 0, na.rm = TRUE))
    warning("validate_matrix: diagonal is not all zero")

  na_frac <- mean(is.na(mat[upper.tri(mat)]))
  if (na_frac > 0.5) warning("validate_matrix: >50% of upper-triangle entries are NA")

  invisible(mat)
}
