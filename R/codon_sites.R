# R/codon_sites.R
# Nei-Gojobori expected synonymous / non-synonymous site counts.
#
# These are the denominators missing from the current pipeline: 09_substitutions.R
# divides raw substitution counts by branch time only, so its "rates" are changes
# per year rather than changes per site per year (Reviewer 1, major comment 2).

#' Build the 64-codon table of expected synonymous / non-synonymous sites.
#'
#' For each codon, each of its 3 nucleotide positions contributes f/3 synonymous
#' sites, where f is the number of the 3 possible single-nucleotide substitutions
#' at that position that leave the encoded amino acid unchanged. Synonymous and
#' non-synonymous sites for a codon therefore sum to 3.
#'
#' Stop codons return NA (they are not scored as coding positions).
#'
#' @return Numeric matrix, 64 rows (codons, uppercase) x 2 columns (S, N).
codon_site_table <- function() {
  gc    <- as.list(Biostrings::GENETIC_CODE)
  bases <- c("A", "C", "G", "T")
  codons <- do.call(paste0, expand.grid(bases, bases, bases,
                                        stringsAsFactors = FALSE))

  one_codon <- function(codon) {
    aa <- gc[[codon]]
    if (is.null(aa) || aa == "*") return(c(S = NA_real_, N = NA_real_))
    S <- 0
    for (i in 1:3) {
      orig <- substr(codon, i, i)
      syn <- sum(vapply(setdiff(bases, orig), function(b) {
        mut <- codon
        substr(mut, i, i) <- b
        aa2 <- gc[[mut]]
        !is.null(aa2) && aa2 == aa
      }, logical(1)))
      S <- S + syn / 3
    }
    c(S = S, N = 3 - S)
  }

  tab <- t(vapply(codons, one_codon, numeric(2)))
  rownames(tab) <- codons
  tab
}

#' Sum expected synonymous / non-synonymous sites over a codon sequence.
#'
#' Codons containing gaps, ambiguity codes, or stop codons are skipped, so the
#' returned counts reflect only positions that could actually be scored for a
#' substitution — this keeps the denominator consistent with the numerator when
#' gappy codons are excluded from counting.
#'
#' @param seqrow    Character vector of nucleotides, length a multiple of 3.
#' @param site_tab  Output of codon_site_table() (pass once, reused across tips).
#' @return Named numeric vector c(S = , N = , n_codons_scored = ).
seq_site_counts <- function(seqrow, site_tab) {
  seqrow <- toupper(as.character(seqrow))
  n <- length(seqrow)
  if (n %% 3 != 0) {
    stop("sequence length (", n, ") is not a multiple of 3 — check reading frame")
  }
  starts <- seq(1, n, by = 3)
  codons <- paste0(seqrow[starts], seqrow[starts + 1], seqrow[starts + 2])
  ok <- codons %in% rownames(site_tab)
  if (!any(ok)) return(c(S = NA_real_, N = NA_real_, n_codons_scored = 0))
  vals <- site_tab[codons[ok], , drop = FALSE]
  keep <- !is.na(vals[, "S"])
  c(S = sum(vals[keep, "S"]),
    N = sum(vals[keep, "N"]),
    n_codons_scored = sum(keep))
}
