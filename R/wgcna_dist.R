#' Topological Overlap Matrix (TOM) Dissimilarity
#'
#' Computes TOM dissimilarity matrix from trait matrix or harmonized object.
#'
#' @param object wide matrix or harmonized long object
#' @param params list of parameters for WGCNA routines
#'
#' @return numeric TOM dissimilarity matrix
#' @export
wgcna_dist <- function(object, params) {

  # Pivot object to have traits in columns and ID in rownames.
  if(!is.matrix(object))
    object <- wgcna_pivot(object)$matrix

  # Check and assign parameters
  params <- wgcna_params(params)

  # Topological overlap (TOM)
  1 - WGCNA::TOMsimilarity(
    WGCNA::adjacency(
      object,
      type = params$signType,
      power = params$power),
    TOMType = params$signType,
    verbose = params$verbose)
}
