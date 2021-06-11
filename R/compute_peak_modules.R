#' Performs clustering analysis on top of a peak intensity matrix
#'
#' @param peak_intensity_matrix matrix of peak intensities. Each column
#' represents intensities of one peak in different samples, each row represent
#' intensitiess of different peaks in one sample
#' @param peak_correlation_matrix correlation matrix of peak intensities
#' @param correlation_threshold see the \emph{cutHeight} parameter of
#' \link{WGCNA::mergeCloseModules}
#' @param deep_split see \link{WGCNA::mergeCloseModules}
#' @param min_cluster_size see \link{WGCNA::blockwiseModules},
#' \link{WGCNA::mergeCloseModules}, and \link{dynamicTreeCut::cutreeDynamic}
#' @param network_type see \link{WGCNA::adjacency}
#'
#' @return data frame with 3 columns: \emph{peak}, \emph{mean_intensity}, and
#' \emph{module}
#' @import WGCNA
compute_peak_modules <- function(
  peak_intensity_matrix,
  peak_correlation_matrix,
  correlation_threshold,
  deep_split,
  min_cluster_size,
  network_type
) {
  peak_correlation_matrix <- 1 - peak_correlation_matrix

  power_estimate <- WGCNA::pickSoftThreshold(peak_intensity_matrix)$powerEstimate
  if (is.na(power_estimate)) {
    warning("Unable to estimate soft threshold. Using fallback value.")
    power_estimate <- 6
  }

  # NOTE: When WGCNA is loaded incorectly then WGCNA::blockwiseModules may fail
  # due to conflicting functions WGCNA::cor and stats::cor.
  modules <- WGCNA::blockwiseModules(
    datExpr = peak_intensity_matrix,
    checkMissingData = FALSE,
    blocks = dynamicTreeCut::cutreeDynamic(
      dendro = flashClust::flashClust(as.dist(peak_correlation_matrix)),
      distM = peak_correlation_matrix,
      minClusterSize = min_cluster_size,
      pamRespectsDendro = FALSE
    ),
    blockSizePenaltyPower = 100,
    power = power_estimate,
    deepSplit = deep_split,
    detectCutHeight = NULL,
    minModuleSize = min_cluster_size,
    pamRespectsDendro = FALSE,
    mergeCutHeight = 1 - correlation_threshold
  )

  if (!modules$MEsOK) {
    warning("Unable to perform peak WGCNA clustering. Using fallback method.")

    distance_matrix <- WGCNA::TOMdist(
      WGCNA::adjacency(
        datExpr = peak_intensity_matrix,
        type = network_type,
        power = power_estimate,
        corOptions = list(use = "p", method = "pearson")
      )
    )

    modules <- WGCNA::mergeCloseModules(
      exprData = peak_intensity_matrix,
      colors = dynamicTreeCut::cutreeDynamic(
        dendro = flashClust::flashClust(as.dist(distance_matrix)),
        distM = distance_matrix,
        deepSplit = deep_split,
        minClusterSize = min_cluster_size,
        pamRespectsDendro = FALSE
      ),
      cutHeight = 1 - correlation_threshold,
      trapErrors = TRUE
    )
  }

  tibble(
    peak = as.integer(colnames(peak_intensity_matrix)),
    mean_intensity = colMeans(peak_intensity_matrix),
    module = as.integer(as.factor(modules$colors)),
  )
}
