#' Optionaly exlude some entries from data frame (fancier dplyr::anti_join)
#'
#' @param data a data frame
#' @param exlude optional vector of values to exclude
#' @param column_name the name of the column where to search for unwanted values
#'
#' @return data frame with excluded entries
optional_exclude <- function(data, exclude, column_name) {
  if (is.null(exclude)) {
    return(data)
  }

  if (is.vector(exclude)) {
    exclude <- tibble::as_tibble_col(exclude, column_name = column_name)
  }

  anti_join(data, exclude, by = column_name)
}


#' Remove pathways which are statisticaly independent on significant annotations
#'
#' Test whether compounds associated with some specific pathway are equaly
#' likely to be marked as significant as compounds not assotiated with the
#' pathway. If \empth{x} is a compound, \emph{P} is some specific pathway and
#' \emph{S} is a group of significant annotations then the data can be
#' represented by the following contingency table:
#'
#' |            | x \in P | x \notin P |
#' |------------|---------|------------|
#' | x \in S    |    a    |      b     |
#' | x \notin S |    c    |      d     |
#'
#' @param pathways data frame of pathway-compound mappings
#' @param significant vector of significant compounds
#'
#' @return currated data frame of pathway-compounds mappings
remove_indipendent_pathways <- function (pathways, significant) {
  n_all_significant_compounds <- n_distinct(significant)
  n_all_compounds_in_pathways <- n_distinct(pathways$compound)

  pathways <- distinct(pathways, .data$pathway, .data$compound)
  pathways <- mutate(pathways, significant = .data$compound %in% !!significant)

  pathways <- with_groups(pathways, .data$pathway, mutate, p_value = {
    a <- sum(.data$significant)
    c <- sum(!.data$significant)
    b <- !!n_all_significant_compounds - a
    d <- !!n_all_compounds_in_pathways - a - b - c

    fisher.test(matrix(c(a, c, b, d), nrow = 2))$p.value
  })

  pathways <- filter(pathways, .data$p_value <= 0.05, .data$significant)
  pathways <- select(pathways, .data$pathway, .data$compound)
  pathways
}


#' For each compound computes most abundant modules
#'
#' For each annotated compound computes a set (ie. more then one) of most
#' abundant (dominant) modules. The abundance of modules is computed from
#' distinct adduct-compound-peak tripplets.
#'
#' @param annotations data frame with annotations
#'
#' @return tibble with \emph{compound}, \emph{dominant_module} columns
compute_dominant_modules <- function(annotations) {
  x <- as_tibble(annotations)
  x <- distinct(x, .data$adduct, .data$compound, .data$peak, .keep_all = TRUE)

  x <- count(x, .data$compound, .data$module)
  x <- with_groups(x, .data$compound, filter, .data$n == max(.data$n))
  x <- select(x, 'compound', dominant_module = 'module')
  x
}


#' For each pathway compute a score
#'
#' @param pathways
#'
#' @return pathway table with a new pathway_score column
compute_pathway_scores <- function (pathways) {
  # TODO: rewrite line 67-144 of pathway_hmdb file
  mutate(pathways, pathway_score = 0)
}


#' Update scores of annotations with pathway scores
#'
#' @param annotations data frame with feature annotations with scores
#' @param pathways pathway mappings with computed pathway scores
#'
#' @return augmented annotations
update_pathway_scores <- function (annotations, pathways) {
  pathways <- group_by(pathways, .data$compound)
  pathways <- summarise(pathways, pathway_score = sum(.data$pathway_score))

  annotations <- left_join(annotations, pathways, by = 'compound')
  annotations <- replace_na(annotations, list('pathway_score' = 0))
  annotations <- mutate(annotations, score = .data$score + .data$pathway_score)
  annotations <- select(annotations, -.data$pathway_score)
  annotations
}


#' Based on pathway data adjust the scores.
#'
#' @param annotations
#' @param pathways
#' @param exluded_pathways
#' @param exluded_pathway_compounds
#' @param score_threshold
#'
#' @return annotations
compute_pathways <- function(
  annotations,
  pathways,
  exluded_pathways = NULL,
  exluded_pathway_compounds = NULL,
  score_threshold = 0.1
) {
  significant_annotations <- filter(annotations, .data$score >= !!score_threshold, !is.na(.data$adduct_weight))
  significant_annotations <- semi_join(significant_annotations, adduct_weights, by = 'adduct')
  significant_compounds <- unique(significant_annotations$compound)

  pathways <- optional_exclude(pathways, exluded_pathways, 'pathway')
  pathways <- optional_exclude(pathways, exluded_pathway_compounds, 'compound')

  pathways <- remove_indipendent_pathways(pathways, significant_compounds)
  pathways <- compute_pathway_scores(pathways)

  annotations <- update_pathway_scores(annotations, pathways)
  annotations <- filter(annotations, .data$score >= !!score_threshold)
  annotations
}


# function (annotations, pathways) {
#   modules <- annotations %>%
#     distinct(.data$peak, .data$adduct, .data$compound, .data$module) %>%
#     count(.data$compound, .data$module)
#   dominant_modules <- modules %>%
#     with_groups(.data$compound, filter, .data$n == max(.data$n)) %>%
#     select(.data$compound, dominant_module = .data$module)
#
#   dominant_modules %>%
#     mutate(module_cardinality, sapply(.data$dominant_module), function (m)
#       n_distinct(modules$compound[modules$module %in% m])
#     )
#
#   pathways <- with_groups(pathways, .data$pathway, mutate, {
#       pathway_compounds -> cur_data()$compound
#
#       u <- is.element(modules$module, dominant_module)
#       v <- is.element(modules$compound, pathway_compounds)
#
#       a <- n_distinct(modules$compound[ u &  v])
#       b <- n_distinct(modules$compound[ u & !v])
#       c <- n_distinct(modules$compound[!u &  v])
#       d <- n_distinct(modules$compound[!u & !v])
#
#       p_score <- sum(modules$n[u])
#       p_value <- fisher.test(matrix(c(a, c, b, d, 2, 2)))$p.value
#
#       tibble(a = a, b = b, c = c, d = d, p_value = p_value)
#   })
# }
