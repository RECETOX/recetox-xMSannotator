init_chemscoremat <- function(chemscoremat) {
  if (is.na(chemscoremat)) {
    chemscoremat <- read.csv("Stage4.csv")
  }
  chemscoremat <- as.data.frame(chemscoremat)
  chemscoremat$mz <- as.numeric(chemscoremat$mz)
  return(chemscoremat)
}

increase_multimatches_score <- function(multimatch_features, adduct_weights) {
  features_adducts <- which(multimatch_features$Adduct %in% adduct_weights[, 1])
  if (length(features_adducts) > 0) {
    multimatch_features$score[features_adducts] <- (multimatch_features$score[features_adducts]) * 100
  }
  return(multimatch_features)
}

reevaluate_multimatches_score <- function(single_molecule_annotation) {
  num_annotations <- nrow(single_molecule_annotation)
  score <- (num_annotations - 1) * max(single_molecule_annotation$score) / num_annotations
  return(score)
}

get_features <- function(mz, match) {
  feature_count <- table(mz)
  features <- switch(match,
    "multiple" = feature_count[which(feature_count > 1)],
    "unique" = feature_count[which(feature_count == 1)]
  )
  features <- names(features)
  return(features)
}

remove_tmp_files <- function(loc) {
  fname <- paste("Stage5_annotation_results", sep = "")
  unlink(fname)
  suppressWarnings(unlink("*.Rda"))

  for (file in c("/stage2/", "/stage3/", "/stage4/", "/stage5/")) {
    loc <- paste(loc, file, sep = "")
    unlink(loc, force = TRUE, recursive = TRUE)
  }

  for (file in c("step1_results.Rda", "plot.pdf", "Rplots.pdf", "Rplots.pdf")) {
    try(unlink(file), silent = TRUE)
  }
}

compute_multimatches <- function(mz, curated_res, adduct_weights) {
  multimatches_idx <- which(curated_res$mz %in% mz)
  multimatch_features <- increase_multimatches_score(curated_res[multimatches_idx, ], adduct_weights)
  max_score_idx <- which(multimatch_features$score == max(multimatch_features$score, na.rm = TRUE))
  multimatches_idx <- multimatches_idx[-max_score_idx]

  for (feature in 1:nrow(multimatch_features)) {
    if (multimatch_features$score[feature] != max(multimatch_features$score)) {
      same_molecule_idx <- which(curated_res$chemical_ID %in% multimatch_features$chemical_ID[feature])
      curated_res$score[same_molecule_idx] <- reevaluate_multimatches_score(curated_res[same_molecule_idx, ])
    }
  }
  return(multimatches_idx)
}

multilevelannotationstep5 <- function(outloc,
                                      adduct_weights = NA,
                                      chemscoremat = NA) {
  setwd(outloc)
  curated_res <- init_chemscoremat(chemscoremat)

  if (is.na(adduct_weights)) {
    adduct_weights <- data.frame(
      Adduct = c("M+H", "M-H"),
      Weight = c(1, 1)
    )
  }

  curated_res <- curated_res[order(curated_res$Confidence,
    curated_res$chemical_ID,
    curated_res$score,
    curated_res$Adduct,
    decreasing = TRUE
  ), ]

  duplicated_features <- get_features(curated_res$mz, "multiple")

  lower_score_multimatches_idx <- unlist(lapply(
    duplicated_features,
    compute_multimatches,
    curated_res,
    adduct_weights
  ))
  
  if (length(lower_score_multimatches_idx) > 0) {
    curated_res <- curated_res[-c(lower_score_multimatches_idx), ]
  }

  unique_features <- get_features(curated_res$mz, "unique")
  curated_res$MatchCategory <- rep("Multiple", nrow(curated_res))
  curated_res$MatchCategory[which(curated_res$mz %in% unique_features)] <- "Unique"

  write.csv(curated_res, file = "Stage5.csv", row.names = FALSE)
  remove_tmp_files(outloc)
  return(curated_res)
}
