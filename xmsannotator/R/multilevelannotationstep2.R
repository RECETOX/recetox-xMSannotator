#' @export
replace_with_module <- function(module_rt_clust) {
  gsub(as.character(module_rt_clust), pattern = "_[0-9]*", replacement = "")
}

#' @export
compute_something <- function(m, curmchemdata, mass_defect_window, isop_res_md) {
  # Get Module ID
  module_rt_group <- replace_with_module(curmchemdata$Module_RTclust[m])

  # Get mass defect
  query_md <- curmchemdata$mz[m] - round(curmchemdata$mz[m])

  md_diff_smaller_window <- abs((isop_res_md$MD) - (query_md)) < mass_defect_window
  is_in_module <- isop_res_md$Module_RTclust == module_rt_group
  extract_indices <- which(md_diff_smaller_window & is_in_module)

  as.data.frame(isop_res_md[extract_indices, ])
}

#' @import plyr
#' @export
compute_chemscore <- function(j,
                              chemids,
                              mchemdata,
                              isop_res_md,
                              mass_defect_window,
                              corthresh,
                              global_cor,
                              mzid,
                              max_diff_rt,
                              adduct_table,
                              adduct_weights,
                              filter.by,
                              max_isp,
                              MplusH.abundance.ratio.check,
                              mass_defect_mode,
                              outloc) {
  chemid <- chemids[j]
  chemscoremat <- {}
  curmchemdata <- mchemdata[which(mchemdata$chemical_ID == chemid), ]
  curmchemdata$mz <- as.numeric(as.character(curmchemdata$mz))
  curmchemdata$time <- as.numeric(as.character(curmchemdata$time))
  curmchemdata <- as.data.frame(curmchemdata)
  curmchemdata$Module_RTclust <- replace_with_module(curmchemdata$Module_RTclust)
  isop_res_md$Module_RTclust <- replace_with_module(isop_res_md$Module_RTclust)

  isp_masses_mz_data <- {}
  isp_masses_mz_data <- lapply(
    1:length(curmchemdata$mz),
    compute_something,
    curmchemdata,
    mass_defect_window,
    isop_res_md
  )
  isp_masses_mz_data <- plyr::ldply(isp_masses_mz_data, rbind)
  colnames(isp_masses_mz_data)[5] <- "AvgIntensity"
  isp_masses_mz_data <- as.data.frame(isp_masses_mz_data)
  isp_masses_mz_data$mz <- as.numeric(as.character(isp_masses_mz_data$mz))
  isp_masses_mz_data$time <- as.numeric(as.character(isp_masses_mz_data$time))
  isp_masses_mz_data <- as.data.frame(isp_masses_mz_data)


  if (is.na(mass_defect_mode) == TRUE) {
    mass_defect_mode <- "pos"
  }

  chem_score <- get_chemscorev1.6.71(
    chemicalid = chemid,
    mchemicaldata = curmchemdata,
    corthresh = corthresh,
    global_cor = global_cor,
    mzid = mzid,
    max_diff_rt = max_diff_rt,
    level_module_isop_annot = isp_masses_mz_data,
    adduct_table = adduct_table,
    adduct_weights = adduct_weights,
    filter.by = filter.by,
    max_isp = max_isp,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
    mass_defect_window = mass_defect_window,
    mass_defect_mode = mass_defect_mode,
    outlocorig = outloc
  )

  if (length(chem_score) > 0) {
    if (chem_score$chemical_score >= (-100)) {
      chem_score$filtdata <- chem_score$filtdata[order(chem_score$filtdata$mz), ]

      if (nrow(chem_score$filtdata) > 0) {
        cur_chem_score <- chem_score$chemical_score
        chemscoremat <- cbind(cur_chem_score, chem_score$filtdata)
      }
      rm(chem_score)
      chemscoremat <- as.data.frame(na.omit(chemscoremat))
    }
  } else {
    rm(chem_score)
  }
  rm("curmchemdata", "isp_masses_mz_data", "chemid")

  return(chemscoremat)
}

#' @import plyr
#' @export
multilevelannotationstep2 <- function(outloc1,
                                      list_number,
                                      max.rt.diff,
                                      chemids_split,
                                      num_sets,
                                      mchemdata,
                                      mass_defect_window,
                                      corthresh,
                                      global_cor,
                                      mzid,
                                      adduct_table,
                                      adduct_weights,
                                      max_isp,
                                      MplusH.abundance.ratio.check,
                                      mass_defect_mode,
                                      chemids,
                                      isop_res_md,
                                      filter.by) {
  setwd(outloc1)

  # load("step1_results.Rda")
  # load("global_cor.Rda")
  # unlink("allmatches_with_isotopes.txt")
  ls()

  outloc <- outloc1
  if (is.na(max.rt.diff) == FALSE) {
    max_diff_rt <- max.rt.diff
  }

  if (list_number > length(chemids_split) | list_number > num_sets) {
    list_number <- length(chemids_split)
    return(0)
  }

  if (anyNA(adduct_weights)) {
    data(adduct_weights)

    adduct_weights1 <- matrix(nrow = 2, ncol = 2, 0)
    adduct_weights1[1, ] <- c("M+H", 1)
    adduct_weights1[2, ] <- c("M-H", 1)
    adduct_weights <- as.data.frame(adduct_weights1)
    colnames(adduct_weights) <- c("Adduct", "Weight")
  }

  outloc1 <- paste(outloc, "/stage2/", sep = "")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1)


  cnames <- colnames(mchemdata)
  cnames[2] <- "time"
  colnames(mchemdata) <- as.character(cnames)
  mchemdata$mz <- as.numeric(as.character(mchemdata$mz))
  mchemdata$time <- as.numeric(as.character(mchemdata$time))

  chem_score <- lapply(
    chemids_split[[list_number]],
    compute_chemscore,
    chemids = chemids,
    mchemdata = mchemdata,
    isop_res_md = isop_res_md,
    mass_defect_window = mass_defect_window,
    corthresh = corthresh,
    global_cor = global_cor,
    mzid = mzid,
    max_diff_rt = max_diff_rt,
    adduct_table = adduct_table,
    adduct_weights = adduct_weights,
    filter.by = filter.by,
    max_isp = max_isp,
    MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
    mass_defect_mode = mass_defect_mode,
    outloc = outloc
  )

  chem_score2 <- chem_score[which(chem_score != "NULL")]

  rm(chem_score)

  curchemscoremat <- plyr::ldply(chem_score2, rbind)
  return(curchemscoremat)
}

#' @export
call_multilevelannotationstep2 <- function(arg1,
                                           outloc,
                                           max.rt.diff,
                                           chemids_split,
                                           num_sets,
                                           mchemdata,
                                           mass_defect_window,
                                           corthresh,
                                           global_cor,
                                           mzid,
                                           adduct_table,
                                           adduct_weights,
                                           max_isp,
                                           MplusH.abundance.ratio.check,
                                           mass_defect_mode,
                                           chemids,
                                           isop_res_md,
                                           filter.by) {
  cur_fname <- paste(outloc, "/stage2/chem_score", arg1, ".Rda", sep = "")
  check_if_exists <- file.exists(cur_fname)


  if (check_if_exists == FALSE) {
    multilevelannotationstep2(
      outloc1 = outloc,
      list_number = arg1,
      max.rt.diff = max.rt.diff,
      chemids_split = chemids_split,
      num_sets = num_sets,
      mchemdata = mchemdata,
      mass_defect_window = mass_defect_window,
      corthresh = corthresh,
      global_cor = global_cor,
      mzid = mzid,
      adduct_table = adduct_table,
      adduct_weights = adduct_weights,
      max_isp = max_isp,
      MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
      mass_defect_mode = mass_defect_mode,
      chemids = chemids,
      isop_res_md = isop_res_md,
      filter.by = filter.by
    )
  } else {
    print(paste("List ", arg1, " already exists.", sep = ""))
  }
}