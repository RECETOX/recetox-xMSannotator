#' @export
do_something <- function(i,
                         mchemicaldata_goodadducts_index,
                         mchemicaldata,
                         level_module_isop_annot,
                         max_diff_rt,
                         mass_defect_window,
                         mass_defect_mode,
                         curformula,
                         exp_isp,
                         max_isp,
                         abund_ratio_vec) {
  m <- mchemicaldata_goodadducts_index[i]
  final_isp_annot_res <- cbind(paste("group", i, sep = ""), mchemicaldata[m, ])

  module_rt_group <- replace_with_module(mchemicaldata$Module_RTclust[m])

  isp_mat_module_rt_group <- as.character(level_module_isop_annot$Module_RTclust)

  query_md <- mchemicaldata$mz[m] - round(mchemicaldata$mz[m])
  query_rt <- mchemicaldata$time[m]
  query_int <- 1 * (mchemicaldata$mean_int_vec[m])

  put_isp_masses_curmz_data <- level_module_isop_annot[which(abs(level_module_isop_annot$time - query_rt) < max_diff_rt & abs((level_module_isop_annot$MD) - (query_md)) < mass_defect_window & isp_mat_module_rt_group == module_rt_group & level_module_isop_annot$AvgIntensity < query_int), ]
  put_isp_masses_curmz_data <- as.data.frame(put_isp_masses_curmz_data)
  put_isp_masses_curmz_data$mz <- as.numeric(as.character(put_isp_masses_curmz_data$mz))
  put_isp_masses_curmz_data$time <- as.numeric(as.character(put_isp_masses_curmz_data$time))
  mchemicaldata <- as.data.frame(mchemicaldata)

  put_isp_masses_curmz_data <- unique(put_isp_masses_curmz_data)
  put_isp_masses_curmz_data <- put_isp_masses_curmz_data[order(put_isp_masses_curmz_data$mz), ]

  if (length(put_isp_masses_curmz_data) < 1) {
    return(final_isp_annot_res)
  }
  int_vec <- put_isp_masses_curmz_data[, c(5)]
  int_vec <- int_vec / mchemicaldata[m, 13]

  curformula <- gsub(curformula, pattern = "Sr", replacement = "")
  curformula <- gsub(curformula, pattern = "Sn", replacement = "")
  curformula <- gsub(curformula, pattern = "Se", replacement = "")
  curformula <- gsub(curformula, pattern = "Sc", replacement = "")
  curformula <- gsub(curformula, pattern = "Sm", replacement = "")

  max_isp_count <- max(exp_isp)

  if (is.na(max_isp_count)) {
    max_isp_count <- 1
  }

  ischeck <- which(int_vec <= max(abund_ratio_vec[-c(1)] + 0.10))

  put_isp_masses_curmz_data[, 2] <- as.numeric(as.character(put_isp_masses_curmz_data[, 2]))

  if (length(ischeck) < 1) {
    return(final_isp_annot_res)
  }

  for (rnum in 1:length(ischeck)) {
    temp_var <- {}
    bool_check <- 1

    isp_v <- ischeck[rnum]
    isp_v <- as.numeric(as.character(isp_v))

    diff_rt <- abs(put_isp_masses_curmz_data[isp_v, 2] - mchemicaldata[m, 2])
    isnum <- (round(put_isp_masses_curmz_data[isp_v, 1]) - round(mchemicaldata[m, 1]))
    bool_check <- 1

    if (mass_defect_mode == "neg" | mass_defect_mode == "both") {
      isnum <- abs(isnum)
    } else {
      if (isnum > 0) {
        bool_check <- 1
      } else {
        bool_check <- 0
      }
    }

    if (diff_rt < max_diff_rt & isnum <= max_isp) {
      max_isp_count <- max_isp

      if (max_isp_count > 0 && isnum <= max_isp_count && bool_check > 0) {
        isnum2 <- (round(put_isp_masses_curmz_data[isp_v, 1]) - round(as.numeric(as.character(mchemicaldata$MonoisotopicMass[m]))))
        isnum2 <- round(isnum2)

        isp_sign <- if (isnum2 <= 0) "-" else "+"

        isnum2 <- abs(isnum2)
        if (isnum2 <= max_isp) {
          form_name <- as.character(paste(mchemicaldata[m, 7], "_[", isp_sign, (isnum2), "]", sep = ""))
          other_inf <- cbind(rep("-", 7))
          temp_var <- cbind(put_isp_masses_curmz_data[isp_v, c(1:2)], t(other_inf), put_isp_masses_curmz_data[isp_v, c(3:4)], put_isp_masses_curmz_data[isp_v, c(2, 5:6)])

          temp_var <- as.data.frame(temp_var)

          colnames(temp_var) <- colnames(mchemicaldata)

          temp_var$Formula <- form_name
          temp_var$Name <- as.character(mchemicaldata[m, 6])
          temp_var$chemical_ID <- as.character(mchemicaldata[m, 5])

          temp_var$Adduct <- paste(mchemicaldata[m, 9], "_[", isp_sign, (abs(isnum)), "]", sep = "")

          temp_var <- as.data.frame(temp_var)
          temp_var <- cbind(paste("group", i, sep = ""), temp_var)
          final_isp_annot_res <- as.data.frame(final_isp_annot_res)

          if (nrow(temp_var) > 0) {
            check_mz <- which(temp_var$mz %in% final_isp_annot_res)

            if (length(check_mz) > 0) {
              temp_var <- temp_var[-c(check_mz), ]
            }
            if (nrow(temp_var) > 0) {
              final_isp_annot_res <- rbind(final_isp_annot_res, temp_var)
            }
          }
        }
      }
    }
  }
  return(final_isp_annot_res)
}

#' @export
do_something_2 <- function(i,
                           mod_names,
                           mchemicaldata) {
  groupA_num <- mod_names[i]

  subdata <- mchemicaldata[which(mchemicaldata$Module_RTclust == groupA_num), ]
  subdata <- subdata[order(subdata$time), ]

  if (nrow(subdata) > 0) {
    groupB <- group_by_rt_histv2(subdata, time_step = 1, max_diff_rt = 10, groupnum = groupA_num)
  } else {
    groupB <- subdata
  }
  rownames(groupB) <- NULL
  return(groupB)
}

#' @export
compute_score <- function(adduct_weights, cur_adducts_with_isotopes) {
  score <- 1 * (10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2]))))
  return(score[1])
}

#' @export
remove_water_adducts <- function(curformula, mchemicaldata) {
  numoxygen <- check_element(curformula, "O")
  water_adducts <- c("M+H-H2O", "M+H-2H2O", "M-H2O-H")
  water_adduct_ind <- which(mchemicaldata$Adduct %in% water_adducts)

  if (numoxygen < 1) {
    if (length(water_adduct_ind) > 0) {
      return(mchemicaldata[-water_adduct_ind, ])
    }
  }

  return(mchemicaldata)
}

#' @export
#' @importFrom Rdisop getMolecule
add_isotopic_peaks <- function(mchemicaldata,
                               adduct_weights,
                               exp_isp,
                               level_module_isop_annot,
                               max_diff_rt,
                               mass_defect_window,
                               mass_defect_mode,
                               max_isp,
                               abund_ratio_vec) {
  mchemicaldata$Module_RTclust <- replace_with_module(mchemicaldata$Module_RTclust)
  mchemicaldata_orig <- mchemicaldata

  adduct_weights_strings <- as.character(adduct_weights[, 1])
  selected_adduct_weights <- which(mchemicaldata_orig$Adduct %in% adduct_weights_strings)

  if (length(selected_adduct_weights) > 0) {
    rel_adduct_data <- mchemicaldata[selected_adduct_weights, ]

    rel_adduct_module <- gsub(rel_adduct_data$Module_RTclust, pattern = "_[0-9]*", replacement = "")
    module_rt_group <- gsub(mchemicaldata$Module_RTclust, pattern = "_[0-9]*", replacement = "")

    mchemicaldata <- mchemicaldata[which(module_rt_group %in% rel_adduct_module), ]
  }

  mchemicaldata <- unique(mchemicaldata)

  curformula <- as.character(mchemicaldata[, 7])

  formula_check <- Rdisop::getMolecule(as.character(curformula[1]))
  exp_isp <- which(formula_check$isotopes[[1]][2, ] >= 0.001)
  abund_ratio_vec <- formula_check$isotopes[[1]][2, exp_isp]

  mchemicaldata <- remove_water_adducts(curformula, mchemicaldata)

  # water check
  if (nrow(mchemicaldata) < 1) {
    chemical_score <- (-100)
    return(list("chemical_score" = chemical_score, "filtdata" = mchemicaldata))
  }

  mchemicaldata$Adduct <- as.character(mchemicaldata$Adduct)

  final_isp_annot_res_all <- mchemicaldata
  level_module_isop_annot <- as.data.frame(level_module_isop_annot)
  level_module_isop_annot$Module_RTclust <- gsub(level_module_isop_annot$Module_RTclust, pattern = "_[0-9]*", replacement = "")

  mchemicaldata_goodadducts_index <- which(mchemicaldata$Adduct %in% as.character(adduct_weights[, 1]))
  final_isp_annot_res_isp <- {}

  if (length(mchemicaldata_goodadducts_index) > 0) {
    final_isp_annot_res_isp <- lapply(
      1:length(mchemicaldata_goodadducts_index),
      do_something,
      mchemicaldata_goodadducts_index,
      mchemicaldata,
      level_module_isop_annot,
      max_diff_rt,
      mass_defect_window,
      mass_defect_mode,
      curformula,
      exp_isp,
      max_isp,
      abund_ratio_vec
    )

    rm(level_module_isop_annot)

    final_isp_annot_res2 <- plyr::ldply(final_isp_annot_res_isp, rbind)[, -c(1)]

    rm(final_isp_annot_res_isp)
    mchemicaldata <- rbind(final_isp_annot_res_all, final_isp_annot_res2) # [,-c(12)]
  }

  # filter for na and sort in ascending mz
  mchemicaldata <- unique(mchemicaldata) %>%
    drop_na(mz) %>%
    arrange(mz)
  mod_names <- unique(mchemicaldata$Module_RTclust)

  diffmatB <- {}
  diffmatB <- lapply(1:length(mod_names), do_something_2, mod_names, mchemicaldata)

  mchemicaldata <- unique(plyr::ldply(diffmatB, rbind))

  rm(diffmatB)

  write.table(mchemicaldata, file = "../Stage2_withisotopes.txt", append = TRUE, sep = "\t", col.names = FALSE)
  return(mchemicaldata)
}

compute_hub_mz_list <- function(check_cor, cur_adducts_with_isotopes, filter.by, check2, adduct_weights) {
  hub_mz_list <- which(check_cor > 0 & (cur_adducts_with_isotopes %in% filter.by == TRUE))

  if (length(hub_mz_list) < 1) {
    hub_mz_list <- which(check_cor > 0 & check2 < 0 & (cur_adducts_with_isotopes %in% adduct_weights[, 1] == TRUE))
  }

  if (length(hub_mz_list) < 1) {
    hub_mz_list <- which(check_cor > 0 & check2 < 0)
  }
  return(hub_mz_list)
}

compute_diff_rt_hubmz <- function(mchemicaldata, hub_mz) {
  diff_rt_hubmz <- apply(mchemicaldata, 1, function(k) {
    curtime <- as.numeric(as.character(k[3]))
    return(abs(mchemicaldata$time[hub_mz] - curtime))
  })
  return(diff_rt_hubmz)
}

compute_hub_mz <- function(hub_mz_list, mchemicaldata, max_diff_rt) {
  hub_mz <- hub_mz_list[which(mchemicaldata$mean_int_vec[hub_mz_list] == max(mchemicaldata$mean_int_vec[hub_mz_list]))[1]]

  max_time_neighbors <- 0

  for (h1 in hub_mz_list) {
    diff_rt_hubmz <- compute_diff_rt_hubmz(mchemicaldata, h1)
    num_time_neighbors <- length(which(diff_rt_hubmz <= max_diff_rt))

    if (num_time_neighbors > max_time_neighbors) {
      hub_mz <- h1
      max_time_neighbors <- num_time_neighbors
    }
  }
  return(hub_mz)
}

calc_base_score <- function(cur_adducts_with_isotopes, adduct_weights, topquant_cor) {
  cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])", replacement = "")
  good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))
  chemical_score <- length(unique(cur_adducts)) * good_adducts_len * (1 * (topquant_cor))

  if (good_adducts_len > 0) {
    chemical_score <- sum(chemical_score * (10^max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))))
    chemical_score <- chemical_score[1]
  }

  check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[1-2]*\\])")

  if (length(which(check2 > 0)) > 0) {
    chemical_score <- 100 * chemical_score
  }
  return(chemical_score)
}

calc_base_score_v2 <- function(cur_adducts_with_isotopes, adduct_weights, topquant_cor) {
  cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")

  good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))
  chemical_score <- length(unique(cur_adducts)) * good_adducts_len * (1 * (topquant_cor))

  if (good_adducts_len > 0) {
    chemical_score <- sum(chemical_score * (as.numeric(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts), 2])))
  }
  return(chemical_score)
}

compute_diff_rt <- function(time) {
  diff_rt <- abs(min(as.numeric(time)) - max(as.numeric(time)))
  return(diff_rt)
}

compute_temp_score_and_data <- function(group_sizes,
                                        good_temp,
                                        mchemicaldata,
                                        time_cor_groups,
                                        max_diff_rt,
                                        adduct_weights,
                                        topquant_cor,
                                        chemicalid) {
  temp_best_score <- (-100)
  temp_best_data <- {}

  for (g2 in 1:length(group_sizes)) {
    if (g2 %in% good_temp) {
      mchemicaldata <- {}
      mchemicaldata <- as.data.frame(rbind(mchemicaldata, time_cor_groups[[g2]]))

      if (length(mchemicaldata) < 1) {
        next
      }

      diff_rt <- round(compute_diff_rt(mchemicaldata$time))

      if (diff_rt <= max_diff_rt) {
        if (dim(mchemicaldata)[1] > 1) {
          chemical_score <- calc_base_score(mchemicaldata$Adduct, adduct_weights, topquant_cor)
        } else {
          chemical_score <- 0
        }
        names(chemical_score) <- chemicalid[1]
      } else {
        d1 <- density(mchemicaldata$time, bw = max_diff_rt, from = min(mchemicaldata$time) - 0.001, to = (0.01 + max(mchemicaldata$time)), na.rm = TRUE)
        s1 <- summary(d1$x)
        iqr1 <- s1[5] - s1[2]

        if (iqr1 > max_diff_rt / 2) {
          iqr1 <- max_diff_rt / 2
        }
        min_val <- s1[2]
        max_val <- s1[5]

        if (length(which(mchemicaldata$time >= min_val & mchemicaldata$time <= max_val)) > 1) {
          mchemicaldata <- mchemicaldata[which(mchemicaldata$time >= (min_val - 1) & mchemicaldata$time <= (max_val - 1)), ]
          if (dim(mchemicaldata)[1] > 1) {
            chemical_score <- calc_base_score(mchemicaldata$Adduct, adduct_weights, topquant_cor)
          }
        } else {
          chemical_score <- calc_base_score_v2(mchemicaldata$Adduct, adduct_weights, topquant_cor)
        }
        names(chemical_score) <- chemicalid[1]
      }

      if (chemical_score > temp_best_score) {
        temp_best_data <- mchemicaldata
        mchemicaldata <- temp_best_data
        temp_best_score <- chemical_score
      }
    }
  }
  return(list("data" = temp_best_data, "score" = temp_best_score))
}

compute_time_cor_groups <- function(mchemicaldata,
                                max_diff_rt,
                                min_val,
                                max_val,
                                iqr1) {
  diff_rt <- compute_diff_rt(mchemicaldata$time)
  breaks <- {}

  if (diff_rt > 2 * max_diff_rt) {
    breaks <- seq(min_val - max_diff_rt, max_val + max_diff_rt, iqr1)
  } else {
    max_val <- max(mchemicaldata$time)
    min_val <- min(mchemicaldata$time)

    if (min_val < max_diff_rt) {
      breaks <- c(0, max_val + 1)
    } else {
      if (diff_rt < max_diff_rt) {
        breaks <- c(min_val - diff_rt, max_val + diff_rt, diff_rt * 2)
      } else {
        breaks <- seq(min_val - diff_rt, max_val + diff_rt, 1 * max_diff_rt)
      }
    }
  }
  time_cor_groups <- sapply(list(myData1 = mchemicaldata), function(x) split(x, cut(mchemicaldata$time, breaks = breaks)))
  return(time_cor_groups)
}

compute_iqr1 <- function(s1) {
  return(min(abs(s1[3] - s1[2]), abs(s1[3] - s1[5])))
}

compute_min_max_iqr_basic <- function(mchemicaldata, max_diff_rt) {
  d1 <- density(mchemicaldata$time, bw = max_diff_rt, from = min(mchemicaldata$time) - 0.001, to = (0.01 + max(mchemicaldata$time)), na.rm = TRUE)
  s1 <- summary(d1$x)
  iqr1 <- IQR(d1$x)

  if (iqr1 > max_diff_rt / 2) {
    iqr1 <- max_diff_rt / 2
  }
  min_val <- s1[2]
  max_val <- s1[5]

  return(list("min_val" = min_val, "max_val" = max_val, "iqr1" = iqr1))
}

compute_min_max_iqr_advanced <- function(mchemicaldata, max_diff_rt) {
  s1 <- summary(mchemicaldata$time)
  iqr1 <- IQR(mchemicaldata$time)
  min_val <- s1[2] - (1.5 * iqr1)
  max_val <- s1[5] + (1.5 * iqr1)

  if (min_val < s1[1] && max_val > s1[6]) {
    iqr1 <- compute_iqr1(s1)
    min_val <- s1[2] - (1.5 * iqr1)
    max_val <- s1[5] + (1.5 * iqr1)
  }

  if (min_val < s1[1]) {
    min_val <- s1[1]
  }
  if (max_val > s1[6]) {
    max_val <- s1[6]
  }

  iqr1 <- compute_iqr1(s1)
  iqr1 <- max(iqr1, max_diff_rt)

  return(list("min_val" = min_val, "max_val" = max_val, "iqr1" = iqr1))
}

compute_topquant_cor <- function(cor_mz) {
  topquant_cor <- max(0, cor_mz[upper.tri(cor_mz)], na.rm = TRUE)
  return(topquant_cor)
}

subset_selected_mz <- function(check_cor,
                               filter.by,
                               check2,
                               adduct_weights,
                               mchemicaldata,
                               MplusH.abundance.ratio.check,
                               max_diff_rt,
                               cor_mz,
                               corthresh,
                               check_cor2) {
  hub_mz_list <- compute_hub_mz_list(check_cor, mchemicaldata$Adduct, filter.by, check2, adduct_weights)
  hub_mz <- compute_hub_mz(hub_mz_list, mchemicaldata, max_diff_rt)
  diff_rt_hubmz <- compute_diff_rt_hubmz(mchemicaldata, hub_mz)

  if (MplusH.abundance.ratio.check) {
    layer_one_associations <- which(cor_mz[hub_mz, ] >= corthresh & mchemicaldata$mean_int_vec < mchemicaldata$mean_int_vec[hub_mz] & diff_rt_hubmz <= max_diff_rt)
  } else {
    layer_one_associations <- which(cor_mz[hub_mz, ] >= corthresh & diff_rt_hubmz <= max_diff_rt)
  }

  if (length(which(check_cor2 >= 0.1)) > 0) {
    selected_mz <- unique(c(hub_mz, layer_one_associations))
    mchemicaldata <- mchemicaldata[selected_mz, ]
  } else {
    mchemicaldata <- mchemicaldata[which(cor_mz[hub_mz, ] >= corthresh), ]
  }
  return(mchemicaldata)
}

get_data_and_score_for_chemical <- function(cor_mz,
                                            corthresh,
                                            adduct_weights,
                                            filter.by,
                                            mchemicaldata,
                                            max_diff_rt,
                                            chemicalid,
                                            MplusH.abundance.ratio.check,
                                            chemical_score,
                                            mchemicaldata_module) {

  cur_adducts_with_isotopes <- mchemicaldata$Adduct
  check2 <- compute_check2(cur_adducts_with_isotopes)

  check_cor <- sapply(1:dim(cor_mz)[1], function(k) {
    count_mz <- length(which(cor_mz[k, ] >= corthresh)) - 1
    return(count_mz)
  })

  check_cor2 <- check_cor / length(check_cor)

  if (length((cur_adducts_with_isotopes %in% adduct_weights[, 1] == TRUE)) > 0) {
    check_cor[check2 > 0 | (cur_adducts_with_isotopes %in% adduct_weights[, 1] == FALSE)] <- 0
  }

  if (length(which(check_cor > 0) == TRUE) > 0) {
    mchemicaldata <- subset_selected_mz(
      check_cor,
      filter.by,
      check2,
      adduct_weights,
      mchemicaldata,
      MplusH.abundance.ratio.check,
      max_diff_rt,
      cor_mz,
      corthresh,
      check_cor2
    )

    mchemicaldata <- na.omit(mchemicaldata)
    if (nrow(mchemicaldata) < 2) {
      #next
      return(list("score" = chemical_score, "data" = mchemicaldata))
    }

    mchemicaldata <- as.data.frame(mchemicaldata)
    diff_rt <- compute_diff_rt(mchemicaldata$time)
    diff_rt <- round(diff_rt)
    topquant_cor <- compute_topquant_cor(cor_mz)

    if (diff_rt <= max_diff_rt) {
      if (dim(mchemicaldata)[1] > 1) {
        chemical_score <- calc_base_score(mchemicaldata$Adduct, adduct_weights, topquant_cor)
        names(chemical_score) <- chemicalid[1]
      }
    }
    else {
      mchemicaldata$Module_RTclust <- gsub(mchemicaldata$Module_RTclust, pattern = "_[0-9]*", replacement = "")
      mchemicaldata <- cbind(mchemicaldata[, c(2:11)], mchemicaldata[, 1], mchemicaldata[, c(12:14)])
      colnames(mchemicaldata) <- c("mz", "time", "MatchCategory", "theoretical.mz", "chemical_ID", "Name", "Formula", "MonoisotopicMass", "Adduct", "ISgroup", "Module_RTclust", "time.y", "mean_int_vec", "MD")

      mchemicaldata <- as.data.frame(mchemicaldata)
      mchemicaldata$time <- as.numeric(as.character(mchemicaldata$time))

      mchemicaldata <- group_by_rt_histv2(
        mchemicaldata,
        time_step = 1,
        max_diff_rt = max_diff_rt,
        groupnum = unique(mchemicaldata$Module_RTclust)
      )

      top_mod_sub <- table(mchemicaldata$Module_RTclust)
      top_mod_sub_names <- names(top_mod_sub)
      max_top_mod <- which(top_mod_sub == max(top_mod_sub))[1]

      mchemicaldata <- mchemicaldata[which(mchemicaldata$Module_RTclust == top_mod_sub_names[max_top_mod]), ]
      mchemicaldata <- mchemicaldata[order(mchemicaldata$mz), ]

      ranges <- compute_min_max_iqr_advanced(mchemicaldata, max_diff_rt)

      if (nrow(mchemicaldata) < 1) {
        #next
        return(list("score" = chemical_score, "data" = mchemicaldata))
      }

      time_cor_groups <- compute_time_cor_groups(
        mchemicaldata,
        max_diff_rt,
        ranges$min_val,
        ranges$max_val,
        ranges$iqr1
      )

      group_sizes <- sapply(time_cor_groups, function(x) {
        dim(as.data.frame(x))[1]
      })

      group_ind_size <- max(group_sizes)[1]

      if (group_ind_size < 2) {
        mchemicaldata <- mchemicaldata_module[order(mchemicaldata_module$mz), ]
        chemical_score <- calc_base_score_v2(mchemicaldata$Adduct, adduct_weights, topquant_cor)
      }

      good_temp <- {}
      for (g1 in 1:length(group_sizes))
      {
        tempdata <- time_cor_groups[[g1]]
        check_reladd <- which(tempdata$Adduct %in% as.character(adduct_weights[, 1]))
        if (length(check_reladd) > 0) {
          good_temp <- c(good_temp, g1)
          group_ind_size <- min(group_ind_size, nrow(tempdata), na.rm = TRUE)
        }
      }

      if (length(which(group_sizes > 1)) > 0) {
        temp <- compute_temp_score_and_data(
          group_sizes,
          good_temp,
          mchemicaldata,
          time_cor_groups,
          max_diff_rt,
          adduct_weights,
          topquant_cor,
          chemicalid
        )
        mchemicaldata <- temp$data
        chemical_score <- temp$score
      }
      else {
        if (dim(mchemicaldata)[1] > 1) {
          basic <- compute_min_max_iqr_basic(mchemicaldata, max_diff_rt)
          min_val <- basic$min_val
          max_val <- basic$max_val

          if (length(which(mchemicaldata$time >= min_val & mchemicaldata$time <= max_val)) > 1) {
            mchemicaldata <- mchemicaldata[which(mchemicaldata$time >= min_val & mchemicaldata$time <= max_val), ]

            if (dim(mchemicaldata)[1] > 1) {
              mchemicaldata$time <- as.numeric(as.character(mchemicaldata$time))

              chemical_score <- 0
              if (length(mchemicaldata) > 0 & nrow(mchemicaldata) > 1) {
                chemical_score <- calc_base_score(mchemicaldata$Adduct, adduct_weights, topquant_cor)
              }
            }
          } else {
            chemical_score <- 0
            if (length(mchemicaldata) > 0 & nrow(mchemicaldata) > 1) {
              diff_rt <- compute_diff_rt(mchemicaldata$time)
              if (diff_rt > max_diff_rt) {
                chemical_score <- calc_base_score_v2(mchemicaldata$Adduct, adduct_weights, topquant_cor)
              }
            }
          }
        }
        names(chemical_score) <- chemicalid[1]
      }
      names(chemical_score) <- chemicalid[1]
      mchemicaldata <- na.omit(mchemicaldata)
    }
  }
  else {
    # no correlation between putative adducts
    good_adducts_len <- length(which(mchemicaldata$Adduct %in% adduct_weights[, 1]))
    if (good_adducts_len > 0) {
      max_adduct_weight <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% mchemicaldata$Adduct), 2])))[1]
      chemical_score <- 10^max_adduct_weight

      good_adduct_index <- which(adduct_weights[, 2] == max_adduct_weight)
      mchemicaldata <- mchemicaldata[which(mchemicaldata$Adduct %in% adduct_weights[good_adduct_index, 1]), ]
    } else {
      chemical_score <- 0
      mchemicaldata <- mchemicaldata_module
      mchemicaldata <- mchemicaldata[order(mchemicaldata$mz), ]
    }
  }
  return(list("score" = chemical_score, "data" = mchemicaldata))
}

compute_cor_mz <- function(mzid_cur, global_cor) {
  cor_mz <- round(global_cor[mzid_cur, mzid_cur], 1)

  if (length(cor_mz) > 1) {
    corrownamesA <- rownames(cor_mz)
    mat_rownames <- strsplit(as.character(corrownamesA), split = "_")

    m1 <- {}
    for (i in 1:length(mat_rownames)) {
      m1 <- rbind(m1, cbind(mat_rownames[[i]][1], mat_rownames[[i]][2]))
    }

    m1 <- as.data.frame(m1)
    colnames(m1) <- c("mz", "time")
    cor_mz2 <- cbind(m1, cor_mz)

    mz_order <- order(cor_mz2$mz)
    cor_mz2 <- cor_mz2[c(mz_order), ]

    cor_mz <- cor_mz2[, -c(1:2)]
    cor_mz <- cor_mz[, mz_order]
  }
  return(cor_mz)
}

get_conf_level <- function(mchemicaldata, adduct_weights) {
  conf_level <- 0
  if (length(mchemicaldata) > 0 && nrow(mchemicaldata) > 0) {
    conf_level <- get_confidence_stage2(curdata = mchemicaldata, adduct_weights = adduct_weights)
    conf_level <- as.numeric(as.character(conf_level))
  }
  return(conf_level)
}

get_min_chemscore <- function(corthresh, max_diff_rt) {
  return(100 * 2 * (1 * (corthresh)) * (1 / ((max_diff_rt * 0.1) + 1)^3))
}

apply_rt_scaling <- function(mchemicaldata, max_diff_rt, chemical_score) {
  if (length(mchemicaldata) > 0) {
    if (nrow(mchemicaldata) > 1) {
      diff_rt <- max(mchemicaldata$time) - min(mchemicaldata$time)

      k_power <- 1
      if (diff_rt > max_diff_rt) {
        k_power <- 10
      }
      chemical_score <- chemical_score * (1 / ((diff_rt * 0.1) + 1)^k_power)
    }
  } else {
    chemical_score <- 0
  }
  return(chemical_score)
}

compute_check2 <- function(cur_adducts_with_isotopes) {
  check2 <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])")
  return(check2)
}

update_best_score_and_data <- function(chemical_score, best_chemical_score, mchemicaldata, best_data, conf_level = NA) {

  condition <- chemical_score > best_chemical_score
  if(!is.na(conf_level)) {
    condition <- condition & conf_level > 0
  }

  if (condition) {
    best_chemical_score <- chemical_score
    best_data <- mchemicaldata
  } else {
    if (chemical_score == best_chemical_score) {
      best_chemical_score <- chemical_score
      best_data <- rbind(best_data, mchemicaldata)
    }
  }
  return(list("score" = best_chemical_score, "data" = best_data))
}

compute_best_score <- function(table_mod,
                               mchemicaldata_orig,
                               adduct_weights,
                               global_cor,
                               corthresh,
                               max_diff_rt,
                               filter.by,
                               chemicalid,
                               MplusH.abundance.ratio.check) {
  top_mod <- names(table_mod)

  best_chemical_score <- (-100)
  best_data <- mchemicaldata_orig

  for (i in 1:length(which(table_mod >= 1))) {
    chemical_score <- (-99999)
    conf_level <- 0
    mchemicaldata <- mchemicaldata_orig %>%
      filter(Module_RTclust == top_mod[i]) %>%
      arrange(mz)

    if (nrow(mchemicaldata) < 2) {
      good_adducts_len <- length(which(mchemicaldata$Adduct %in% adduct_weights[, 1]))
      if (good_adducts_len > 0) {
        chemical_score <- compute_score(adduct_weights, mchemicaldata$Adduct)
        best <- update_best_score_and_data(chemical_score, best_chemical_score, mchemicaldata, best_data)
        best_chemical_score <- best$score
        best_data <- best$data
      }
      # next;
    }

    mzid_cur <- paste(mchemicaldata$mz, mchemicaldata$time, sep = "_")
    cor_mz <- compute_cor_mz(mzid_cur, global_cor)

    if (length(cor_mz) > 1) {
      intermediate_result <- get_data_and_score_for_chemical(
        cor_mz,
        corthresh,
        adduct_weights,
        filter.by,
        mchemicaldata,
        max_diff_rt,
        chemicalid,
        MplusH.abundance.ratio.check,
        chemical_score,
        mchemicaldata_orig[which(mchemicaldata_orig$Module_RTclust == top_mod[i]), ]
      )

      mchemicaldata <- intermediate_result$data
      chemical_score <- intermediate_result$score
    }

    check2 <- compute_check2(mchemicaldata$Adduct)

    if (length(check2) > 0) {
      for (a1 in 1:length(check2)) {
        strlength <- attr(check2[[a1]], "match.length")
        if (strlength[1] > (-1)) {
          cur_adducts <- gsub(mchemicaldata$Adduct, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
          count_abundant_form <- length(which(cur_adducts %in% cur_adducts[a1]))
          if (count_abundant_form < 2) {
            mchemicaldata <- mchemicaldata[-a1, ]
          }
        }
      }
    }

    conf_level <- get_conf_level(mchemicaldata, adduct_weights)
    chemical_score <- apply_rt_scaling(mchemicaldata, max_diff_rt, chemical_score)
    min_chemical_score <- get_min_chemscore(corthresh, max_diff_rt)
    chemical_score <- if (chemical_score > min_chemical_score) chemical_score * (conf_level^conf_level) else 0

    best <- update_best_score_and_data(chemical_score, best_chemical_score, mchemicaldata, best_data, conf_level)
    best_chemical_score <- best$score
    best_data <- best$data
  }
  return(list("score" = best_chemical_score, "data" = best_data))
}

compute_table_mod <- function(data) {
  table_mod <- table(data)
  table_mod <- table_mod[table_mod > 0]
  table_mod <- table_mod[order(table_mod, decreasing = TRUE)]
  return(table_mod)
}

#' @export
compute_chemical_score <- function(mchemicaldata,
                                   adduct_weights,
                                   global_cor,
                                   corthresh,
                                   filter.by,
                                   max_diff_rt,
                                   chemicalid,
                                   MplusH.abundance.ratio.check) {
  table_mod <- compute_table_mod(mchemicaldata$Module_RTclust)

  mchemicaldata_orig <- mchemicaldata

  if (length(which(table_mod >= 1)) > 0) {
    best <- compute_best_score(
      table_mod,
      mchemicaldata_orig,
      adduct_weights,
      global_cor,
      corthresh,
      max_diff_rt,
      filter.by,
      chemicalid,
      MplusH.abundance.ratio.check
    )

    ## for loop complete

    if (best$score > 0) {
      chemical_score <- best$score
      mchemicaldata <- best$data
      names(chemical_score) <- chemicalid[1]
    } else {
      chemical_score <- 0
    }
  }
  ####### add code for only correlation criteria here
  if (chemical_score <= 1) {
    mchemicaldata <- mchemicaldata_orig
    cur_adducts_with_isotopes <- mchemicaldata$Adduct
    good_adducts_len <- length(which(cur_adducts_with_isotopes %in% adduct_weights[, 1]))

    if (good_adducts_len > 0) {
      max_adduct_weight <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% cur_adducts_with_isotopes), 2])))[1]
      chemical_score <- ((10^max_adduct_weight))
      chemical_score <- chemical_score[1] - 1
      chemical_score <- chemical_score[1]

      good_adduct_index <- which(adduct_weights[, 2] == max_adduct_weight)
      mchemicaldata <- mchemicaldata[which(cur_adducts_with_isotopes %in% adduct_weights[good_adduct_index, 1]), ]
    }
  }

  if (nrow(mchemicaldata) > 0) {
    mchemicaldata <- unique(mchemicaldata)
  } else {
    chemical_score <- 0
  }

  result <- list("chemical_score" = chemical_score, "filtdata" = mchemicaldata)
  return(result)
}


#' @import tidyr
#' @import dplyr
#' @import plyr
#' @importFrom magrittr %>%
#' @export
get_chemscorev1.6.71 <- function(chemicalid,
                                 mchemicaldata,
                                 corthresh,
                                 global_cor,
                                 mzid, max_diff_rt = 10,
                                 level_module_isop_annot,
                                 adduct_table,
                                 adduct_weights,
                                 filter.by = c("M+H"),
                                 max_isp = 100,
                                 MplusH.abundance.ratio.check = TRUE,
                                 mass_defect_window = 0.01,
                                 mass_defect_mode = "pos",
                                 outlocorig) {
  setwd(outlocorig)

  outloc1 <- paste(outlocorig, "/stage2/", sep = "")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1)

  if (length(mchemicaldata$mz) < 1) stop("No mz data found!")

  mchemicaldata <- add_isotopic_peaks(
    mchemicaldata,
    adduct_weights,
    exp_isp,
    level_module_isop_annot,
    max_diff_rt,
    mass_defect_window,
    mass_defect_mode,
    max_isp,
    abund_ratio_vec
  )

  result <- compute_chemical_score(
    mchemicaldata,
    adduct_weights,
    global_cor,
    corthresh,
    filter.by,
    max_diff_rt,
    chemicalid,
    MplusH.abundance.ratio.check
  )

  #rm("mzid", "global_cor", "temp_global_cor")
  setwd("..")

  return(result)
}