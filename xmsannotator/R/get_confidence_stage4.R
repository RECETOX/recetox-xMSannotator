filter_clusters <- function(curdata, table_names) {
  curdata <- curdata[which(curdata$Module_RTclust == table_names[1]), ]
  return(curdata)
}

create_cluster_table <- function(curdata) {
  cluster_table <- table(curdata$Module_RTclust)
  cluster_table <- cluster_table[which(cluster_table > 0)]
  cluster_table <- cluster_table[which(cluster_table == max(cluster_table)[1])]
  return(cluster_table)
}

compute_delta_rt <- function(curdata) {
  delta_rt <- max(curdata$time) - min(curdata$time)
  delta_rt <- round(delta_rt)
  return(delta_rt)
}

assign_conf <- function(curdata, filter.by, delta_rt, max_diff_rt, conf_level) {
  if (curdata$score[1] > 0) {
    if (length(which(curdata$Adduct %in% filter.by)) > 0) {
      if (nrow(curdata) > 1 && length(unique(curdata$Adduct)) > 1 && delta_rt < max_diff_rt) {
        conf_level <- "High"
      } else {
        conf_level <- "Medium"
      }
    } else {
      conf_level <- "Low"
    }
  }
  return(conf_level)
}

get_confidence_stage4 <-function(curdata,
                                 max_diff_rt,
                                 adduct_weights = NA,
                                 filter.by = NA,
                                 min_ions_perchem = 1,
                                 max_isp = 5) {
    curdata <- curdata[order(curdata$Adduct), ]
    
    # (I) Score is greater than zero.
    if (curdata$score[1] < 0.1) {
      score_level <- 0
      final_res <- cbind(score_level, curdata)
      final_res <- as.data.frame(final_res)
      return(final_res)
    }
    
    curdata_all <- curdata
    cur_adducts_with_isotopes <- curdata$Adduct
    
    cur_adducts <-gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
    
    if (is.na(adduct_weights)) {
      adduct_weights <- data.frame(Adduct = c("M+H", "M-H"), Weight = c(1, 1))
    }
  
    adduct_monoisot <- data.frame(cbind("M", 1, 1, 0, "-", "S"))
    colnames(adduct_monoisot) <- colnames(adduct_table)
    adduct_table <- rbind(adduct_table, adduct_monoisot)
    adduct_table <- adduct_table[order(adduct_table$Adduct), ]
    
    curdata <- cbind(curdata, cur_adducts)
    curdata <- merge(curdata, adduct_table, by.x = "cur_adducts", by.y = "Adduct")
    curdata <- curdata[, -c(1)]
    
    chemscoremat_conf_levels <- "High"
    if (length(which(adduct_table$Adduct %in% cur_adducts)) < 1) {
      chemscoremat_conf_levels <- "None"
    }
    
    if (length(unique(curdata$Adduct)) < 2) {
      if (curdata$score < 10 && length(which(cur_adducts %in% filter.by)) < 1) {
        chemscoremat_conf_levels <- "None"
      } else if (length(which(cur_adducts %in% filter.by)) > 0) {
        chemscoremat_conf_levels <- "Low"
      }
     
      if (nrow(curdata) > 1) {
        if (nrow(curdata) == 2) {
          chemscoremat_conf_levels <- "Medium"
        }
        
        if (length(which(cur_adducts %in% filter.by)) < 1) {
          chemscoremat_conf_levels <- "None"
        } else {
          curdata <- curdata[which(cur_adducts %in% filter.by), ]
          final_res <- cbind(2, curdata)
          return(final_res)
        }
      } 
    }
    
    curdata$time <- as.numeric(as.character(curdata$time))
    delta_rt <- compute_delta_rt(curdata)
    
    if (delta_rt > max_diff_rt) {
      chemscoremat_conf_levels <- "None"
      
      groupB <- group_by_rt(curdata, time_step = 3, max.rt.diff = max_diff_rt, groupnum = "")
      groupB <- groupB[order(groupB$mz), ]
      
      cname1 <- paste(groupB$mz, groupB$time, sep = "_")
      if (length(which(duplicated(cname1) == TRUE)) > 0) {
        groupB <- groupB[-which(duplicated(cname1) == TRUE), ]
      }
      curdata <- curdata[order(curdata$mz), ]

      module_clust <- groupB[, 1]
      curdata$Module_RTclust <- module_clust
      
      if (length(which(adduct_weights[, 1] %in% cur_adducts)) > 0 && curdata$score[1] > 0.1) {
        if (is.na(filter.by)) {
          good_mod <- curdata$Module_RTclust[which(curdata$Adduct %in% adduct_weights[, 1])]
        } else {
          good_mod <- curdata$Module_RTclust[which(curdata$Adduct %in% filter.by)]
        }
        
        curdata <- curdata[which(curdata$Module_RTclust %in% good_mod), ]
        cluster_table <- create_cluster_table(curdata)
        curdata <- filter_clusters(curdata, names(cluster_table))
        delta_rt <- compute_delta_rt(curdata)
        
        if (is.na(filter.by)) {
          if (curdata$score[1] > 0 && nrow(curdata) > 1 && length(unique(curdata$Adduct)) > 1 && delta_rt < max_diff_rt) {
            chemscoremat_conf_levels <- "High"
          } else {
            chemscoremat_conf_levels <- "Low"
          }
        } else if (length(which(cluster_table > 1)) < 1) {
            chemscoremat_conf_levels <- assign_conf(curdata, filter.by, delta_rt, max_diff_rt, chemscoremat_conf_levels)
        }
      }
      module_clust <- paste(curdata$module_num, curdata$Module_RTclust, sep = "")
      curdata$Module_RTclust <- module_clust
    }
    curdata <- as.data.frame(curdata)
    
    temp_curdata <- curdata
    temp_curdata$Module_RTclust <- gsub(temp_curdata$Module_RTclust, pattern = "(_[0-9]*)", replacement = "")
    table_modules <- table(temp_curdata$Module_RTclust)
    rm(temp_curdata)
    module_names <- names(table_modules[which(table_modules > 0)])
    
    if (length(module_names) > 1) {
      if (curdata$score < 10) {
        chemscoremat_conf_levels <- "None"
      } else if (curdata$score > 10 && length(which(cur_adducts %in% filter.by)) > 0) {
        chemscoremat_conf_levels <- "Medium"
      }
    }
    
    if (nrow(curdata) < 2) {
      if (length(which(cur_adducts %in% adduct_weights[, 1])) < 1) {
        chemscoremat_conf_levels <- "Low"
      } else {
        # matches an M+H and score is greater than 10
        if (curdata$score > 10 && length(which(cur_adducts %in% filter.by)) > 0) {
          chemscoremat_conf_levels <- "Medium"
        } else {
          chemscoremat_conf_levels <- "None"
        }
      }
    } else {
      min_molecules <- min(curdata$num_molecules)
      # number of molecules check
      if (min_molecules > 1) {
        chemscoremat_conf_levels <- "Low"
      } else {
        # (V) Abundance ratio checks for isotopes, multimers (dimers and trimers), 
        # and multiply charged adducts with respect to the singly charged adducts 
        # and ions according to heuristic rules.
        
        check_abundance <- gregexpr(text = cur_adducts, pattern = "([2-3]+M)")
        if (length(check_abundance) > 0) {
          min_mol_ind <- which(curdata$num_molecules == min_molecules)
          max_int_min_mol <- max(curdata$mean_int_vec[min_mol_ind])
          bad_ind_status <- rep(0, length(check_abundance))
          min_molecules <- min(adduct_table[which(adduct_table$Adduct %in% cur_adducts), 2])
          
          if (min_molecules < 2) {
            # check pattern of each adduct
            for (abundant in 1:length(check_abundance)) {
              strlength <- attr(check_abundance[[abundant]], "match.length")
              if (strlength[1] > (-1)) {
                abund_ratio <- curdata$mean_int_vec[abundant] / max_int_min_mol
                
                if (!is.na(abund_ratio)) {
                  if (abund_ratio > 1) {
                    bad_ind_status[abundant] <- 1
                    if (chemscoremat_conf_levels == "High") {
                      chemscoremat_conf_levels <- "Medium"
                    } else {
                      chemscoremat_conf_levels <- "Low"
                    }
                  }
                } else {
                  chemscoremat_conf_levels <- "Low"
                  bad_ind_status[abundant] <- 1
                }
              }
            }
            
            if (length(which(bad_ind_status == 1)) > 0) {
              bad_adducts <- cur_adducts[which(bad_ind_status == 1)]
            }
            
            if (length(nrow(curdata) > 0)) {
              cur_adducts_with_isotopes <- curdata$Adduct
              
              cur_adducts <-gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
            } else {
              chemscoremat_conf_levels <- "None"
            }
          } else {
            chemscoremat_conf_levels <- "Low"
          }
          adduct_charges <- curdata$charge
          min_charge <- min(curdata$charge)
          min_charge_ind <- which(curdata$charge == min_charge)
          
          high_charge_ind <- which(curdata$charge > 1)
          max_int_min_charge <- max(curdata$mean_int_vec)
          
          if (length(high_charge_ind) > 0) {
            abund_ratio_min_maxcharge <-
              max(curdata$mean_int_vec[min_charge_ind])[1] / max(curdata$mean_int_vec[high_charge_ind])[1]
            
            if (abund_ratio_min_maxcharge < 1) {
              chemscoremat_conf_levels <- "Low"
            }
          }
          
          bad_ind_status <- rep(1, length(check_abundance))
          
          if (min_charge > 1) {
            chemscoremat_conf_levels <- "Low"
          }
          
          # isotope based check for +1 and +2 to make sure that the abundant form is present
          check_abundance <- gregexpr(text = cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])")
          
          if (length(check_abundance) > 0) {
            for (abundant in 1:length(check_abundance)) {
              strlength <- attr(check_abundance[[abundant]], "match.length")
              
              if (strlength[1] > (-1)) {
                count_abundant_form <- length(which(cur_adducts %in% cur_adducts[abundant]))
                
                if (count_abundant_form < 2) {
                  curdata <- curdata[-c(abundant), ]
                }
              }
            }
          }
          
          cur_adducts_with_isotopes <- curdata$Adduct

          cur_adducts <- gsub(cur_adducts_with_isotopes, pattern = "(_\\[(\\+|\\-)[0-9]*\\])", replacement = "")
          formula_vec <- curdata$Formula
          curformula <- as.character(formula_vec[1])
        }
      }
    }
  
    # (IV) Hydrogen/carbon ratio check.
    formula_vec <- curdata$Formula
    curformula <- as.character(formula_vec[1])
    curformula <- gsub(curformula, pattern = "Ca|Cl|Cd|Cr", replacement = "")
    
    numcarbons <- check_element(curformula, "C")
    if (numcarbons < 1) {
      chemscoremat_conf_levels <- "None"
    }
    
    if (length(unique(curdata$Adduct)) < 2) {
      if (nrow(curdata) > 1) {
        if (nrow(curdata) == 2) {
          chemscoremat_conf_levels <- "Medium"
        }
        
        if (length(which(cur_adducts %in% filter.by)) < 1) {
          chemscoremat_conf_levels <- "None"
        }
      } else {
        if (curdata$score < 10) {
          if (length(which(cur_adducts %in% filter.by)) < 1) {
            chemscoremat_conf_levels <- "None"
          } else {
            chemscoremat_conf_levels <- "Low"
          }
        } else if (curdata$score > 10 && length(which(cur_adducts %in% filter.by)) > 0) {
          chemscoremat_conf_levels <- "Medium"
        }
      }
    }
    
    if (nrow(curdata) < 1) {
      score_level <- 0
      curdata <- curdata_all
    } else {
      # 3 -> High; 2: -> Medium; 1 -> Low; 0 -> None
      name_to_score = data.frame(row.names = c("High", "Medium", "Low", "None"), 
                                 val=c(3, 2, 1, 0))
      
      score_level <- name_to_score[chemscoremat_conf_levels,]
    }
    
    if (is.na(score_level[1])) {
      stop(curdata)
    }
    
    final_res <- cbind(score_level, curdata)
    final_res <- as.data.frame(final_res)
    
    if (length(which(is.na(final_res$score_level))) > 0) {
      final_res$score_level <-
        replace(final_res$score_level, which(is.na(final_res$score_level)), 0)
    }
    
    num_uniq_adducts <- length(unique(curdata$Adduct))
    uniq_adducts <- unique(curdata$Adduct)
    num_good_adducts <- length(which(uniq_adducts %in% filter.by))
    
    if (num_good_adducts > 0) {
      final_res$score <- final_res$score * num_good_adducts
    } else {
      final_res$score <- 0
    }
    
    if (nrow(final_res) < min_ions_perchem) {
      final_res$score_level <- 0
    }

    final_res <- as.data.frame(final_res)
    
    rm(
      list = c(
        "curdata_all",
        "cur_adducts_with_isotopes",
        "cur_adducts",
        "adduct_monoisot",
        "curdata",
        "temp_curdata",
        "table_modules",
        "chemscoremat_conf_levels",
        "score_level",
        "groupB",
        "good_mod",
        "module_clust",
        "uniq_adducts",
        "num_good_adducts",
        "num_uniq_adducts"
      )
    )
    
    return(final_res)
  }
