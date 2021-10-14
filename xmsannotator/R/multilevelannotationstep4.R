compute_confidence_levels <- function(c,
                                      chemids,
                                      chemscoremat,
                                      filter.by,
                                      max.rt.diff,
                                      adduct_weights,
                                      max_isp,
                                      min_ions_perchem) {
    cur_chemid <- chemids[c]

    curdata <- chemscoremat[which(chemscoremat$chemical_ID == cur_chemid), ]
    curdata <- curdata[order(curdata$Adduct), ]
    
    bool_check <- 1

    if (!is.na(filter.by)) {
        check_adduct <- which(curdata$Adduct %in% filter.by)
        if (length(check_adduct) <= 0) {
            bool_check <- 0
        }
    }

    Confidence <- 0
    if (bool_check == 1) {
        curdata <- get_confidence_stage4(
                        curdata,
                        max.rt.diff,
                        adduct_weights = adduct_weights,
                        filter.by = filter.by,
                        max_isp = max_isp,
                        min_ions_perchem = min_ions_perchem
                    )
        if (!is.na(curdata[1, 1])) {
            Confidence <- as.numeric(as.character(curdata[, 1]))
            if (Confidence < 2) {
                if (length(which(curdata$Adduct %in% adduct_weights[which(as.numeric(adduct_weights[, 2]) > 0), 1])) > 0) {
                    if (curdata$score > 10) {
                        mnum <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% curdata$Adduct), 2])))[1]
                        curdata <- curdata[which(curdata$Adduct %in% adduct_weights[which(as.numeric(as.character(adduct_weights[, 2])) >= mnum), 1]), ]
                        Confidence <- 2
                    }
                }
            }
        }
    } else {
        if (length(which(curdata$Adduct %in% adduct_weights[, 1])) > 0) {
            if (curdata$score >= 10) {
                mnum <- max(as.numeric(as.character(adduct_weights[which(adduct_weights[, 1] %in% curdata$Adduct), 2])))[1]
                if (length(which(curdata$Adduct %in% filter.by)) > 0) {
                    curdata <- curdata[which(curdata$Adduct %in% filter.by), ]
                    Confidence <- 2
                }
            }
        }
    }

    if (nrow(curdata) > 1) {
        if (curdata$score < 10) {
            if (length(unique(curdata$Adduct)) < 2) {
                Confidence <- 0
            }
        }
    }

    curdata <- cbind(Confidence, curdata)
    curdata <- as.data.frame(curdata)
    curdata <- curdata[, c("Confidence", "chemical_ID")]
    curdata <- unique(curdata)
    return(curdata)
}

compute_delta_ppm <- function(chemscoremat_with_confidence) {
    # this is fishy but necessary 
    chemscoremat_with_confidence$mz <- as.numeric(as.character(chemscoremat_with_confidence$mz))
    chemscoremat_with_confidence$theoretical.mz <- as.numeric(as.character(chemscoremat_with_confidence$theoretical.mz))
    
    chemscoremat_with_confidence_temp <- chemscoremat_with_confidence[, c("mz", "theoretical.mz")]
    chemscoremat_with_confidence_temp <- apply(chemscoremat_with_confidence_temp, 1, as.numeric)
    chemscoremat_with_confidence_temp <- t(chemscoremat_with_confidence_temp)
    chemscoremat_with_confidence_temp <- as.data.frame(chemscoremat_with_confidence_temp)
    
    delta_ppm <- apply(chemscoremat_with_confidence_temp, 1, function(x) {
        return(10^6 * abs(x[2] - x[1]) / (x[2]))
    })
    delta_ppm <- round(delta_ppm, 2)
    
    chemscoremat_with_confidence <- cbind(chemscoremat_with_confidence[, 1:8], delta_ppm, chemscoremat_with_confidence[, 9:dim(chemscoremat_with_confidence)[2]])
    chemscoremat_with_confidence <- chemscoremat_with_confidence[order(chemscoremat_with_confidence$Confidence, decreasing = TRUE), ]
    return(chemscoremat_with_confidence)
}

boost_confidence_of_IDs <- function(chemscoremat_with_confidence, boostIDs, max.mz.diff, max.rt.diff) {
    cnames_boost <- colnames(boostIDs)
    
    if (length(cnames_boost) > 1) {
        chemscoremat_with_confidence_mzrt <- chemscoremat_with_confidence[, c("mz", "time")]
        validated_mzrt <- boostIDs[, c("mz", "time")]
        
        ghilicpos <- getVenn(chemscoremat_with_confidence_mzrt,
                             name_a = "exp", validated_mzrt, name_b = "boost", mz.thresh = max.mz.diff, time.thresh = max.rt.diff,
                             alignment.tool = NA, xMSanalyzer.outloc = getwd(), use.unique.mz = FALSE, plotvenn = FALSE
        )
        
        save(ghilicpos, file = "ghilicpos.Rda")
        
        g1 <- ghilicpos$common
        rm(ghilicpos)
        
        if (!is.na(max.rt.diff)) {
            t1 <- table(g1$index.B)
            ind_names <- names(t1)
            parent_bad_ind <- {}
        }
        
        t1 <- table(chemscoremat_with_confidence$Confidence, chemscoremat_with_confidence$chemical_ID)
        cnames <- colnames(t1)
        cnames <- cnames[which(cnames %in% boostIDs$ID)]
        
        good_ind_1 <- {}
        
        for (ind2 in 1:dim(g1)[1]) {
            temp_ind1 <- g1$index.A[ind2]
            temp_ind2 <- g1$index.B[ind2]
            
            
            if (chemscoremat_with_confidence$chemical_ID[temp_ind1] %in% boostIDs$ID[temp_ind2]) {
                good_ind_1 <- c(good_ind_1, g1$index.A[ind2])
            }
        }
        
        overlap_mz_time_id <- good_ind_1
        
        chemscoremat_with_confidence$Confidence[overlap_mz_time_id] <- 4
        chemscoremat_with_confidence$score[overlap_mz_time_id] <- chemscoremat_with_confidence$score[overlap_mz_time_id] * 100
        t1 <- table(chemscoremat_with_confidence$Confidence[overlap_mz_time_id], chemscoremat_with_confidence$chemical_ID[overlap_mz_time_id])
        
        cnames1 <- colnames(t1)
        cnames2 <- cnames1[which(t1 > 0)]
    } else {
        good_ind <- which(chemscoremat_with_confidence$chemical_ID %in% boostIDs)
        if (length(good_ind) > 0) {
            chemscoremat_with_confidence$Confidence[good_ind] <- 4
            chemscoremat_with_confidence$score[good_ind] <- chemscoremat_with_confidence$score[good_ind] * 100
        }
    }
    return(chemscoremat_with_confidence)
}

#' @importFrom foreach foreach %do% %dopar%
multilevelannotationstep4 <- function(outloc,
                                      chemscoremat,
                                      max.mz.diff = 5,
                                      max.rt.diff = 30,
                                      adduct_weights = NA,
                                      filter.by = NA,
                                      min_ions_perchem = 1,
                                      boostIDs = NA,
                                      max_isp = 5,
                                      dbAllinf = NA) {
    setwd(outloc)

    chemids <- unique(chemscoremat$chemical_ID)

    data(adduct_table)
    adduct_table <- adduct_table[order(adduct_table$Adduct), ]
    
    if (is.na(adduct_weights)) {
        adduct_weights <- data.frame(Adduct = c("M+H", "M-H"), Weight = c(1, 1))
    }

    # assign confidence level
    chemscoremat_conf_levels <- lapply(
        seq_len(length(chemids)),
        compute_confidence_levels,
        chemids,
        chemscoremat,
        filter.by,
        max.rt.diff,
        adduct_weights,
        max_isp,
        min_ions_perchem
    )
    
    chemscoremat_conf_levels <- plyr::ldply(chemscoremat_conf_levels, rbind)
    chemscoremat_conf_levels <- as.data.frame(chemscoremat_conf_levels)

    chemscoremat_with_confidence <- merge(chemscoremat_conf_levels, unique(chemscoremat), by = "chemical_ID")
    chemscoremat_with_confidence <- as.data.frame(chemscoremat_with_confidence)
    
    cnames1 <- colnames(chemscoremat_with_confidence)
    
    chemscoremat_with_confidence <- compute_delta_ppm(chemscoremat_with_confidence)

    # (II) Presence of required adducts/forms specified by the user 
    # for assignment to high confidence categories (e.g., M + H).
    if (!is.na(boostIDs)) {
        chemscoremat_with_confidence <- boost_confidence_of_IDs(chemscoremat_with_confidence, boostIDs, max.mz.diff, max.rt.diff)
    }
    
    t2 <- table(chemscoremat_with_confidence$mz)
    uniquemz <- names(which(t2 == 1))

    # assign match category
    chemscoremat_with_confidence$MatchCategory <- rep("Multiple", dim(chemscoremat_with_confidence)[1])
    chemscoremat_with_confidence$MatchCategory[which(chemscoremat_with_confidence$mz %in% uniquemz)] <- "Unique"

    write.csv(chemscoremat_with_confidence, file = "Stage4.csv", row.names = FALSE)

    chemscoremat_with_confidence <- as.data.frame(chemscoremat_with_confidence)
    chemscoremat_with_confidence <- chemscoremat_with_confidence[order(chemscoremat_with_confidence$Confidence, decreasing = TRUE), ]

    print("Stage 4 confidence level distribution for unique chemical/metabolite IDs")
    print(table(chemscoremat_with_confidence$Confidence[-which(duplicated(chemscoremat_with_confidence$chemical_ID) == TRUE)]))

    print("Stage 4 confidence level distribution for unique chemical/metabolite formulas")
    print(table(chemscoremat_with_confidence$Confidence[-which(duplicated(chemscoremat_with_confidence$Formula) == TRUE)]))

    return(chemscoremat_with_confidence)
}