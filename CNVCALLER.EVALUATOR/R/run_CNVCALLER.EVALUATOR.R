run_CNVCALLER.EVALUATOR <- function(calls,
                                    refs,
                                    targets,
                                    seg_dups,
                                    parameters){
  TP <- 0
  TP_for_FP <- 0
  FP <- 0
  TN <- 0
  FN <- 0
  calls[,"chr"] <- as.character(calls[,"chr"])
  refs[,"chr"] <- as.character(refs[,"chr"])
  targets[,"chr"] <- as.character(targets[,"chr"])
  seg_dups[,"chr"] <- as.character(seg_dups[,"chr"])
  min_frequency <- as.double(parameters$min_frequency)
  max_frequency <- as.double(parameters$max_frequency)
  min_length <- strtoi(parameters$min_length)
  max_length <- strtoi(parameters$max_length)
  min_targets <- strtoi(parameters$min_targets)
  max_targets <- strtoi(parameters$max_targets)
  seg_dups_filter_enable <- strtoi(parameters$seg_dups_filter_enable)
  cnv_type <- parameters$cnv_type
  samples <- unique(refs[,"sample_name"])
  num_of_original_samples_in_refs <- length(samples)
  chromosomes <- c(1:22, "X", "Y", paste0("chr",c(1:22, "X", "Y")))
  for(chromosome in chromosomes) {
    print(paste("Processing chr: ", chromosome, sep=""))
    calls_for_chr <- subset(calls, chr == chromosome) # calls
    refs_for_chr <- subset(refs, chr == chromosome) # refs
    if (nrow(calls_for_chr) == 0 && nrow(refs_for_chr) == 0) {
      next()
    }
    targets_for_chr <- subset(targets, chr == chromosome)
    seg_dups_for_chr <- subset(seg_dups, chr == chromosome)
    ## transfomating seg_dups to targets overlapped with seg_dups
    indexes_to_delete <- c()
    for (i in 1:nrow(targets_for_chr)) {
      num_of_seg_dups <- nrow(subset(seg_dups_for_chr, chr == targets_for_chr[i,'chr'] & !((targets_for_chr[i,'ed_bp'] < st_bp) | (ed_bp < targets_for_chr[i,'st_bp']))))
      if (num_of_seg_dups > 0) {
        indexes_to_delete <- c(indexes_to_delete, i)
      }
    }
    if (length(indexes_to_delete) != 0) {
      seg_dups_for_chr <- targets_for_chr[-indexes_to_delete,]
    }

    refs_for_chr_for_TP <- refs_for_chr
    refs_for_chr_for_TP <- filter_cnvs_by_frequency(refs_for_chr_for_TP, min_frequency, max_frequency)
    refs_for_chr_for_TP <- filter_cnvs_by_targets(refs_for_chr_for_TP, targets_for_chr, min_targets, max_targets)
    calls_for_chr_for_TP <- calls_for_chr

    # TODO filter by lratio from inst/ file
    #refs_for_chr <- filter_cnvs_by_cnv_type(refs_for_chr, cnv_type)
    #calls_for_chr <- filter_cnvs_by_cnv_type(calls_for_chr, cnv_type)
    # TODO delete freq from ref and after that also from calls (cnvs for the same positions)
    #refs_for_chr <- filter_cnvs_by_frequency(refs_for_chr, min_frequency, max_frequency)
    calls_for_chr <- filter_cnvs_by_frequency(calls_for_chr, min_frequency, max_frequency)
    # TODO delete length from ref and after that also from calls (cnvs for the same positions)
    #refs_for_chr <- filter_cnvs_by_length(refs_for_chr, min_length, max_length)
    #calls_for_chr <- filter_cnvs_by_length(calls_for_chr, min_length, max_length)
    # TODO delete targets from ref and after that also from calls (cnvs for the same positions)
    #refs_for_chr <- filter_cnvs_by_targets(refs_for_chr, targets_for_chr, min_targets, max_targets)
    calls_for_chr <- filter_cnvs_by_targets(calls_for_chr, targets_for_chr, min_targets, max_targets)
    # TODO delete seg_dups from ref and after that also from calls (cnvs for the same positions)
    #refs_for_chr <- filter_cnvs_by_seg_dups(refs_for_chr, seg_dups_for_chr, seg_dups_filter_enable)
    calls_for_chr <- filter_cnvs_by_seg_dups(calls_for_chr, seg_dups_for_chr, seg_dups_filter_enable)
    if (nrow(calls_for_chr) == 0 && nrow(refs_for_chr) == 0) {  # TODO
      next()
    }
    for(sample in samples) {
      calls_for_chr_for_sample <- subset(calls_for_chr, sample_name == sample)
      refs_for_chr_for_sample <- subset(refs_for_chr, sample_name == sample)
      calls_for_chr_for_TP_for_sample <- subset(calls_for_chr_for_TP, sample_name == sample)
      refs_for_chr_for_TP_for_sample <- subset(refs_for_chr_for_TP, sample_name == sample)
      if (nrow(calls_for_chr_for_sample) == 0 && nrow(refs_for_chr_for_sample) == 0) {  # TODO
        # next()
      } else {
        intersection_matrix <- build_intersection_matrix(calls_for_chr_for_sample, refs_for_chr_for_sample)
        intersection_matrix <- filter_intersection_matrix_by_overlap_factor(intersection_matrix, parameters$min_overlap_factor)
        num_of_original_targets_in_refs <- 10 #nrow(targets_in_refs[!duplicated(targets_in_refs[,c("chr", "st_bp", "ed_bp")]),])
        confusion_matrix <- calc_confusion_matrix(intersection_matrix, num_of_original_targets_in_refs, num_of_original_samples_in_refs)
        TP_for_FP <- TP_for_FP + confusion_matrix$TP
        FP <- FP + confusion_matrix$FP
        TN <- TN + confusion_matrix$TN
        FN <- FN + confusion_matrix$FN
      }

      if (nrow(calls_for_chr_for_TP_for_sample) == 0 && nrow(refs_for_chr_for_TP_for_sample) == 0) {  # TODO
        # next()
      } else {
        intersection_matrix <- build_intersection_matrix(calls_for_chr_for_TP_for_sample, refs_for_chr_for_TP_for_sample)
        intersection_matrix <- filter_intersection_matrix_by_overlap_factor(intersection_matrix, parameters$min_overlap_factor)
        num_of_original_targets_in_refs <- 10 #nrow(targets_in_refs[!duplicated(targets_in_refs[,c("chr", "st_bp", "ed_bp")]),])
        confusion_matrix <- calc_confusion_matrix(intersection_matrix, num_of_original_targets_in_refs, num_of_original_samples_in_refs)
        TP <- TP + (nrow(refs_for_chr_for_TP_for_sample) - confusion_matrix$FN)
      }
    }
    confusion_matrix$TP <- TP
    confusion_matrix$FP <- nrow(calls_for_chr) - TP_for_FP
    confusion_matrix$TN <- TN
    confusion_matrix$FN <- FN
    print(confusion_matrix)
    print(calc_quality_statistics(TP, FP, TN, FN))

    TP <- TP + confusion_matrix$TP
    FP <- FP + confusion_matrix$FP
    TN <- TN + confusion_matrix$TN
    FN <- FN + confusion_matrix$FN
  }
  quality_statistics <- calc_quality_statistics(TP, FP, TN, FN)
  print(quality_statistics)

  return(list(TP=TP,
              FP=FP,
              TN=TN,
              FN=FN,
              sensitivity=round(quality_statistics$sensitivity, digits=3), 
              specificity=round(quality_statistics$specificity, digits=3), 
              precision=round(quality_statistics$precision, digits=3), 
              accuracy=round(quality_statistics$accuracy, digits=3)))
}


#######################################################################
####################### calls or refs matrix ##########################
#######################################################################
# "sample_name","chr","cnv","st_bp","ed_bp","copy_no"
# NA0123         1   del   10000   20000   1
