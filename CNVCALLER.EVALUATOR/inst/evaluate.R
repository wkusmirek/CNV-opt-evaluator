#!/usr/bin/env Rscript

library('devtools')
library('CNVCALLER.EVALUATOR')
library('optparse')
library('stringr')

option_list <- list(
  make_option("--calls_table", default="public.wiktor_calls_k_to",
              help="Calls table. [default %default]"),
  make_option("--refs_table", default="public.ref_calls_hapmap3",
              help="Reference cnv table. [default %default]"),
  make_option("--targets_table", default="public.targets",
              help="Targets table. [default %default]"),
  make_option("--seg_dups_table", default="public.ref_seg_dups",
              help="Segmental duplications table. [default %default]"),
  make_option("--scenario_id", default="1",
              help="Scenario id. [default %default]"),
  make_option("--min_lratio", default="0.00",
              help="Minimum value of lratio (only for CODEX). [default %default]"),
  make_option("--max_lratio", default="100000.00",
              help="Maximum value of lratio (only for CODEX). [default %default]"),
  make_option("--min_frequency", default="0.00",
              help="Minimum value of cnv frequency (0.00 - 1.00). [default %default]"),
  make_option("--max_frequency", default="1.00",
              help="Maximum value of cnv frequency (0.00 - 1.00). [default %default]"),
  make_option("--min_length", default="0",
              help="Minimum value of cnv length [bp]. [default %default]"),
  make_option("--max_length", default="1000000000",
              help="Maximum value of cnv length [bp]. [default %default]"),
  make_option("--min_targets", default="0",
              help="Minimum number of targets in considered cnv. [default %default]"),
  make_option("--max_targets", default="1000000000",
              help="Maximum number of targets in considered cnv. [default %default]"),
  make_option("--seg_dups_filter_enable", default="1",
              help="Enable segmental duplication filter (0 or 1). [default %default]"),
  make_option("--cnv_type", default="all",
              help="Considered type of cnvs (all, del or dup). [default %default]"),
  make_option("--min_overlap_factor", default="0",
              help="Minimum overlap factor. [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

read_parameters <- function(tabName, id, conn){
  calls_table <- opt$calls_table
  refs_table <- opt$refs_table
  targets_table <- opt$targets_table
  seg_dups_table <- opt$seg_dups_table
  scenario_id <- opt$scenario_id
  min_lratio <- opt$min_lratio
  max_lratio <- opt$max_lratio
  min_frequency <- opt$min_frequency
  max_frequency <- opt$max_frequency
  min_length <- opt$min_length
  max_length <- opt$max_length
  min_targets <- opt$min_targets
  max_targets <- opt$max_targets
  seg_dups_filter_enable <- opt$seg_dups_filter_enable
  min_overlap_factor <- opt$min_overlap_factor
  cnv_type <- opt$cnv_type
  return(list(calls_table=calls_table, 
              refs_table=refs_table,
              targets_table=targets_table,
              seg_dups_table=seg_dups_table,
              scenario_id=scenario_id,
              min_lratio=min_lratio,
              max_lratio=max_lratio,
              min_frequency=min_frequency,
              max_frequency=max_frequency,
              min_length=min_length,
              max_length=max_length,
              min_targets=min_targets,
              max_targets=max_targets,
              seg_dups_filter_enable=seg_dups_filter_enable,
              cnv_type=cnv_type,
              min_overlap_factor=min_overlap_factor))
}

read_cnv_table <- function(tabName, conn, scenario_id, min_lratio, max_lratio){
  cnvs <- read.csv(tabName)
  # map names of columns: from cnvs table to CNVCALLER.EVALUATOR package
  colnames(cnvs)[colnames(cnvs) == 'sample_name'] <- 'sample_name'
  colnames(cnvs)[colnames(cnvs) == 'chr'] <- 'chr'
  colnames(cnvs)[colnames(cnvs) == 'cnv'] <- 'cnv'
  colnames(cnvs)[colnames(cnvs) == 'st_bp'] <- 'st_bp'
  colnames(cnvs)[colnames(cnvs) == 'ed_bp'] <- 'ed_bp'
  colnames(cnvs)[colnames(cnvs) == 'copy_no'] <- 'copy_no'
  cnvs[,'chr'] <- as.character(cnvs[,'chr'])
  cnvs
}

run_evaluator <- function(calls, refs, targets, seg_dups, parameters){
  statistics <- run_CNVCALLER.EVALUATOR(calls,
                                        refs,
                                        targets,
                                        seg_dups,
                                        parameters)
  statistics
}

parameters <- read_parameters(opt$paramsTabName, opt$id, conn_psql)
calls <- read_cnv_table(parameters$calls_table, conn_psql, opt$scenario_id, opt$min_lratio, opt$max_lratio)
print(calls[1:5,])
refs <- read_cnv_table(parameters$refs_table, conn_psql, opt$scenario_id, opt$min_lratio, opt$max_lratio)
print(refs[1:5,])
targets <- read.delim(parameters$targets_table)
print(targets[1:5,])
seg_dups <- read.delim(parameters$seg_dups_table)
print(seg_dups[1:5,])
statistics <- run_evaluator(calls, refs, targets, seg_dups, parameters)
print(statistics)

