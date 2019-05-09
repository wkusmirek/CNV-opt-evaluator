Usage example:

docker run --rm -it -v /home/vboxsmyrna/CNV:/home/vboxsmyrna/CNV -w /home/vboxsmyrna/CNV wkusmirek/cnv-opt-evaluator bash -c "Rscript CNVCALLER.EVALUATOR/inst/evaluate_cnvcaller_without_table_data_from_files.R --calls_table=cnvkitcov_2525_samples_ref_by_kmeans_3/calls_1.csv --refs_table=ref/ref_calls_1000_2525_samples_chr1.csv --seg_dups_filter_enable=0 --max_frequency=0.05 --max_targets=2"

