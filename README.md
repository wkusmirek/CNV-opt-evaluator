Usage example:

docker run --rm -it -v /home/vboxsmyrna/CNV:/home/vboxsmyrna/CNV -w /home/vboxsmyrna/CNV wkusmirek/cnv-opt-evaluator bash -c "Rscript /usr/local/lib/R/site-library/CNVCALLER.EVALUATOR/evaluate.R --targets_table=/home/vboxsmyrna/CNV/ref/ref_exons.tsv --seg_dups_table=/home/vboxsmyrna/CNV/ref/ref_seg_dups.tsv--calls_table=/home/vboxsmyrna/CNV/codexcov_861_samples_ref_by_kmeans_10/calls_1.csv --refs_table=/home/vboxsmyrna/CNV/ref/ref_calls_1000_861_samples_chr1.csv --seg_dups_filter_enable=0 --max_frequency=0.05 --max_targets=2"

