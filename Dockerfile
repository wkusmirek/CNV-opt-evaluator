FROM ubuntu:18.04
MAINTAINER kusmirekwiktor@gmail.com

RUN Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"

RUN Rscript -e "library(devtools);install_github('wkusmirek/CNV-opt-evaluator/CNVCALLER.EVALUATOR', build_vignettes = FALSE)"
