FROM ubuntu:18.04
MAINTAINER kusmirekwiktor@gmail.com

RUN apt-get update && apt-get install -y gnupg2

RUN apt-get install -y software-properties-common

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN apt update

RUN apt-get install -y tzdata

RUN apt-get install -y r-base

RUN apt-get -y install libcurl4-openssl-dev libssl-dev

RUN Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"

RUN Rscript -e "install.packages('optparse', repos = 'http://cran.us.r-project.org')"

RUN Rscript -e "install.packages('stringr', repos = 'http://cran.us.r-project.org')"

RUN Rscript -e "library(devtools);install_github('wkusmirek/CNV-opt-evaluator/CNVCALLER.EVALUATOR', build_vignettes = FALSE)"
