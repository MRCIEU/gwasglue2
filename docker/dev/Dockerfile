
# FROM bioconductor/bioconductor_docker:devel
FROM eddelbuettel/r2u:jammy

# Update this date to re-run apt-get.
# We need to install git with jammy
LABEL package.date=2023-02-03
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
     sqlite3 libcurl4-openssl-dev \
     git


WORKDIR /project

# Get reference datasets
RUN wget http://fileserve.mrcieu.ac.uk/vcf/annotations.vcf.gz.rsidx && \
    wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz && \
    mkdir -p data/reference/ld && \
    mv annotations.vcf.gz.rsidx data/reference/ && \
    mv 1kg.v3.tgz data/reference/ld/ && \
    cd data/reference/ld/ && tar xzvf 1kg.v3.tgz && rm 1kg.v3.tgz && cd -

# Get bcftools
ENV BCFTOOLS_VERSION 1.16
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    make && \
    mv bcftools /bin/ && \
    cd ../ && \
    rm -r bcftools-${BCFTOOLS_VERSION} bcftools-${BCFTOOLS_VERSION}.tar.bz2

# Get plink2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20220814.zip && \
    unzip plink2_linux_x86_64_20220814.zip && \
    mv plink2 /bin/ && \
    rm plink2_linux_x86_64_20220814.zip

# Get plink1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip && \
    mkdir plinktemp && mv plink_linux_x86_64_20220402.zip plinktemp && cd plinktemp && \
    unzip plink_linux_x86_64_20220402.zip && \
    mv plink /bin/ && \
    cd ../ && rm -r plinktemp


# R environment
ENV RENV_VERSION 0.16.0

# Install R packages for development
RUN R -e "install.packages(c('devtools','methods'))"

# Install R packages
RUN R -e "install.packages(c('tidyverse', 'remotes', 'renv', 'coloc', 'susieR', \
        'rmarkdown', 'knitr', 'commonmark', 'markdown','Rfast')); \
        remotes::install_github('mrcieu/ieugwasr'); \
        remotes::install_github('mrcieu/TwoSampleMR'); \
        remotes::install_github('mrcieu/gwasvcf'); \
        remotes::install_github('cnfoley/hyprcoloc', build_opts = c('--resave-data', '--no-manual'), build_vignettes = FALSE)"
# hypercoloc needs pandoc to build vignettes


# Useful tools for vscode
RUN R -e "install.packages('languageserver')" && \
    R -e 'devtools::install_github("ManuelHentschel/vscDebugger")'

# # Install Dentist
# # RUN wget -O DENTIST_1.2.0.0.gz https://www.dropbox.com/s/cq5mmfonocvdiwh/DENTIST.1.2.0.0.gz?dl=0



# Install gwasglue2
# RUN git clone https://github.com/MRCIEU/gwasglue2

ENTRYPOINT [ "/bin/bash" ]


