
# FROM bioconductor/bioconductor_docker:devel
FROM eddelbuettel/r2u:jammy

# Update this date to re-run apt-get.
# We need to install git with jammy
LABEL package.date = 2023-07-25
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
     sqlite3 libcurl4-openssl-dev \
     git \
     pandoc \
     qpdf -y \
     vim


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
ENV RENV_VERSION 1.0.0

# Install R packages for development
RUN R -e "install.packages(c('devtools','methods','covr'))"

# Install bioconductor packages (need to be installed before gwasvcf)
RUN R -e "install.packages(c('BiocGenerics', 'Biostrings', 'GenomeInfoDb', 'GenomicRanges', \
        'IRanges', ' Rsamtools', 'S4Vectors', 'VariantAnnotation', 'SummarizedExperiment','rtracklayer'))"


# Install R packages
RUN R -e "install.packages(c('tidyverse', 'remotes', 'renv', 'coloc', \
        'rmarkdown', 'knitr', 'bookdown', 'commonmark', 'markdown','Rfast', \
        'httpgd','crayon','hash', 'data.table', 'stringr', ' RCurl', 'rlang', \
        ' utils', 'RSQLite', 'digest','graphics', 'stats','pander')); \       
        remotes::install_github('mrcieu/ieugwasr'); \
        remotes::install_github('mrcieu/TwoSampleMR');  \
        remotes::install_github('mrcieu/gwasvcf'); \
        remotes::install_github('mrcieu/GwasDataImport'); \
        remotes::install_github('ritarasteiro/susieR'); \
        remotes::install_github('ritarasteiro/hyprcoloc', build_opts = c('--resave-data', '--no-manual'), build_vignettes = FALSE)"
# hypercoloc needs pandoc to build vignettes


# Useful tools for vscode
RUN R -e "install.packages('languageserver')" && \
    R -e 'devtools::install_github("ManuelHentschel/vscDebugger")'

# # Install Dentist
# # RUN wget -O DENTIST_1.2.0.0.gz https://www.dropbox.com/s/cq5mmfonocvdiwh/DENTIST.1.2.0.0.gz?dl=0



# Install gwasglue2
# RUN git clone https://github.com/MRCIEU/gwasglue2

ENTRYPOINT [ "/bin/bash" ]


