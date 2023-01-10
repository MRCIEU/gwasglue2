FROM bioconductor/bioconductor_docker:devel

# RUN sudo apt clean
# RUN sudo apt update
# RUN sudo apt install sqlite3

WORKDIR /project
RUN git clone https://github.com/MRCIEU/gwasglue2 -b finemap_coloc

ENTRYPOINT [ "/bin/bash" ]

#This is a test
