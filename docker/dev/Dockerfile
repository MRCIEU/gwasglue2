FROM mrcieu/gwasglue2-base:0.1.0

# Useful tools for vscode
RUN R -e "install.packages('languageserver')" && \
    R -e 'devtools::install_github("ManuelHentschel/vscDebugger")'

ENTRYPOINT [ "/bin/bash" ]


