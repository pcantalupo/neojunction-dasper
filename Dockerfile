# Docker inheritance
FROM bioconductor/bioconductor_docker:RELEASE_3_12

# Install required Bioconductor package
RUN R -e 'BiocManager::install("dasper")'
RUN R -e 'install.packages(c("tidyverse"))'
RUN R -e 'install.packages(c("optparse"))'


