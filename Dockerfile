# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest


## create directories
RUN mkdir -p src
## copy files
COPY src/install_packages.R src/install_packages.R

## install R-packages
RUN Rscript src/install_packages.R

#### IT HAS TO BE FINISHED ####
