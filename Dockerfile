# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

MAINTAINER "Bartosz Czech" bartosz.czech@upwr.edu.pl

## create directories
RUN mkdir -p src
## copy files
COPY src/install_packages.R src/install_packages.R
COPY src/drawFusions.R drawFusions.R
COPY src src
## install R-packages
RUN Rscript src/install_packages.R
CMD ["Rscript", "src/install_packages.R"]

RUN chown root:staff drawFusions.R \
	&& chmod 777 drawFusions.R
