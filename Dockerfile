# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

RUN apt-get update && \
	apt-get install -y libcurl4-openssl-dev libxml2-dev

## copy files
COPY . /code

## install R-packages
RUN Rscript /code/src/install_packages.R

