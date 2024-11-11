#!/bin/bash
set -e

SHINY_SERVER_VERSION=${1:-${SHINY_SERVER_VERSION:-latest}}

## build ARGs
NCPUS=${NCPUS:--1}

# a function to install apt packages only if they are not installed
function apt_install() {
    if ! dpkg -s "$@" >/dev/null 2>&1; then
        if [ "$(find /var/lib/apt/lists/* | wc -l)" = "0" ]; then
            apt-get update
        fi
        apt-get install -y --allow-downgrades --no-install-recommends "$@"
    fi
}

apt_install \
    sudo \
    gdebi-core \
    libcairo2=1.18.0-1+b1 \
	libcairo-script-interpreter2=1.18.0-1+b1 \
    lsb-release \
    libcurl4-openssl-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget \
    default-jdk \
    libxml2-dev \
    libssl-dev \
    libudunits2-dev \
    librsvg2-dev \
    libmagick++-dev \
    r-cran-ncdf4 \
    libz-dev \
    libnss-winbind \
    winbind \
    dirmngr \
    gnupg \
    apt-transport-https \
    ca-certificates \
    software-properties-common \
    libfftw3-dev \
    texlive \
    libgdal-dev \
    ghostscript \
    g++
    

#Rscript install_step1.R
#Rscript install_step2.R

# Run dependency scripts
/rocker_scripts/install_s6init.sh
/rocker_scripts/install_pandoc.sh

# Install Shiny server

if [ "$SHINY_SERVER_VERSION" = "latest" ]; then
    SHINY_SERVER_VERSION=$(wget -qO- https://download3.rstudio.org/ubuntu-18.04/x86_64/VERSION)
fi

wget --no-verbose "https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-${SHINY_SERVER_VERSION}-amd64.deb" -O ss-latest.deb
gdebi -n ss-latest.deb
rm ss-latest.deb

# Get R packages
install2.r --error --skipinstalled -n "$NCPUS" shiny rmarkdown
