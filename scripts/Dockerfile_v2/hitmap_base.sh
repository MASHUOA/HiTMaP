#!/bin/bash
set -e
function apt_install() {
    if ! dpkg -s "$@" >/dev/null 2>&1; then
        if [ "$(find /var/lib/apt/lists/* | wc -l)" = "0" ]; then
            apt-get update
        fi
        apt-get install -y --allow-downgrades --no-install-recommends "$@"
    fi
}
apt-get update
apt_install \
    sudo \
    gdebi-core \
    libcairo2 \
	libcairo-script-interpreter2 \
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
    libharfbuzz-dev \
    libfribidi-dev \
	libtool \
	automake \
	autoconf \
	m4 \
	perl \
	clang

