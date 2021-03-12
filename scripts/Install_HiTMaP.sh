sudo apt-get update
sudo apt-get install libfftw3-dev
sudo apt-get install librsvg2-dev
sudo apt-get install tcl-dev tk-dev
sudo apt-get install r-cran-ncdf4
apt install openjdk-11-jre-headless
apt install openjdk-11-jdk-headless
apt-get install libz-dev
sudo apt install libxml2-dev
sudo apt install libssl-dev
sudo apt install libcurl4-openssl-dev
sudo apt-get install libnss-winbind winbind
sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt-cache policy r-base
apt-get purge r-base
sudo apt-get install r-base-core="4.0.2-1.2004.0"
sudo apt-get install libmagick++-dev
apt-get install libfftw3-dev
sudo apt-get install r-base-dev texlive-full
sudo apt-get install libudunits2-dev
sudo apt-get install libgdal-dev
sudo R CMD javareconf
sudo Rscript -e "HiTMaP::HiTMaP_GUI(wd='~/',port = 8787)"