#install the git package
install.packages("remotes")
#library(devtools)
library(remotes)
install_github("MASHUOA/HiTMaP",auth_token ="cf6d877f8b6ada1865987b13f9e9996c1883014a",force=T)
3
no
#Update all dependencies
BiocManager::install(ask = F)
yes
library(HiTMaP)