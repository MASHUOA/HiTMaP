docker pull mashuoa/hitmap
docker tag mashuoa/hitmap:1.02 mashuoa/hitmap:latest

docker ps
docker commit cf8d5717bbbf mashuoa/hitmap:latest
docker image ls
docker image push mashuoa/hitmap:latest
docker build --no-cache -t mashuoa/hitmap:base1 -f C:\Users\admgguo484\Documents\GitHub\HiTMaP\dockerfiles_shiny_hitmap
docker build --no-cache -t mashuoa/hitmap:base1 .
docker run --rm -p 80:3838 -v %userprofile%\Documents\expdata:/root/expdata mashuoa/hitmap

ssh gguo484@md-cer00340.its.auckland.ac.nz
password:2authentoken

docker run --name hitmap -v %userprofile%\Documents\expdata:/root/expdata -a stdin -a stdout -i -t mashuoa/hitmap:base_nodata /bin/bash 
docker commit cf8d5717bbbf mashuoa/hitmap:latest
docker image push mashuoa/hitmap:latest
#sudo apt-get install vim
#vim etc/ImageMagick-6/policy.xml
#change the pixel limit of width and height
#vi etc/ImageMagick-6/policy.xml
#cp etc/ImageMagick-6/policy.xml /root/exp/
#cp /root/expdata/policy.xml etc/ImageMagick-6/ 
docker tag mashuoa/hitmap:base_R405 mashuoa/hitmap:natcomms

docker tag mashuoa/hitmap:base_R405 mashuoa/hitmap:latest

docker image push mashuoa/hitmap:latest

docker image push mashuoa/hitmap:natcomms

docker tag mashuoa/hitmap:latest mashuoa/hitmap:largefile

docker image push mashuoa/hitmap:largefile

sudo docker stop hitmap
sudo docker rm hitmap

sudo docker pull mashuoa/hitmap:shiny_server

sudo iptables -I INPUT 1 -p tcp --dport 80 -j ACCEPT

sudo docker stop hitmap
sudo docker rm hitmap
sudo docker run --name hitmap -p 80:3838 -v ~/expdata:/root/expdata -a stdin -a stdout -i -t mashuoa/hitmap:gui_latest_deploy /bin/bash 
R
library(HiTMaP)
HiTMaP:::HiTMaP_GUI()


docker run --rm -p 80:3838 -v %userprofile%\Documents\expdata:/root/expdata mashuoa/hitmap


docker run --entrypoint "/bin/bash" --rm  -a stdin -a stdout -i -t rocker/r-ver:4.3.1
wsl --list -v
wsl --export docker-desktop-data "G:\docker-desktop-data.tar"
wsl --unregister docker-desktop-data
wsl --import docker-desktop-data "G:\Docker\data" "G:\docker-desktop-data.tar" --version 2


docker run --name alphafold -a stdin -a stdout -i -t mashuoa/alpha_fold:latest_run /bin/bash

docker run --name hitmapstudio -e PASSWORD=rstudio -p 8787:8787 -v G:\Documents\:/home/rstudio -a stdin -a stdout -i -t mashuoa/hitmap:studio
docker run --name hitmap_docker -p 80:3838 -v %userprofile%\Documents\expdata:/root/expdata -a stdin -a stdout -i -t mashuoa/hitmap:gui_latest /bin/bash 

docker run --entrypoint "/bin/bash" --user root -v G:\Documents\GitHub\HiTMaP\scripts\Docker_file:/go -a stdin -a stdout -i -t python
pip3 install spython
cd /go
spython recipe Dockerfile_base_latest &> Singularity_base.def

docker pull quay.io/singularity/singularity:v4.0.1
docker run --entrypoint "/bin/bash" --name singularity_cmd -v G:\Documents\GitHub\HiTMaP\scripts\Docker_file:/go -a stdin -a stdout -i -t quay.io/singularity/singularity:v4.0.1
sudo singularity build Singularity_base.sif Singularity_base.def

docker run --entrypoint ["/bin/bash"] --name singularity_cmd -v %userprofile%\Documents\go:/go -a stdin -a stdout -i -t quay.io/singularity/singularity:v4.0.1 /bin/bash
docker run --name singularity -v %userprofile%\Documents\go:/go -a stdin -a stdout -i -t quay.io/singularity/singularity:v4.0.1 pull docker://mashuoa/hitmap:latest
docker run -v %userprofile%\Documents\go:/go -a stdin -a stdout -i -t quay.io/singularity/singularity:v4.0.1 shell hitmap_latest.sif





