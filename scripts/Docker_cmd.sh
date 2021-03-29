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

#sudo apt-get install vim
#vim etc/ImageMagick-6/policy.xml
#change the pixel limit of width and height
#vi etc/ImageMagick-6/policy.xml
#cp etc/ImageMagick-6/policy.xml /root/exp/
#cp /root/expdata/policy.xml etc/ImageMagick-6/ 